#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "Bank.h"
#include "qmfunctions/Orbital.h"

namespace mrchem {

using namespace Eigen;

int metadata_block[3]; // can add more metadata in future
int const size_metadata = 3;

CentralBank::~CentralBank() {
    // delete all data and accounts
}

struct Blockdata_struct {
    std::vector<double *> data; // to store the incoming data
    std::vector<bool> deleted;  // to indicate if it has been deleted already
    MatrixXd BlockData;         // to put the final block
    // eigen matrix are per default stored column-major (one can add columns at the end)
    std::vector<int> N_rows;
    std::map<int, int> id2data; // internal index of the data in the block
    std::vector<int> id;        // the id of each column. Either nodeid, or orbid
};
std::map<int, std::map<int, Blockdata_struct *> *>
    get_nodeid2block; // to get block from its nodeid (all coeff for one node)
std::map<int, std::map<int, Blockdata_struct *> *> get_orbid2block; // to get block from its orbid

void CentralBank::open() {
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    char safe_data1;
    int deposit_size = sizeof(deposit);
    int n_chunks, ix;
    int messages[message_size];
    int datasize = -1;
    std::map<int, int> get_numberofclients;

    bool printinfo = false;
    int id_shift = max_tag / 2; // to ensure that nodes, orbitals and functions do not collide
    int max_account_id = -1;
    // The bank never goes out of this loop until it receives a close message!
    while (true) {
        MPI_Recv(&messages, message_size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_bank, &status);
        if (printinfo)
            std::cout << world_rank << " got message " << messages[0] << " from " << status.MPI_SOURCE << " account "
                      << messages[1] << std::endl;
        int message = messages[0];
        int account = messages[1];
        if (message == NEW_ACCOUNT) {
            // we just have to pick out a number that is not already assigned
            account = (max_account_id + 1) % 1000000000;
            while (get_deposits.count(account)) account = (account + 1) % 1000000000; // improbable this is used
            max_account_id = account;
            // create default content
            get_deposits[account] = new std::vector<deposit>;
            get_deposits[account]->resize(1);
            get_id2ix[account] = new std::map<int, int>;
            get_id2qu[account] = new std::map<int, int>;
            get_queue[account] = new std::vector<queue_struct>;
            get_queue[account]->resize(1);
            get_orbid2block[account] = new std::map<int, Blockdata_struct *>;
            get_nodeid2block[account] = new std::map<int, Blockdata_struct *>;
            get_numberofclients[account] = messages[1];
            currentsize[account] = 0;
            MPI_Send(&account, 1, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
        }
        std::vector<deposit> &deposits = *get_deposits[account];
        std::map<int, int> &id2ix = *get_id2ix[account]; // gives zero if id is not defined
        std::map<int, int> &id2qu = *get_id2qu[account];
        std::vector<queue_struct> &queue = *get_queue[account];
        std::map<int, Blockdata_struct *> &orbid2block = *get_orbid2block[account];
        std::map<int, Blockdata_struct *> &nodeid2block = *get_nodeid2block[account];

        if (message == CLOSE_ACCOUNT) {
            get_numberofclients[account]--;
            if (get_numberofclients[account] == 0) {
                // all clients have closed the account. We remove the account.
                totcurrentsize -= currentsize[account];
                clear_account(account);
            }
        }

        if (message == CLOSE_BANK) {
            if (is_bank and printinfo) std::cout << "Bank is closing" << std::endl;
            this->clear_bank();
            break; // close bank, i.e stop listening for incoming messages
        }
        if (message == CLEAR_BANK) {
            this->clear_bank();
            for (auto const &block : nodeid2block) {
                if (block.second == nullptr) continue;
                for (int i = 0; i < block.second->data.size(); i++) {
                    if (not block.second->deleted[i]) {
                        currentsize[account] -= block.second->N_rows[i] / 128; // converted into kB
                        totcurrentsize -= block.second->N_rows[i] / 128;       // converted into kB
                        delete[] block.second->data[i];
                    }
                }
                delete block.second;
            }
            nodeid2block.clear();
            orbid2block.clear();
            // send message that it is ready (value of message is not used)
            MPI_Ssend(&message, 1, MPI_INT, status.MPI_SOURCE, 77, comm_bank);
        }
        if (message == CLEAR_BLOCKS) {
            // clear only blocks whith id less than status.MPI_TAG.
            std::vector<int> toeraseVec; // it is dangerous to erase an iterator within its own loop
            for (auto const &block : nodeid2block) {
                if (block.second == nullptr) toeraseVec.push_back(block.first);
                if (block.second == nullptr) continue;
                if (block.first >= status.MPI_TAG and status.MPI_TAG != 0) continue;
                for (int i = 0; i < block.second->data.size(); i++) {
                    if (not block.second->deleted[i]) {
                        currentsize[account] -= block.second->N_rows[i] / 128; // converted into kB
                        totcurrentsize -= block.second->N_rows[i] / 128;       // converted into kB
                        delete[] block.second->data[i];
                    }
                }
                currentsize[account] -= block.second->BlockData.size() / 128; // converted into kB
                totcurrentsize -= block.second->BlockData.size() / 128;       // converted into kB
                block.second->BlockData.resize(0, 0); // NB: the matrix does not clear itself otherwise
                assert(currentsize[account] >= 0);
                this->currentsize[account] = std::max(0ll, currentsize[account]);
                toeraseVec.push_back(block.first);
            }
            for (int ierase : toeraseVec) { nodeid2block.erase(ierase); }
            toeraseVec.clear();
            std::vector<int> datatoeraseVec;
            for (auto const &block : orbid2block) {
                if (block.second == nullptr) toeraseVec.push_back(block.first);
                if (block.second == nullptr) continue;
                datatoeraseVec.clear();
                for (int i = 0; i < block.second->data.size(); i++) {
                    if (block.second->id[i] < status.MPI_TAG or status.MPI_TAG == 0) datatoeraseVec.push_back(i);
                    if (block.second->id[i] < status.MPI_TAG or status.MPI_TAG == 0) block.second->data[i] = nullptr;
                }
                std::sort(datatoeraseVec.begin(), datatoeraseVec.end());
                std::reverse(datatoeraseVec.begin(), datatoeraseVec.end());
                for (int ierase : datatoeraseVec) {
                    block.second->id.erase(block.second->id.begin() + ierase);
                    block.second->data.erase(block.second->data.begin() + ierase);
                    block.second->N_rows.erase(block.second->N_rows.begin() + ierase);
                }
                if (block.second->data.size() == 0) toeraseVec.push_back(block.first);
            }
            for (int ierase : toeraseVec) { orbid2block.erase(ierase); }

            if (status.MPI_TAG == 0) orbid2block.clear();
            // could have own clear for data?
            for (int ix = 1; ix < deposits.size(); ix++) {
                if (deposits[ix].id >= id_shift) {
                    if (deposits[ix].hasdata) delete deposits[ix].data;
                    if (deposits[ix].hasdata) id2ix[deposits[ix].id] = 0; // indicate that it does not exist
                    deposits[ix].hasdata = false;
                }
            }
            // send message that it is ready (value of message is not used)
            MPI_Ssend(&message, 1, MPI_INT, status.MPI_SOURCE, 78, comm_bank);
        }
        if (message == GETMAXTOTDATA) {
            int maxsize_int = maxsize / 1024; // convert into MB
            MPI_Send(&maxsize_int, 1, MPI_INT, status.MPI_SOURCE, 1171, comm_bank);
        }
        if (message == GETTOTDATA) {
            int maxsize_int = totcurrentsize / 1024; // convert into MB
            MPI_Send(&maxsize_int, 1, MPI_INT, status.MPI_SOURCE, 1172, comm_bank);
        }

        if (message == GET_NODEDATA or message == GET_NODEBLOCK) {
            // NB: has no queue system yet
            int nodeid = status.MPI_TAG; // which block to fetch from
            if (nodeid2block.count(nodeid) and nodeid2block[nodeid] != nullptr) {
                Blockdata_struct *block = nodeid2block[nodeid];
                int dataindex = 0; // internal index of the data in the block
                int size = 0;
                if (message == GET_NODEDATA) {
                    int orbid = messages[3];           // which part of the block to fetch
                    dataindex = block->id2data[orbid]; // column of the data in the block
                    size = block->N_rows[dataindex];   // number of doubles to fetch
                    if (size != messages[4]) std::cout << "ERROR nodedata has wrong size" << std::endl;
                } else {
                    // send entire block. First make one contiguous superblock
                    // Prepare the data as one contiguous block
                    if (block->data.size() == 0)
                        std::cout << "Zero size blockdata! " << nodeid << " " << block->N_rows.size() << std::endl;
                    block->BlockData.resize(block->N_rows[0], block->data.size());
                    size = block->N_rows[0] * block->data.size();
                    if (printinfo)
                        std::cout << " rewrite into superblock " << block->data.size() << " " << block->N_rows[0]
                                  << " tag " << status.MPI_TAG << std::endl;
                    for (int j = 0; j < block->data.size(); j++) {
                        for (int i = 0; i < block->N_rows[j]; i++) { block->BlockData(i, j) = block->data[j][i]; }
                    }
                    // repoint to the data in BlockData
                    for (int j = 0; j < block->data.size(); j++) {
                        if (block->deleted[j] == true) std::cout << "ERROR data already deleted " << std::endl;
                        assert(block->deleted[j] == false);
                        delete[] block->data[j];
                        block->deleted[j] = true;
                        block->data[j] = block->BlockData.col(j).data();
                    }
                    dataindex = 0; // start from first column
                    // send info about the size of the superblock
                    metadata_block[0] = status.MPI_TAG;     // nodeid
                    metadata_block[1] = block->data.size(); // number of columns
                    metadata_block[2] = size;               // total size = rows*columns
                    MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, nodeid, comm_bank);
                    // send info about the id of each column
                    MPI_Send(block->id.data(), metadata_block[1], MPI_INT, status.MPI_SOURCE, nodeid + 1, comm_bank);
                }
                double *data_p = block->data[dataindex];
                if (size > 0) MPI_Send(data_p, size, MPI_DOUBLE, status.MPI_SOURCE, nodeid + 2, comm_bank);
            } else {
                if (printinfo) std::cout << " block " << nodeid << " does not exist " << std::endl;
                // Block with this id does not exist.
                if (message == GET_NODEDATA) {
                    int size = messages[4]; // number of doubles to send
                    if (size == 0) {
                        std::cout << "WARNING: GET_NODEDATA asks for zero size data" << std::endl;
                        metadata_block[2] = size;
                        MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, nodeid, comm_bank);
                    } else {
                        std::vector<double> zero(size, 0.0); // send zeroes
                        MPI_Ssend(zero.data(), size, MPI_DOUBLE, status.MPI_SOURCE, nodeid + 2, comm_bank);
                    }
                } else {
                    metadata_block[0] = status.MPI_TAG; // nodeid
                    metadata_block[1] = 0;              // number of columns
                    metadata_block[2] = 0;              // total size = rows*columns
                    MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, metadata_block[0], comm_bank);
                }
            }
        }
        if (message == GET_ORBBLOCK) {
            // NB: BLOCKDATA has no queue system yet
            int orbid = status.MPI_TAG; // which block to fetch from

            if (orbid2block.count(orbid) and orbid2block[orbid] != nullptr) {
                Blockdata_struct *block = orbid2block[orbid];
                int dataindex = 0; // internal index of the data in the block
                int size = 0;
                // send entire block. First make one contiguous superblock
                // Prepare the data as one contiguous block
                if (block->data.size() == 0)
                    std::cout << "Zero size blockdata! C " << orbid << " " << block->N_rows.size() << std::endl;
                size = 0;
                for (int j = 0; j < block->data.size(); j++) size += block->N_rows[j];

                std::vector<double> coeff(size);
                int ij = 0;
                for (int j = 0; j < block->data.size(); j++) {
                    for (int i = 0; i < block->N_rows[j]; i++) { coeff[ij++] = block->data[j][i]; }
                }
                // send info about the size of the superblock
                metadata_block[0] = orbid;
                metadata_block[1] = block->data.size(); // number of columns
                metadata_block[2] = size;               // total size = rows*columns
                MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, orbid, comm_bank);
                MPI_Send(block->id.data(), metadata_block[1], MPI_INT, status.MPI_SOURCE, orbid + 1, comm_bank);
                MPI_Send(coeff.data(), size, MPI_DOUBLE, status.MPI_SOURCE, orbid + 2, comm_bank);
            } else {
                // it is possible and allowed that the block has not been written
                if (printinfo)
                    std::cout << " block does not exist " << orbid << " " << orbid2block.count(orbid) << std::endl;
                // Block with this id does not exist.
                metadata_block[0] = orbid;
                metadata_block[1] = 0; // number of columns
                metadata_block[2] = 0; // total size = rows*columns
                MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, orbid, comm_bank);
            }
        }

        if (message == GET_ORBITAL or message == GET_ORBITAL_AND_WAIT or message == GET_ORBITAL_AND_DELETE or
            message == GET_FUNCTION or message == GET_DATA) {
            // withdrawal
            int ix = id2ix[status.MPI_TAG];
            if (id2ix.count(status.MPI_TAG) == 0 or ix == 0) {
                if (printinfo)
                    std::cout << world_rank << " not found " << status.MPI_TAG << " " << message << std::endl;
                if (message == GET_ORBITAL or message == GET_ORBITAL_AND_DELETE) {
                    // do not wait for the orbital to arrive
                    int found = 0;
                    if (printinfo) std::cout << world_rank << " sending found 0 to " << status.MPI_SOURCE << std::endl;
                    MPI_Send(&found, 1, MPI_INT, status.MPI_SOURCE, 117, comm_bank);
                } else {
                    // the id does not exist. Put in queue and Wait until it is defined
                    if (printinfo) std::cout << world_rank << " queuing " << status.MPI_TAG << std::endl;
                    if (id2qu[status.MPI_TAG] == 0) {
                        queue.push_back({status.MPI_TAG, {status.MPI_SOURCE}});
                        id2qu[status.MPI_TAG] = queue.size() - 1;
                    } else {
                        // somebody is already waiting for this id. queue in queue
                        queue[id2qu[status.MPI_TAG]].clients.push_back(status.MPI_SOURCE);
                    }
                }
            } else {
                int ix = id2ix[status.MPI_TAG];
                if (deposits[ix].id != status.MPI_TAG) std::cout << ix << " Bank accounting error " << std::endl;
                if (message == GET_ORBITAL or message == GET_ORBITAL_AND_WAIT or message == GET_ORBITAL_AND_DELETE) {
                    if (message == GET_ORBITAL or message == GET_ORBITAL_AND_DELETE) {
                        int found = 1;
                        MPI_Send(&found, 1, MPI_INT, status.MPI_SOURCE, 117, comm_bank);
                    }
                    send_orbital(*deposits[ix].orb, status.MPI_SOURCE, deposits[ix].id, comm_bank);
                    if (message == GET_ORBITAL_AND_DELETE) {
                        currentsize[account] -= deposits[ix].orb->getSizeNodes(NUMBER::Total);
                        totcurrentsize -= deposits[ix].orb->getSizeNodes(NUMBER::Total);
                        deposits[ix].orb->free(NUMBER::Total);
                        id2ix[status.MPI_TAG] = 0;
                    }
                }
                if (message == GET_FUNCTION) {
                    send_function(*deposits[ix].orb, status.MPI_SOURCE, deposits[ix].id, comm_bank);
                }
                if (message == GET_DATA) {
                    MPI_Send(deposits[ix].data,
                             deposits[ix].datasize,
                             MPI_DOUBLE,
                             status.MPI_SOURCE,
                             deposits[ix].id,
                             comm_bank);
                }
            }
        }
        if (message == SAVE_NODEDATA) {
            int nodeid = messages[2]; // which block to write (should = status.MPI_TAG)
            int orbid = messages[3];  // which part of the block
            int size = messages[4];   // number of doubles

            // test if the block exists already
            if (printinfo) std::cout << world_rank << " save data nodeid " << nodeid << " size " << size << std::endl;
            if (nodeid2block.count(nodeid) == 0 or nodeid2block[nodeid] == nullptr) {
                if (printinfo) std::cout << world_rank << " block does not exist yet  " << std::endl;
                // the block does not exist yet, create it
                Blockdata_struct *block = new Blockdata_struct;
                nodeid2block[nodeid] = block;
            }
            if (orbid2block.count(orbid) == 0 or orbid2block[orbid] == nullptr) {
                // the block does not exist yet, create it
                Blockdata_struct *orbblock = new Blockdata_struct;
                orbid2block[orbid] = orbblock;
            }
            // append the incoming data
            Blockdata_struct *block = nodeid2block[nodeid];
            block->id2data[orbid] = nodeid2block[nodeid]->data.size(); // internal index of the data in the block
            double *data_p = new double[size];
            currentsize[account] += size / 128; // converted into kB
            totcurrentsize += size / 128;       // converted into kB
            this->maxsize = std::max(totcurrentsize, this->maxsize);
            block->data.push_back(data_p);
            block->deleted.push_back(false);
            block->id.push_back(orbid);
            block->N_rows.push_back(size);

            Blockdata_struct *orbblock = orbid2block[orbid];
            orbblock->id2data[nodeid] = orbblock->data.size(); // internal index of the data in the block
            orbblock->data.push_back(data_p);
            orbblock->deleted.push_back(false);
            orbblock->id.push_back(nodeid);
            orbblock->N_rows.push_back(size);

            MPI_Recv(data_p, size, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, comm_bank, &status);
            if (printinfo)
                std::cout << " written block " << nodeid << " id " << orbid << " subblocks "
                          << nodeid2block[nodeid]->data.size() << std::endl;
        }
        if (message == SAVE_ORBITAL or message == SAVE_FUNCTION or message == SAVE_DATA) {
            // make a new deposit
            int exist_flag = 0;
            if (id2ix[status.MPI_TAG]) {
                std::cout << "WARNING: id " << status.MPI_TAG << " exists already"
                          << " " << status.MPI_SOURCE << " " << message << " " << std::endl;
                ix = id2ix[status.MPI_TAG]; // the deposit exist from before. Will be overwritten
                exist_flag = 1;
                if (message == SAVE_DATA and !deposits[ix].hasdata) {
                    exist_flag = 0;
                    deposits[ix].data = new double[datasize];
                    deposits[ix].hasdata = true;
                }
            } else {
                ix = deposits.size(); // NB: ix is now index of last element + 1
                deposits.resize(ix + 1);
                if (message == SAVE_ORBITAL or message == SAVE_FUNCTION) deposits[ix].orb = new Orbital(0);
                if (message == SAVE_DATA) {
                    deposits[ix].data = new double[datasize];
                    deposits[ix].hasdata = true;
                }
            }
            deposits[ix].id = status.MPI_TAG;
            id2ix[deposits[ix].id] = ix;
            deposits[ix].source = status.MPI_SOURCE;
            if (message == SAVE_ORBITAL) {
                recv_orbital(*deposits[ix].orb, deposits[ix].source, deposits[ix].id, comm_bank);
                if (exist_flag == 0) {
                    currentsize[account] += deposits[ix].orb->getSizeNodes(NUMBER::Total);
                    totcurrentsize += deposits[ix].orb->getSizeNodes(NUMBER::Total);
                    this->maxsize = std::max(totcurrentsize, this->maxsize);
                }
            }
            if (message == SAVE_FUNCTION) {
                recv_function(*deposits[ix].orb, deposits[ix].source, deposits[ix].id, comm_bank);
            }
            if (message == SAVE_DATA) {
                deposits[ix].datasize = datasize;
                MPI_Recv(
                    deposits[ix].data, datasize, MPI_DOUBLE, deposits[ix].source, deposits[ix].id, comm_bank, &status);
                currentsize[account] += datasize / 128; // converted into kB
                totcurrentsize += datasize / 128;       // converted into kB
                this->maxsize = std::max(totcurrentsize, this->maxsize);
            }
            if (id2qu[deposits[ix].id] != 0) {
                // someone is waiting for those data. Send to them
                int iq = id2qu[deposits[ix].id];
                if (deposits[ix].id != queue[iq].id) std::cout << ix << " Bank queue accounting error " << std::endl;
                for (int iqq : queue[iq].clients) {
                    if (message == SAVE_ORBITAL) { send_orbital(*deposits[ix].orb, iqq, queue[iq].id, comm_bank); }
                    if (message == SAVE_FUNCTION) { send_function(*deposits[ix].orb, iqq, queue[iq].id, comm_bank); }
                    if (message == SAVE_DATA) {
                        MPI_Send(deposits[ix].data, datasize, MPI_DOUBLE, iqq, queue[iq].id, comm_bank);
                    }
                }
                queue[iq].clients.clear(); // cannot erase entire queue[iq], because that would require to shift all the
                                           // id2qu value larger than iq
                queue[iq].id = -1;
                id2qu.erase(deposits[ix].id);
            }
        }
        if (message == SET_DATASIZE) {
            int datasize_new = messages[2];
            if (datasize_new != datasize) {
                // make sure that all old data arrays are deleted
                for (int ix = 1; ix < deposits.size(); ix++) {
                    if (deposits[ix].hasdata) {
                        delete deposits[ix].data;
                        deposits[ix].hasdata = false;
                    }
                }
            }
            datasize = datasize_new;
        }
    }
#endif
}

// Ask to close the Bank
void CentralBank::close() {
#ifdef MRCHEM_HAS_MPI
    for (int i = 0; i < bank_size; i++) { MPI_Send(&CLOSE_BANK, 1, MPI_INT, bankmaster[i], 0, comm_bank); }
#endif
}

void CentralBank::clear_bank() {
#ifdef MRCHEM_HAS_MPI
    for (auto account : accounts) { clear_account(account); }
#endif
}

int CentralBank::clearAccount(int account, int iclient, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    closeAccount(account);
    return openAccount(iclient, comm);
#endif
}
void CentralBank::clear_account(int account) {
#ifdef MRCHEM_HAS_MPI
    std::vector<deposit> &deposits = *get_deposits[account];
    for (int ix = 1; ix < deposits.size(); ix++) {
        if (deposits[ix].orb != nullptr) deposits[ix].orb->free(NUMBER::Total);
        if (deposits[ix].hasdata) delete deposits[ix].data;
        deposits[ix].hasdata = false;
    }
    deposits.clear();
    delete get_queue[account];
    get_queue.erase(account);
    delete get_id2ix[account];
    delete get_id2qu[account];
    get_id2ix.erase(account);
    get_id2qu.erase(account);
    get_deposits.erase(account);
    currentsize.erase(account);

    std::map<int, Blockdata_struct *> &nodeid2block = *get_nodeid2block[account];
    std::map<int, Blockdata_struct *> &orbid2block = *get_orbid2block[account];

    std::vector<int> toeraseVec; // it is dangerous to erase an iterator within its own loop
    for (auto const &block : nodeid2block) {
        if (block.second == nullptr) toeraseVec.push_back(block.first);
        if (block.second == nullptr) continue;
        for (int i = 0; i < block.second->data.size(); i++) {
            if (not block.second->deleted[i]) {
                currentsize[account] -= block.second->N_rows[i] / 128; // converted into kB
                totcurrentsize -= block.second->N_rows[i] / 128;       // converted into kB
                delete[] block.second->data[i];
            }
        }
        currentsize[account] -= block.second->BlockData.size() / 128; // converted into kB
        totcurrentsize -= block.second->BlockData.size() / 128;       // converted into kB
        block.second->BlockData.resize(0, 0);                         // NB: the matrix does not clear itself otherwise
        assert(currentsize[account] >= 0);
        toeraseVec.push_back(block.first);
    }
    for (int ierase : toeraseVec) { nodeid2block.erase(ierase); }
    toeraseVec.clear();
    std::vector<int> datatoeraseVec; // it is dangerous to erase an iterator within its own loop
    for (auto const &block : orbid2block) {
        if (block.second == nullptr) toeraseVec.push_back(block.first);
        if (block.second == nullptr) continue;
        datatoeraseVec.clear();
        for (int i = 0; i < block.second->data.size(); i++) {
            datatoeraseVec.push_back(i);
            block.second->data[i] = nullptr;
        }
        std::sort(datatoeraseVec.begin(), datatoeraseVec.end());
        std::reverse(datatoeraseVec.begin(), datatoeraseVec.end());
        for (int ierase : datatoeraseVec) {
            block.second->id.erase(block.second->id.begin() + ierase);
            block.second->data.erase(block.second->data.begin() + ierase);
            block.second->N_rows.erase(block.second->N_rows.begin() + ierase);
        }
        if (block.second->data.size() == 0) toeraseVec.push_back(block.first);
    }
    for (int ierase : toeraseVec) { orbid2block.erase(ierase); }

    orbid2block.clear();
    for (int ix = 1; ix < deposits.size(); ix++) {
        if (deposits[ix].id >= max_tag / 2) {
            if (deposits[ix].hasdata) delete deposits[ix].data;
            if (deposits[ix].hasdata) (*get_id2ix[account])[deposits[ix].id] = 0; // indicate that it does not exist
            deposits[ix].hasdata = false;
        }
    }
    delete get_nodeid2block[account];
    delete get_orbid2block[account];
    get_nodeid2block.erase(account);
    get_orbid2block.erase(account);

#endif
}

int CentralBank::openAccount(int iclient, MPI_Comm comm) {
// NB: this is a collective call, since we need all the accounts to be synchronized
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = NEW_ACCOUNT;
    int size;
    MPI_Comm_size(comm, &size);
    messages[1] = size;
    int account_id = -1;
    if (iclient == 0) {
        for (int i = 0; i < bank_size; i++) {
            int account_id_i;
            MPI_Send(messages, message_size, MPI_INT, bankmaster[i], 0, comm_bank);
            MPI_Recv(&account_id_i, 1, MPI_INT, bankmaster[i], 1, comm_bank, &status);
            if (i > 0 and account_id_i != account_id) MSG_ABORT("Account id mismatch!");
            account_id = account_id_i;
        }
        MPI_Bcast(&account_id, 1, MPI_INT, 0, comm);
    } else {
        MPI_Bcast(&account_id, 1, MPI_INT, 0, comm);
    }
    return account_id;
#endif
}

void CentralBank::closeAccount(int account_id) {
// The account will in reality not be removed before everybody has sent a close message
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = CLOSE_ACCOUNT;
    messages[1] = account_id;
    for (int i = 0; i < bank_size; i++) { MPI_Send(messages, message_size, MPI_INT, bankmaster[i], 0, comm_bank); }
#endif
}

int CentralBank::get_maxtotalsize() {
    int maxtot = 0;
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int datasize;
    for (int i = 0; i < bank_size; i++) {
        MPI_Send(&GETMAXTOTDATA, 1, MPI_INT, bankmaster[i], 0, comm_bank);
        MPI_Recv(&datasize, 1, MPI_INT, bankmaster[i], 1171, comm_bank, &status);
        maxtot = std::max(maxtot, datasize);
    }
#endif
    return maxtot;
}

std::vector<int> CentralBank::get_totalsize() {
    std::vector<int> tot;
#ifdef HAVE_MPI
    MPI_Status status;
    int datasize;
    for (int i = 0; i < bank_size; i++) {
        MPI_Send(&GETTOTDATA, 1, MPI_INT, bankmaster[i], 0, comm_bank);
        MPI_Recv(&datasize, 1, MPI_INT, bankmaster[i], 1172, comm_bank, &status);
        tot.push_back(datasize);
    }
#endif
    return tot;
}

// Accounts: (clients)

// save orbital in Bank with identity id
int BankAccount::put_orb(int id, Orbital &orb) {
#ifdef MRCHEM_HAS_MPI
    // for now we distribute according to id
    if (id > max_tag) MSG_ABORT("Bank id must be less than max allowed tag ");
    int messages[message_size];
    messages[0] = SAVE_ORBITAL;
    messages[1] = account_id;
    MPI_Send(messages, message_size, MPI_INT, bankmaster[id % bank_size], id, comm_bank);
    send_orbital(orb, bankmaster[id % bank_size], id, comm_bank);
#endif
    return 1;
}

// get orbital with identity id.
// If wait=0, return immediately with value zero if not available (default)
// else, wait until available
int BankAccount::get_orb(int id, Orbital &orb, int wait) {
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[1] = account_id;
    if (wait == 0) {
        messages[0] = GET_ORBITAL;
        MPI_Send(messages, message_size, MPI_INT, bankmaster[id % bank_size], id, comm_bank);
        int found;
        MPI_Recv(&found, 1, MPI_INT, bankmaster[id % bank_size], 117, comm_bank, &status);
        if (found != 0) {
            recv_orbital(orb, bankmaster[id % bank_size], id, comm_bank);
            return 1;
        } else {
            return 0;
        }
    } else {
        messages[0] = GET_ORBITAL_AND_WAIT;
        MPI_Send(messages, message_size, MPI_INT, bankmaster[id % bank_size], id, comm_bank);
        recv_orbital(orb, bankmaster[id % bank_size], id, comm_bank);
    }
#endif
    return 1;
}

// get orbital with identity id, and delete from bank.
// return immediately with value zero if not available
int BankAccount::get_orb_del(int id, Orbital &orb) {
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = GET_ORBITAL_AND_DELETE;
    messages[1] = account_id;
    MPI_Send(messages, message_size, MPI_INT, bankmaster[id % bank_size], id, comm_bank);
    int found;
    MPI_Recv(&found, 1, MPI_INT, bankmaster[id % bank_size], 117, comm_bank, &status);
    if (found != 0) {
        recv_orbital(orb, bankmaster[id % bank_size], id, comm_bank);
        return 1;
    } else {
        return 0;
    }
#endif
    return 1;
}

// save function in Bank with identity id
int BankAccount::put_func(int id, QMFunction &func) {
#ifdef MRCHEM_HAS_MPI
    // for now we distribute according to id
    if (id > max_tag / 2) MSG_ABORT("Bank id must be less than max allowed tag / 2");
    id += max_tag / 2;
    int messages[message_size];
    messages[0] = SAVE_FUNCTION;
    messages[1] = account_id;
    MPI_Send(messages, message_size, MPI_INT, bankmaster[id % bank_size], id, comm_bank);
    send_function(func, bankmaster[id % bank_size], id, comm_bank);
#endif
    return 1;
}

// get function with identity id
int BankAccount::get_func(int id, QMFunction &func) {
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    id += max_tag / 2;
    int messages[message_size];
    messages[0] = GET_FUNCTION;
    messages[1] = account_id;
    MPI_Send(messages, message_size, MPI_INT, bankmaster[id % bank_size], id, comm_bank);
    recv_function(func, bankmaster[id % bank_size], id, comm_bank);
#endif
    return 1;
}

// set the size of the data arrays (in size of doubles) to be sent/received later
void BankAccount::set_datasize(int datasize, int iclient, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    if (iclient == 0) {
        for (int i = 0; i < bank_size; i++) {
            int messages[message_size];
            messages[0] = SET_DATASIZE;
            messages[1] = account_id;
            messages[2] = datasize;
            MPI_Send(messages, message_size, MPI_INT, bankmaster[i], 0, comm_bank);
        }
    }
    MPI_Barrier(comm);
#endif
}

// save data in Bank with identity id . datasize MUST have been set already. NB:not tested
int BankAccount::put_data(int id, int size, double *data) {
#ifdef MRCHEM_HAS_MPI
    // for now we distribute according to id
    if (id > max_tag) MSG_ABORT("Bank id must be less than max allowed tag");
    int messages[message_size];
    messages[0] = SAVE_DATA;
    messages[1] = account_id;
    MPI_Send(messages, message_size, MPI_INT, bankmaster[id % bank_size], id, comm_bank);
    MPI_Send(data, size, MPI_DOUBLE, bankmaster[id % bank_size], id, comm_bank);
#endif
    return 1;
}

// get data with identity id
int BankAccount::get_data(int id, int size, double *data) {
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = GET_DATA;
    messages[1] = account_id;
    MPI_Send(messages, message_size, MPI_INT, bankmaster[id % bank_size], id, comm_bank);
    MPI_Recv(data, size, MPI_DOUBLE, bankmaster[id % bank_size], id, comm_bank, &status);
#endif
    return 1;
}

// save data in Bank with identity id as part of block with identity nodeid.
int BankAccount::put_nodedata(int id, int nodeid, int size, double *data) {
#ifdef MRCHEM_HAS_MPI
    // for now we distribute according to nodeid
    if (id > max_tag) MSG_ABORT("Bank id must be less than max allowed tag");
    metadata_block[0] = nodeid; // which block
    metadata_block[1] = id;     // id within block
    metadata_block[2] = size;   // size of this data
    int messages[message_size];
    messages[0] = SAVE_NODEDATA;
    messages[1] = account_id;
    messages[2] = nodeid; // which block
    messages[3] = id;     // id within block
    messages[4] = size;   // size of this data
    MPI_Send(messages, message_size, MPI_INT, bankmaster[nodeid % bank_size], nodeid, comm_bank);
    MPI_Send(data, size, MPI_DOUBLE, bankmaster[nodeid % bank_size], nodeid, comm_bank);
#endif
    return 1;
}

// get data with identity id
int BankAccount::get_nodedata(int id, int nodeid, int size, double *data, std::vector<int> &idVec) {
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    // get the column with identity id
    metadata_block[0] = nodeid; // which block
    metadata_block[1] = id;     // id within block.
    metadata_block[2] = size;   // expected size of data
    int messages[message_size];
    messages[0] = GET_NODEDATA;
    messages[1] = account_id;
    messages[2] = nodeid; // which block
    messages[3] = id;     // id within block.
    messages[4] = size;   // expected size of data
    MPI_Send(messages, message_size, MPI_INT, bankmaster[nodeid % bank_size], nodeid, comm_bank);
    MPI_Recv(data, size, MPI_DOUBLE, bankmaster[nodeid % bank_size], nodeid + 2, comm_bank, &status);
#endif
    return 1;
}

// get all data for nodeid
int BankAccount::get_nodeblock(int nodeid, double *data, std::vector<int> &idVec) {
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    // get the entire superblock and also the id of each column
    int messages[message_size];
    messages[0] = GET_NODEBLOCK;
    messages[1] = account_id;
    MPI_Send(messages, message_size, MPI_INT, bankmaster[nodeid % bank_size], nodeid, comm_bank);
    MPI_Recv(metadata_block, size_metadata, MPI_INT, bankmaster[nodeid % bank_size], nodeid, comm_bank, &status);
    idVec.resize(metadata_block[1]);
    int size = metadata_block[2];
    if (size > 0)
        MPI_Recv(
            idVec.data(), metadata_block[1], MPI_INT, bankmaster[nodeid % bank_size], nodeid + 1, comm_bank, &status);
    if (size > 0) MPI_Recv(data, size, MPI_DOUBLE, bankmaster[nodeid % bank_size], nodeid + 2, comm_bank, &status);
#endif
    return 1;
}

// get all data with identity orbid
int BankAccount::get_orbblock(int orbid, double *&data, std::vector<int> &nodeidVec, int bankstart) {
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int nodeid = orb_rank + bankstart;
    // get the entire superblock and also the nodeid of each column
    int messages[message_size];
    messages[0] = GET_ORBBLOCK;
    messages[1] = account_id;
    MPI_Send(messages, message_size, MPI_INT, bankmaster[nodeid % bank_size], orbid, comm_bank);
    MPI_Recv(metadata_block, size_metadata, MPI_INT, bankmaster[nodeid % bank_size], orbid, comm_bank, &status);
    nodeidVec.resize(metadata_block[1]);
    int totsize = metadata_block[2];
    if (totsize > 0)
        MPI_Recv(nodeidVec.data(),
                 metadata_block[1],
                 MPI_INT,
                 bankmaster[nodeid % bank_size],
                 orbid + 1,
                 comm_bank,
                 &status);
    data = new double[totsize];
    if (totsize > 0) MPI_Recv(data, totsize, MPI_DOUBLE, bankmaster[nodeid % bank_size], orbid + 2, comm_bank, &status);
#endif
    return 1;
}

// remove all blockdata with nodeid < nodeidmax
// NB:: collective call. All clients must call this
void BankAccount::clear_blockdata(int iclient, int nodeidmax, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    // 1) wait until all clients are ready
    MPI_Barrier(comm);
    // master send signal to bank
    if (iclient == 0) {
        int messages[message_size];
        messages[1] = account_id;
        messages[0] = CLEAR_BLOCKS;
        for (int i = 0; i < bank_size; i++) {
            MPI_Send(messages, message_size, MPI_INT, bankmaster[i], nodeidmax, comm_bank);
        }
        for (int i = 0; i < bank_size; i++) {
            // wait until Bank is finished and has sent signal
            MPI_Status status;
            int message;
            MPI_Recv(&message, 1, MPI_INT, bankmaster[i], 78, comm_bank, &status);
        }
    }
    MPI_Barrier(comm);
#endif
}

// creator. NB: collective
BankAccount::BankAccount(int iclient, MPI_Comm comm) {
    this->account_id = dataBank.openAccount(iclient, comm);
#ifdef MRCHEM_HAS_MPI
    MPI_Barrier(comm);
#endif
}

// destructor
BankAccount::~BankAccount() {
    // The account will in reality not be removed before everybody has sent a delete message
    dataBank.closeAccount(this->account_id);
}

// closes account and reopen a new empty account. NB: account_id will change
void BankAccount::clear(int iclient, MPI_Comm comm) {
    this->account_id = dataBank.clearAccount(this->account_id, iclient, comm);
}

} // namespace mrchem
