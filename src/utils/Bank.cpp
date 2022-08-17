/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "Bank.h"
#include "qmfunctions/Orbital.h"

namespace mrchem {

using namespace Eigen;
using namespace std;

int metadata_block[3]; // can add more metadata in future
int const size_metadata = 3;

Bank::~Bank() {
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
std::map<int, std::map<int, Blockdata_struct *> *> get_nodeid2block; // to get block from its nodeid (all coeff for one node)
std::map<int, std::map<int, Blockdata_struct *> *> get_orbid2block;  // to get block from its orbid

void Bank::open() {
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
    int next_task = 0;
    int tot_ntasks = 0;
    std::map<int, std::vector<int>> readytasks;
    // The bank never goes out of this loop until it receives a close message!
    while (true) {
        MPI_Recv(messages, message_size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_bank, &status);
        if (printinfo) std::cout << world_rank << " got message " << messages[0] << " from " << status.MPI_SOURCE << " account " << messages[1] << " last account " << max_account_id << std::endl;
        int message = messages[0];

        // can be called directly:
        if (message == CLOSE_BANK) {
            if (is_bank and printinfo) std::cout << "Bank is closing" << std::endl;
            this->clear_bank();
            break; // close bank, i.e stop listening for incoming messages
        } else if (message == GET_MAXTOTDATA) {
            int maxsize_int = maxsize / 1024; // convert into MB
            MPI_Send(&maxsize_int, 1, MPI_INT, status.MPI_SOURCE, 1171, comm_bank);
            continue;
        } else if (message == GET_TOTDATA) {
            int maxsize_int = totcurrentsize / 1024; // convert into MB
            MPI_Send(&maxsize_int, 1, MPI_INT, status.MPI_SOURCE, 1172, comm_bank);
            continue;
        } else if (message == NEW_ACCOUNT) {
            // we just have to pick out a number that is not already assigned
            int account = (max_account_id + 1) % 1000000000;
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
            get_readytasks[account] = new std::map<int, std::vector<int>>;
            currentsize[account] = 0;
            MPI_Send(&account, 1, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
            continue;
        }

        // the following is only accessible through an account

        int account = messages[1];
        auto it_dep = get_deposits.find(account);
        if (it_dep == get_deposits.end() || it_dep->second == nullptr) {
            cout << "ERROR, my dep account does not exist!! " << account << " " << message << endl;
            MSG_ABORT("Account error");
        }
        std::vector<deposit> &deposits = *get_deposits[account];
        std::map<int, int> &id2ix = *get_id2ix[account]; // gives zero if id is not defined
        std::map<int, int> &id2qu = *get_id2qu[account];
        std::vector<queue_struct> &queue = *get_queue[account];
        std::map<int, Blockdata_struct *> &orbid2block = *get_orbid2block[account];
        std::map<int, Blockdata_struct *> &nodeid2block = *get_nodeid2block[account];
        auto it_tasks = get_readytasks.find(account);
        if (it_tasks == get_readytasks.end() || it_tasks->second == nullptr) {
            cout << "ERROR, my account does not exist!! " << account << " " << message << endl;
            MSG_ABORT("Account error");
        }
        std::map<int, std::vector<int>> &readytasks = *get_readytasks[account];

        if (message == CLOSE_ACCOUNT) {
            get_numberofclients[account]--;
            if (get_numberofclients[account] == 0) {
                // all clients have closed the account. We remove the account.
                totcurrentsize -= currentsize[account];
                remove_account(account);
            }
        }

        else if (message == CLEAR_BANK) {
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
        } else if (message == CLEAR_BLOCKS) {
            // clear only blocks whith id less than idmax.
            int idmax = messages[2];
            std::vector<int> toeraseVec; // it is dangerous to erase an iterator within its own loop
            for (auto const &block : nodeid2block) {
                if (block.second == nullptr) toeraseVec.push_back(block.first);
                if (block.second == nullptr) continue;
                if (block.first >= idmax and idmax != 0) continue;
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
                //assert(currentsize[account] >= 0);
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
                    if (block.second->id[i] < idmax or idmax == 0) datatoeraseVec.push_back(i);
                    if (block.second->id[i] < idmax or idmax == 0) block.second->data[i] = nullptr;
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

            if (idmax == 0) orbid2block.clear();
            // could have own clear for data?
            for (int ix = 1; ix < deposits.size(); ix++) {
                if (deposits[ix].hasdata) delete deposits[ix].data;
                if (deposits[ix].hasdata) id2ix[deposits[ix].id] = 0; // indicate that it does not exist
                deposits[ix].hasdata = false;
            }
            // send message that it is ready (value of message is not used)
            MPI_Ssend(&message, 1, MPI_INT, status.MPI_SOURCE, 78, comm_bank);
        }

        else if (message == GET_NODEDATA or message == GET_NODEBLOCK) {
            // NB: has no queue system yet
            int nodeid = messages[2]; // which block to fetch from
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
                    if (block->data.size() == 0) std::cout << "Zero size blockdata! " << nodeid << " " << block->N_rows.size() << std::endl;
                    block->BlockData.resize(block->N_rows[0], block->data.size());
                    size = block->N_rows[0] * block->data.size();
                    if (printinfo) std::cout << " rewrite into superblock " << block->data.size() << " " << block->N_rows[0] << " nodeid " << nodeid << std::endl;
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
                    metadata_block[0] = nodeid;             // nodeid
                    metadata_block[1] = block->data.size(); // number of columns
                    metadata_block[2] = size;               // total size = rows*columns
                    MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
                    // send info about the id of each column
                    MPI_Send(block->id.data(), metadata_block[1], MPI_INT, status.MPI_SOURCE, 2, comm_bank);
                }
                double *data_p = block->data[dataindex];
                if (size > 0) MPI_Send(data_p, size, MPI_DOUBLE, status.MPI_SOURCE, 3, comm_bank);
            } else {
                if (printinfo) std::cout << " block " << nodeid << " does not exist " << std::endl;
                // Block with this id does not exist.
                if (message == GET_NODEDATA) {
                    int size = messages[4]; // number of doubles to send
                    if (size == 0) {
                        std::cout << "WARNING: GET_NODEDATA asks for zero size data" << std::endl;
                        metadata_block[2] = size;
                        MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, 3, comm_bank);
                    } else {
                        std::vector<double> zero(size, 0.0); // send zeroes
                        MPI_Ssend(zero.data(), size, MPI_DOUBLE, status.MPI_SOURCE, 3, comm_bank);
                    }
                } else {
                    metadata_block[0] = nodeid;
                    metadata_block[1] = 0; // number of columns
                    metadata_block[2] = 0; // total size = rows*columns
                    MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, 3, comm_bank);
                }
            }
        } else if (message == GET_ORBBLOCK) {
            // NB: BLOCKDATA has no queue system yet
            int orbid = messages[2]; // which block to fetch from

            if (orbid2block.count(orbid) and orbid2block[orbid] != nullptr) {
                Blockdata_struct *block = orbid2block[orbid];
                int dataindex = 0; // internal index of the data in the block
                int size = 0;
                // send entire block. First make one contiguous superblock
                // Prepare the data as one contiguous block
                if (block->data.size() == 0) std::cout << "Zero size blockdata! C " << orbid << " " << block->N_rows.size() << std::endl;
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
                MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
                MPI_Send(block->id.data(), metadata_block[1], MPI_INT, status.MPI_SOURCE, 2, comm_bank);
                MPI_Send(coeff.data(), size, MPI_DOUBLE, status.MPI_SOURCE, 3, comm_bank);
            } else {
                // it is possible and allowed that the block has not been written
                if (printinfo) std::cout << " block does not exist " << orbid << " " << orbid2block.count(orbid) << std::endl;
                // Block with this id does not exist.
                metadata_block[0] = orbid;
                metadata_block[1] = 0; // number of columns
                metadata_block[2] = 0; // total size = rows*columns
                MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
            }
        }

        else if (message == GET_ORBITAL or message == GET_ORBITAL_AND_WAIT or message == GET_ORBITAL_AND_DELETE or message == GET_FUNCTION or message == GET_DATA) {
            // withdrawal
            int id = messages[2];
            int ix = id2ix[id];
            if (id2ix.count(id) == 0 or ix == 0) {
                if (printinfo) std::cout << world_rank << " not found " << id << " " << message << std::endl;
                if (message == GET_ORBITAL or message == GET_ORBITAL_AND_DELETE) {
                    // do not wait for the orbital to arrive
                    int found = 0;
                    if (printinfo) std::cout << world_rank << " sending found 0 to " << status.MPI_SOURCE << std::endl;
                    MPI_Send(&found, 1, MPI_INT, status.MPI_SOURCE, 117, comm_bank);
                } else {
                    // the id does not exist. Put in queue and Wait until it is defined
                    if (printinfo) std::cout << world_rank << " queuing " << id << std::endl;
                    if (id2qu[id] == 0) {
                        queue.push_back({id, {status.MPI_SOURCE}});
                        id2qu[id] = queue.size() - 1;
                    } else {
                        // somebody is already waiting for this id. queue in queue
                        queue[id2qu[id]].clients.push_back(status.MPI_SOURCE);
                    }
                }
            } else {
                int ix = id2ix[id];
                if (deposits[ix].id != id) std::cout << ix << " Bank accounting error " << std::endl;
                if (message == GET_ORBITAL or message == GET_ORBITAL_AND_WAIT or message == GET_ORBITAL_AND_DELETE) {
                    if (message == GET_ORBITAL or message == GET_ORBITAL_AND_DELETE) {
                        int found = 1;
                        MPI_Send(&found, 1, MPI_INT, status.MPI_SOURCE, 117, comm_bank);
                    }
                    send_orbital(*deposits[ix].orb, status.MPI_SOURCE, 1, comm_bank);
                    if (message == GET_ORBITAL_AND_DELETE) {
                        currentsize[account] -= deposits[ix].orb->getSizeNodes(NUMBER::Total);
                        totcurrentsize -= deposits[ix].orb->getSizeNodes(NUMBER::Total);
                        deposits[ix].orb->free(NUMBER::Total);
                        id2ix[id] = 0;
                    }
                }
                if (message == GET_FUNCTION) { send_function(*deposits[ix].orb, status.MPI_SOURCE, 1, comm_bank); }
                if (message == GET_DATA) { MPI_Send(deposits[ix].data, deposits[ix].datasize, MPI_DOUBLE, status.MPI_SOURCE, 1, comm_bank); }
            }
        } else if (message == SAVE_NODEDATA) {
            int nodeid = messages[2]; // which block to write
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

            MPI_Recv(data_p, size, MPI_DOUBLE, status.MPI_SOURCE, 1, comm_bank, &status);
            if (printinfo) std::cout << " written block " << nodeid << " id " << orbid << " subblocks " << nodeid2block[nodeid]->data.size() << std::endl;
        } else if (message == SAVE_ORBITAL or message == SAVE_FUNCTION or message == SAVE_DATA) {
            // make a new deposit
            int exist_flag = 0;
            int id = messages[2];
            if (id2ix[id]) {
                std::cout << "WARNING: id " << id << " exists already"
                          << " " << status.MPI_SOURCE << " " << message << " " << std::endl;
                ix = id2ix[id]; // the deposit exist from before. Will be overwritten
                exist_flag = 1;
                if (message == SAVE_DATA and !deposits[ix].hasdata) {
                    datasize = messages[3];
                    exist_flag = 0;
                    deposits[ix].data = new double[datasize];
                    deposits[ix].hasdata = true;
                }
            } else {
                ix = deposits.size(); // NB: ix is now index of last element + 1
                deposits.resize(ix + 1);
                if (message == SAVE_ORBITAL or message == SAVE_FUNCTION) deposits[ix].orb = new Orbital(0);
                if (message == SAVE_DATA) {
                    datasize = messages[3];
                    deposits[ix].data = new double[datasize];
                    deposits[ix].hasdata = true;
                }
            }
            deposits[ix].id = id;
            id2ix[deposits[ix].id] = ix;
            deposits[ix].source = status.MPI_SOURCE;
            if (message == SAVE_ORBITAL) {
                recv_orbital(*deposits[ix].orb, deposits[ix].source, 1, comm_bank);
                if (exist_flag == 0) {
                    currentsize[account] += deposits[ix].orb->getSizeNodes(NUMBER::Total);
                    totcurrentsize += deposits[ix].orb->getSizeNodes(NUMBER::Total);
                    this->maxsize = std::max(totcurrentsize, this->maxsize);
                }
            }
            if (message == SAVE_FUNCTION) { recv_function(*deposits[ix].orb, deposits[ix].source, 1, comm_bank); }
            if (message == SAVE_DATA) {
                datasize = messages[3];
                deposits[ix].datasize = datasize;
                MPI_Recv(deposits[ix].data, datasize, MPI_DOUBLE, deposits[ix].source, 1, comm_bank, &status);
                currentsize[account] += datasize / 128; // converted into kB
                totcurrentsize += datasize / 128;       // converted into kB
                this->maxsize = std::max(totcurrentsize, this->maxsize);
            }
            if (id2qu[deposits[ix].id] != 0) {
                // someone is waiting for those data. Send to them
                int iq = id2qu[deposits[ix].id];
                if (deposits[ix].id != queue[iq].id) std::cout << ix << " Bank queue accounting error " << std::endl;
                for (int iqq : queue[iq].clients) {
                    if (message == SAVE_ORBITAL) { send_orbital(*deposits[ix].orb, iqq, 1, comm_bank); }
                    if (message == SAVE_FUNCTION) { send_function(*deposits[ix].orb, iqq, 1, comm_bank); }
                    if (message == SAVE_DATA) { MPI_Send(deposits[ix].data, messages[3], MPI_DOUBLE, iqq, 1, comm_bank); }
                }
                queue[iq].clients.clear(); // cannot erase entire queue[iq], because that would require to shift all the
                                           // id2qu value larger than iq
                queue[iq].id = -1;
                id2qu.erase(deposits[ix].id);
            }

            // Task manager members:
        } else if (message == INIT_TASKS) {
            tot_ntasks = messages[2];
            next_task = 0;
        } else if (message == GET_NEXTTASK) {
            int task = next_task;
            if (next_task >= tot_ntasks) task = -1; // flag to show all tasks are assigned
            MPI_Send(&task, 1, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
            next_task++;
        } else if (message == PUT_READYTASK) {
            readytasks[messages[2]].push_back(messages[3]);
        }
        if (message == DEL_READYTASK) {
            for (int i = 0; i < readytasks[messages[2]].size(); i++) { // we expect small sizes
                if (readytasks[messages[2]][i] == messages[3]) {
                    readytasks[messages[2]].erase(readytasks[messages[2]].begin() + i);
                    break;
                }
            }
        } else if (message == GET_READYTASK) {
            int nready = 0;
            if (readytasks.count(messages[2]) > 0) nready = readytasks[messages[2]].size();
            MPI_Send(&nready, 1, MPI_INT, status.MPI_SOURCE, 844, mpi::comm_bank);
            if (nready > 0) MPI_Send(readytasks[messages[2]].data(), nready, MPI_INT, status.MPI_SOURCE, 845, mpi::comm_bank);
        } else if (message == GET_READYTASK_DEL) {
            int nready = 0;
            if (readytasks.count(messages[2]) > 0) { nready = readytasks[messages[2]].size(); }
            MPI_Send(&nready, 1, MPI_INT, status.MPI_SOURCE, 844, mpi::comm_bank);
            if (nready > 0) MPI_Send(readytasks[messages[2]].data(), nready, MPI_INT, status.MPI_SOURCE, 845, mpi::comm_bank);
            if (nready > 0) readytasks[messages[2]].resize(0);
        }
    }
#endif
}

// Ask to close the Bank
void Bank::close() {
#ifdef MRCHEM_HAS_MPI
    int messages[message_size];
    messages[0] = CLOSE_BANK;
    for (int i = 0; i < bank_size; i++) MPI_Send(messages, 1, MPI_INT, bankmaster[i], 0, comm_bank);
    if (tot_bank_size > bank_size) MPI_Send(messages, 1, MPI_INT, task_bank, 0, comm_bank);
#endif
}

void Bank::clear_bank() {
#ifdef MRCHEM_HAS_MPI
    for (auto account : accounts) { remove_account(account); }
#endif
}

int Bank::clearAccount(int account, int iclient, MPI_Comm comm) {
#ifdef MRCHEM_HAS_MPI
    closeAccount(account);
    return openAccount(iclient, comm);
#else
    return 1;
#endif
}
void Bank::remove_account(int account) {
#ifdef MRCHEM_HAS_MPI
    auto it = get_deposits.find(account);
    if (it == get_deposits.end() || it->second == nullptr) {
        cout << "ERROR, my account depositsdoes not exist!! " << account << " " << endl;
        MSG_ABORT("depositsAccount error");
    }
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
    delete get_readytasks[account];
    get_id2ix.erase(account);
    get_id2qu.erase(account);
    get_deposits.erase(account);
    currentsize.erase(account);
    get_readytasks.erase(account);

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
        //assert(currentsize[account] >= 0);
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
        if (deposits[ix].hasdata) delete deposits[ix].data;
        if (deposits[ix].hasdata) (*get_id2ix[account])[deposits[ix].id] = 0; // indicate that it does not exist
        deposits[ix].hasdata = false;
    }
    delete get_nodeid2block[account];
    delete get_orbid2block[account];
    get_nodeid2block.erase(account);
    get_orbid2block.erase(account);

#endif
}

int Bank::openAccount(int iclient, MPI_Comm comm) {
    // NB: this is a collective call, since we need all the accounts to be synchronized
    int account_id[1] = {-1};
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = NEW_ACCOUNT;
    messages[1] = 0;
    int size;
    MPI_Comm_size(comm, &size);
    messages[1] = size;
    if (iclient == 0) {
        for (int i = 0; i < bank_size; i++) {
            int account_id_i[1];
            MPI_Send(messages, message_size, MPI_INT, bankmaster[i], 0, comm_bank);
            MPI_Recv(account_id_i, 1, MPI_INT, bankmaster[i], 1, comm_bank, &status);
            if (i > 0 and account_id_i[0] != account_id[0]) MSG_ABORT("Account id mismatch!");
            account_id[0] = account_id_i[0];
        }
        MPI_Bcast(account_id, 1, MPI_INT, 0, comm);
    } else {
        MPI_Bcast(account_id, 1, MPI_INT, 0, comm);
    }
#endif
    return account_id[0];
}

int Bank::openTaskManager(int ntasks, int iclient, MPI_Comm comm) {
    // NB: this is a collective call, since we need all the accounts to be synchronized
    int account_id = -1;
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = NEW_ACCOUNT;
    int size;
    MPI_Comm_size(comm, &size);
    messages[1] = size;
    if (iclient == 0) {
        MPI_Send(messages, 2, MPI_INT, task_bank, 0, comm_bank);
        MPI_Recv(&account_id, 1, MPI_INT, task_bank, 1, comm_bank, &status);
        if (tot_bank_size == bank_size) {
            // make a dummy account so that all account_id are synchronized
            int account_id_i;
            for (int i = 0; i < bank_size; i++) {
                if (bankmaster[i] != task_bank) {
                    MPI_Send(messages, 2, MPI_INT, bankmaster[i], 0, comm_bank);
                    MPI_Recv(&account_id_i, 1, MPI_INT, bankmaster[i], 1, comm_bank, &status);
                    if (i > 0 and account_id_i != account_id) MSG_ABORT("Account id mismatch!");
                }
            }
        }
        messages[0] = INIT_TASKS;
        messages[1] = account_id;
        messages[2] = ntasks;
        MPI_Send(messages, 3, MPI_INT, task_bank, 2, comm_bank);
        MPI_Bcast(&account_id, 1, MPI_INT, 0, comm);
    } else {
        MPI_Bcast(&account_id, 1, MPI_INT, 0, comm);
    }
#endif
    return account_id;
}

void Bank::closeAccount(int account_id) {
// The account will in reality not be removed before everybody has sent a close message
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = CLOSE_ACCOUNT;
    messages[1] = account_id;
    for (int i = 0; i < bank_size; i++) MPI_Send(messages, 2, MPI_INT, bankmaster[i], 0, comm_bank);
#endif
}

void Bank::closeTaskManager(int account_id) {
// The account will in reality not be removed before everybody has sent a close message
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = CLOSE_ACCOUNT;
    messages[1] = account_id;
    MPI_Send(messages, 2, MPI_INT, task_bank, 0, comm_bank);
#endif
}

int Bank::get_maxtotalsize() {
    int maxtot = 0;
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    int datasize;
    int messages[message_size];
    messages[0] = GET_MAXTOTDATA;
    for (int i = 0; i < bank_size; i++) {
        MPI_Send(messages, 1, MPI_INT, bankmaster[i], 0, comm_bank);
        MPI_Recv(&datasize, 1, MPI_INT, bankmaster[i], 1171, comm_bank, &status);
        maxtot = std::max(maxtot, datasize);
    }
#endif
    return maxtot;
}

std::vector<int> Bank::get_totalsize() {
    std::vector<int> tot;
#ifdef HAVE_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = GET_TOTDATA;
    int datasize;
    for (int i = 0; i < bank_size; i++) {
        MPI_Send(messages, 1, MPI_INT, bankmaster[i], 0, comm_bank);
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
    int messages[message_size];
    messages[0] = SAVE_ORBITAL;
    messages[1] = account_id;
    messages[2] = id;
    MPI_Send(messages, 3, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    send_orbital(orb, bankmaster[id % bank_size], 1, comm_bank);
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
    messages[2] = id;
    if (wait == 0) {
        messages[0] = GET_ORBITAL;
        MPI_Send(messages, 3, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
        int found;
        MPI_Recv(&found, 1, MPI_INT, bankmaster[id % bank_size], 117, comm_bank, &status);
        if (found != 0) {
            recv_orbital(orb, bankmaster[id % bank_size], 1, comm_bank);
            return 1;
        } else {
            return 0;
        }
    } else {
        messages[0] = GET_ORBITAL_AND_WAIT;
        MPI_Send(messages, 3, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
        recv_orbital(orb, bankmaster[id % bank_size], 1, comm_bank);
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
    messages[2] = id;
    MPI_Send(messages, 3, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    int found;
    MPI_Recv(&found, 1, MPI_INT, bankmaster[id % bank_size], 117, comm_bank, &status);
    if (found != 0) {
        recv_orbital(orb, bankmaster[id % bank_size], 1, comm_bank);
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
    messages[2] = id;
    MPI_Send(messages, 3, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    send_function(func, bankmaster[id % bank_size], 1, comm_bank);
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
    messages[2] = id;
    MPI_Send(messages, 3, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    recv_function(func, bankmaster[id % bank_size], 1, comm_bank);
#endif
    return 1;
}

// save data in Bank with identity id . datasize MUST have been set already. NB:not tested
int BankAccount::put_data(int id, int size, double *data) {
#ifdef MRCHEM_HAS_MPI
    // for now we distribute according to id
    if (id > max_tag) MSG_ABORT("Bank id must be less than max allowed tag");
    int messages[message_size];
    messages[0] = SAVE_DATA;
    messages[1] = account_id;
    messages[2] = id;
    messages[3] = size;
    MPI_Send(messages, 4, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    MPI_Send(data, size, MPI_DOUBLE, bankmaster[id % bank_size], 1, comm_bank);
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
    messages[2] = id;
    MPI_Send(messages, 3, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    MPI_Recv(data, size, MPI_DOUBLE, bankmaster[id % bank_size], 1, comm_bank, &status);
#endif
    return 1;
}

// save data in Bank with identity id as part of block with identity nodeid.
int BankAccount::put_nodedata(int id, int nodeid, int size, double *data) {
#ifdef MRCHEM_HAS_MPI
    // for now we distribute according to nodeid
    if (id > max_tag) MSG_ABORT("Bank id must be less than max allowed tag");
    int messages[message_size];
    messages[0] = SAVE_NODEDATA;
    messages[1] = account_id;
    messages[2] = nodeid; // which block
    messages[3] = id;     // id within block
    messages[4] = size;   // size of this data
    MPI_Send(messages, 5, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Send(data, size, MPI_DOUBLE, bankmaster[nodeid % bank_size], 1, comm_bank);
#endif
    return 1;
}

// get data with identity id
int BankAccount::get_nodedata(int id, int nodeid, int size, double *data, std::vector<int> &idVec) {
#ifdef MRCHEM_HAS_MPI
    MPI_Status status;
    // get the column with identity id
    int messages[message_size];
    messages[0] = GET_NODEDATA;
    messages[1] = account_id;
    messages[2] = nodeid; // which block
    messages[3] = id;     // id within block.
    messages[4] = size;   // expected size of data
    MPI_Send(messages, 5, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Recv(data, size, MPI_DOUBLE, bankmaster[nodeid % bank_size], 3, comm_bank, &status);
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
    messages[2] = nodeid;

    MPI_Send(messages, 3, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Recv(metadata_block, size_metadata, MPI_INT, bankmaster[nodeid % bank_size], 1, comm_bank, &status);
    idVec.resize(metadata_block[1]);
    int size = metadata_block[2];
    if (size > 0) MPI_Recv(idVec.data(), metadata_block[1], MPI_INT, bankmaster[nodeid % bank_size], 2, comm_bank, &status);
    if (size > 0) MPI_Recv(data, size, MPI_DOUBLE, bankmaster[nodeid % bank_size], 3, comm_bank, &status);
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
    messages[2] = orbid;
    MPI_Send(messages, 3, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Recv(metadata_block, size_metadata, MPI_INT, bankmaster[nodeid % bank_size], 1, comm_bank, &status);
    nodeidVec.resize(metadata_block[1]);
    int totsize = metadata_block[2];
    if (totsize > 0) MPI_Recv(nodeidVec.data(), metadata_block[1], MPI_INT, bankmaster[nodeid % bank_size], 2, comm_bank, &status);
    data = new double[totsize];
    if (totsize > 0) MPI_Recv(data, totsize, MPI_DOUBLE, bankmaster[nodeid % bank_size], 3, comm_bank, &status);
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
        messages[0] = CLEAR_BLOCKS;
        messages[1] = account_id;
        messages[2] = nodeidmax;
        for (int i = 0; i < bank_size; i++) { MPI_Send(messages, 3, MPI_INT, bankmaster[i], 0, comm_bank); }
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

// creator. NB: collective
TaskManager::TaskManager(int ntasks, int iclient, MPI_Comm comm) {
    this->n_tasks = ntasks;
    if (bank_size == 0) return;
    this->account_id = dataBank.openTaskManager(ntasks, iclient, comm);
#ifdef MRCHEM_HAS_MPI
    MPI_Barrier(comm);
#endif
}

// destructor
TaskManager::~TaskManager() {
    // The account will in reality not be removed before everybody has sent a delete message
    if (this->account_id < 0) return;
    dataBank.closeTaskManager(this->account_id);
}

int TaskManager::next_task() {
    int nexttask = 0;
#ifdef MRCHEM_HAS_MPI
    if (this->account_id >= 0) {
        MPI_Status status;
        int messages[message_size];
        messages[0] = GET_NEXTTASK;
        messages[1] = account_id;
        MPI_Send(messages, message_size, MPI_INT, task_bank, 0, comm_bank);
        MPI_Recv(&nexttask, 1, MPI_INT, task_bank, 1, comm_bank, &status);
        return nexttask;
    }
#endif
    nexttask = this->task++;
    if (nexttask >= this->n_tasks) nexttask = -1;
    return nexttask;
}

void TaskManager::put_readytask(int i, int j) {
#ifdef MRCHEM_HAS_MPI
    if (this->account_id < 0) return;
    MPI_Status status;
    int messages[message_size];
    messages[0] = PUT_READYTASK;
    messages[1] = account_id;
    messages[2] = i;
    messages[3] = j;
    MPI_Send(messages, message_size, MPI_INT, task_bank, 0, comm_bank);
#endif
}

void TaskManager::del_readytask(int i, int j) {
#ifdef MRCHEM_HAS_MPI
    if (this->account_id < 0) return;
    MPI_Status status;
    int messages[message_size];
    messages[0] = DEL_READYTASK;
    messages[1] = account_id;
    messages[2] = i;
    messages[3] = j;
    MPI_Send(messages, message_size, MPI_INT, task_bank, 0, comm_bank);
#endif
}

std::vector<int> TaskManager::get_readytask(int i, int del) {
    std::vector<int> readytasks;
#ifdef MRCHEM_HAS_MPI
    if (this->account_id < 0) return readytasks;
    MPI_Status status;
    int messages[message_size];
    messages[0] = GET_READYTASK;
    if (del == 1) messages[0] = GET_READYTASK_DEL;
    messages[1] = account_id;
    messages[2] = i;
    MPI_Send(messages, message_size, MPI_INT, task_bank, 0, comm_bank);
    int nready;
    MPI_Recv(&nready, 1, MPI_INT, task_bank, 844, comm_bank, &status);
    if (nready > 0) {
        readytasks.resize(nready);
        MPI_Recv(readytasks.data(), nready, MPI_INT, task_bank, 845, comm_bank, &status);
    }
#endif
    return readytasks;
}

} // namespace mrchem
