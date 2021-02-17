#pragma once

#include "MRCPP/Parallel"

#include "mrchem.h"
#include "parallel.h"
#include "qmfunctions/qmfunction_fwd.h"

namespace mrchem {

using namespace mpi;

struct deposit {
    Orbital *orb;
    double *data; // for pure data arrays
    bool hasdata;
    int datasize;
    int id = -1; // to identify what is deposited
    int source;  // mpi rank from the source of the data
};

struct queue_struct {
    int id;
    std::vector<int> clients;
};
int const CLOSE_BANK = 1;
int const CLEAR_BANK = 2;
int const NEW_ACCOUNT = 3;
int const CLOSE_ACCOUNT = 4;
int const GET_ORBITAL = 5;
int const GET_ORBITAL_AND_WAIT = 6;
int const GET_ORBITAL_AND_DELETE = 7;
int const SAVE_ORBITAL = 8;
int const GET_FUNCTION = 9;
int const SAVE_FUNCTION = 10;
int const SET_DATASIZE = 11;
int const GET_DATA = 12;
int const SAVE_DATA = 13;
int const SAVE_NODEDATA = 14;
int const GET_NODEDATA = 15;
int const GET_NODEBLOCK = 16;
int const GET_ORBBLOCK = 17;
int const CLEAR_BLOCKS = 18;
int const GETMAXTOTDATA = 19;
int const GETTOTDATA = 20;

class CentralBank {
public:
    CentralBank() = default;
    ~CentralBank();
    void open();
    void close();
    int openAccount(int iclient, MPI_Comm comm);
    int clearAccount(int account, int iclient, MPI_Comm comm); // closes and open fresh account
    void closeAccount(int account_id);                         // remove the account
    long long totcurrentsize = 0ll;                            // number of kB used by all accounts
    int get_maxtotalsize();

private:
    std::vector<int> accounts;                          // open bank accounts
    std::map<int, std::vector<deposit> *> get_deposits; // gives deposits of an account
    std::map<int, std::map<int, int> *> get_id2ix;
    std::map<int, std::map<int, int> *> get_id2qu;
    std::map<int, std::vector<queue_struct> *> get_queue; // gives deposits of an account
    std::map<int, long long> currentsize;                 // total deposited data size (without containers)
    long long maxsize = 0;                                // max total deposited data size (without containers)
    void clear_bank();
    void clear_account(int account); // remove the content of the account
    std::vector<int> get_totalsize();
};

class BankAccount {
public:
    BankAccount(int iclient = orb_rank, MPI_Comm comm = comm_orb);
    ~BankAccount();
    int account_id = -1;
    void clear(int i = orb_rank, MPI_Comm comm = comm_orb);
    int put_orb(int id, Orbital &orb);
    int get_orb(int id, Orbital &orb, int wait = 0);
    int get_orb_del(int id, Orbital &orb);
    int put_func(int id, QMFunction &func);
    int get_func(int id, QMFunction &func);
    void set_datasize(int datasize, int iclient = orb_rank, MPI_Comm comm = comm_orb);
    int put_data(int id, int size, double *data);
    int get_data(int id, int size, double *data);
    int put_nodedata(int id, int nodeid, int size, double *data);
    int get_nodedata(int id, int nodeid, int size, double *data, std::vector<int> &idVec);
    int get_nodeblock(int nodeid, double *data, std::vector<int> &idVec);
    int get_orbblock(int orbid, double *&data, std::vector<int> &nodeidVec, int bankstart);
    void clear_blockdata(int i = orb_rank, int nodeidmax = 0, MPI_Comm comm = comm_orb);
};

int const message_size = 5;

} // namespace mrchem
