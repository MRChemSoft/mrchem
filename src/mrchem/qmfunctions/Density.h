#pragma once

#include "constants.h"
#include "parallel.h"

template<int D> class FunctionTree;

class Density {
public:
    Density(bool spin, bool shared);
    Density(bool s = false);
    Density(const Density &rho);
    virtual ~Density();
    void clear();

    void setIsSpinDensity(bool s) { this->is_spin = s; }
    bool isSpinDensity() const { return this->is_spin; }
    int getNNodes(int type = Density::Total) const;
    void setIsShared(bool s) { this->is_shared = s; }
    bool isShared() const { return this->is_shared; }

    void setDensity(int s, FunctionTree<3> *rho);

    void allocTotal();
    void allocSpin();
    void allocAlpha();
    void allocBeta();

    bool hasTotal() const { if (this->dens_t == 0) return false; return true; }
    bool hasSpin() const { if (this->dens_s == 0) return false; return true; }
    bool hasAlpha() const { if (this->dens_a == 0) return false; return true; }
    bool hasBeta() const { if (this->dens_b == 0) return false; return true; }

    bool IsShared() const { if (this->is_shared == 0) return false; return true; }
    SharedMemory* shMem;

    FunctionTree<3> &total() { return *this->dens_t; }
    FunctionTree<3> &spin() { return *this->dens_s; }
    FunctionTree<3> &alpha() { return *this->dens_a; }
    FunctionTree<3> &beta() { return *this->dens_b; }

    const FunctionTree<3> &total() const { return *this->dens_t; }
    const FunctionTree<3> &spin() const { return *this->dens_s; }
    const FunctionTree<3> &alpha() const { return *this->dens_a; }
    const FunctionTree<3> &beta() const { return *this->dens_b; }

    void send_Density(int dest, int tag);
    void Rcv_Density(int source, int tag);
#ifdef HAVE_MPI
    MPI_Win MPI_win = MPI_WIN_NULL;
#endif
    enum Spin { Total, Spin, Alpha, Beta };

protected:
    bool is_spin;
    bool is_shared;
    FunctionTree<3> *dens_t;
    FunctionTree<3> *dens_s;
    FunctionTree<3> *dens_a;
    FunctionTree<3> *dens_b;
};


