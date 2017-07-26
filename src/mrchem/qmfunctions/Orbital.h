#ifndef ORBITAL_H
#define ORBITAL_H

#include <complex>

#include "constants.h"

#include "QMFunction.h"

class Orbital : public QMFunction {
public:
    Orbital(int occ, int s);
    Orbital(const Orbital &orb);
    Orbital &operator=(const Orbital &orb) { NOT_IMPLEMENTED_ABORT;}
    virtual ~Orbital() { clear(); }
    void clear(bool free = true);

    int getSpin() const { return this->spin; }
    int getOccupancy() const { return this->occupancy; }
    double getError() const { return this->error; }

    void setSpin(int s) { this->spin = s; }
    void setOccupancy(int occ) { this->occupancy = occ; }
    void setError(double err) { this->error = err; }

    bool isConverged(double prec) const;

    void compare(const Orbital &orb) const;
    int compareSpin(const Orbital &orb) const;
    int compareOccupancy(const Orbital &orb) const;

    std::complex<double> dot(Orbital &ket);
    /*! determines the exchange factor to be used in the calculation of the exact exchange
     *
     * \param [in] orb input orbital to which K is applied
     *
     * The factor is computed in terms of the occupancy of the two orbitals and in terms of the spin
     * 0.5 factors are used in order to preserve occupancy of the set of doubly occupied orbitals
     * this-> is the orbital defining the operator whereas the input orbital (orb) is the one              
     * the operator is applied to
     *
     * Occupancy: Single/Double
     * Spin: alpha/beta
     *
     * K (this->) | orb (input) | factor
     * alpha      | alpha       | 1.0       
     * alpha      | beta        | 0.0       
     * alpha      | double      | 0.5      
     * -------------------------------
     * beta       | alpha       | 0.0       
     * beta       | beta        | 1.0       
     * beta       | double      | 0.5 
     * -------------------------------
     * double     | alpha       | 1.0       
     * double     | beta        | 1.0       
     * double     | double      | 1.0       
     *
     */
    double getExchangeFactor(const Orbital &orb) const;

    char printSpin() const {
        char sp = 'u';
        if (this->spin == Alpha) sp = 'a';
        if (this->spin == Beta) sp = 'b';
        return sp;
    }

    void send_Orbital(int dest, int tag);
    void Rcv_Orbital(int source, int tag);
    void Isend_Orbital(int dest, int tag);
    void IRcv_Orbital(int source, int tag);

    friend std::ostream& operator<<(std::ostream &o, Orbital &orb) {
        o << std::setw(25) << orb.getSquareNorm();
        o << std::setw(3) << orb.getOccupancy();
        o << std::setw(4) << orb.printSpin();
        o << std::setw(24) << orb.getError() << std::endl;
        return o;
    }

    enum Spin { Paired, Alpha, Beta };

protected:
    int spin;
    int occupancy;
    double error;
};

#endif // ORBITAL_H
