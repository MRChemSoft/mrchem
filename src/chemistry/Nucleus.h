#pragma once

#include <vector>
#include <iostream>

#include "PeriodicTable.h"
#include "Element.h"

namespace mrchem {

class Nucleus final {
public:
    Nucleus(const Element &elm, const double *r = 0)
            : charge(elm.getZ()),
              element(&elm) {
        setCoord(r);
    }
    Nucleus(const Nucleus &nuc)
            : charge(nuc.charge),
              element(nuc.element) {
        setCoord(nuc.coord);
    }
    const Nucleus &operator=(const Nucleus &nuc) {
        if (this != &nuc) {
            this->charge = nuc.charge;
            this->element = nuc.element;
            setCoord(nuc.coord);
        }
        return *this;
    }

    ~Nucleus() { }

    void setCharge(double z) { this->charge = z; }
    void setCoord(const double *r) {
        for (int d = 0; d < 3; d++) {
            if (r != 0) {
                this->coord[d] = r[d];
            } else {
                this->coord[d] = 0.0;
            }
        }
    }

    double getCharge() const { return this->charge; }
    const double *getCoord() const { return this->coord; }
    const Element &getElement() const { return *this->element; }

    friend std::ostream& operator<<(std::ostream &o, const Nucleus &nuc) {
        o << std::endl << *nuc.element << "   Z = ";
        o << nuc.charge << std::endl << "   ";
        for (int i = 0; i < 3; i++) {
            o << nuc.coord[i] << " ";
        }
        o << std::endl;
        return o;
    }

private:
    double charge;
    double coord[3];
    const Element *element;
};

class Nuclei : public std::vector<Nucleus> {
public:
    void push_back(const Nucleus &nuc) {
        std::vector<Nucleus>::push_back(nuc);
    }
    void push_back(const char *sym, const double *r) {
        PeriodicTable pt;
        Nucleus nuc(pt.getElement(sym), r);
        std::vector<Nucleus>::push_back(nuc);
    }
};

} //namespace mrchem
