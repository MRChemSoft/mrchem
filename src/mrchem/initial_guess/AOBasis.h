#pragma once

#include <vector>
#include "AOContraction.h"

class AOBasis {
public:
    AOBasis();
    AOBasis(const AOBasis &bas);
    virtual ~AOBasis();

    void append(const AOContraction &ctr);

    GaussExp<3> getAO(int n, const double *center) const;
    GaussExp<3> getBasis(const double *center) const;
    GaussExp<3> getNormBasis(const double *center) const;

    AOContraction &getContraction(int n) { return *this->ctrs[n]; }
    const AOContraction &getContraction(int n) const { return *this->ctrs[n]; }

    int size() const { return this->ctrs.size(); }
    int getNFunc() const { return this->nFunc; }

    // This should print shell by shell
    friend std::ostream& operator<<(std::ostream &o, const AOBasis &b) {
        o << "    nFunc " << b.nFunc << std::endl;
        for (unsigned int i = 0; i < b.ctrs.size(); i++) {
            o << "    " << *b.ctrs[i];
        }
        return o;
    }
private:
    int nPrim; ///< Total number of primitives in set
    int nFunc; ///< Total number of functions (all l-components included)
    std::vector<AOContraction *> ctrs;
};

