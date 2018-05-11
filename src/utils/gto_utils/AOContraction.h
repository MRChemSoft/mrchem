#pragma once

#include <iostream>
#include <vector>

#include "MRCPP/Gaussians"

namespace mrchem {
namespace gto_utils {

static const int MAX_L = 3; // Max f-functions

/** Contracted atomic orbital.
 *
 * AOContraction defines a singel, contracted atomic orbital as a GaussExp
 * with coefficients.
 *
 */
class AOContraction final {
public:
    AOContraction(int l = 0);

    void append(double e, double c);
    mrcpp::GaussExp<3> getNormContraction(int m, const double *center) const;
    mrcpp::GaussExp<3> getContraction(int m, const double *center) const;

    int getNComp() const { return this->nComp; }
    int getMoment() const { return this->L; }
    void setMoment(int l) { this->L = l; }

    friend std::ostream& operator<<(std::ostream &o, const AOContraction &c) {
        o << " " << c.L << " " << c.coefs.size() << std::endl;
        for (unsigned int i = 0; i < c.expo.size(); i++) {
            o << "    " << c.expo[i] << "   " << c.coefs[i] << std::endl;
        }
        return o;
    }
protected:
    int L;
    int nComp;
    std::vector<double> expo;
    std::vector<double> coefs;
};

} //namespace gto_utils
} //namespace mrchem
