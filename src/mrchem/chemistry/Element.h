#pragma once

#include "TelePrompter.h"

/** Basic chemical element data. */
class Element {
public:
    Element(int z, const char *s, const char *n, double m, double vdw, double cov, double g)
        : Z(z), symbol(s), name(n), mass(m), r_vdw(vdw), r_cov(cov), g_val(g) { }
    virtual ~Element() { }

    const std::string &getSymbol() const { return this->symbol; }
    const std::string &getName() const { return this->name; }

    int getZ() const { return this->Z; }
    double getMass() const { return this->mass; }
    double getVdw() const { return this->r_vdw; }
    double getCov() const { return this->r_cov; }
    double getGValue() const { if (this->g_val < 0.0) MSG_ERROR("g-value unavailable"); return this->g_val; }

    friend std::ostream& operator<<(std::ostream &o, const Element &e) {
        o << e.symbol;
        return o;
    }

protected:
    const int Z;                /** atomic number */
    const std::string symbol;   /** atomic symbol */
    const std::string name;     /** element name */
    const double mass;          /** atomic mass */
    const double r_vdw;         /** van der waals radius */
    const double r_cov;         /** covalent radius */
    const double g_val;         /** g-value */
};

