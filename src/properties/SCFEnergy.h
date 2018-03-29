#pragma once

#include "MRCPP/Printer"

#include "mrchem.h"

namespace mrchem {

class SCFEnergy final {
public:
    SCFEnergy(double nuc = 0.0, double el = 0.0,
              double orb = 0.0, double kin = 0.0,
              double en = 0.0,  double ee = 0.0,
              double xc = 0.0,  double x = 0.0) :
        E_nuc(nuc), E_el(el), E_orb(orb), E_kin(kin),
        E_en(en), E_ee(ee), E_x(x), E_xc(xc){ }
    ~SCFEnergy() { }

    double getTotalEnergy() const { return this->E_nuc + this->E_el; }
    double getNuclearEnergy() const { return this->E_nuc; }
    double getElectronicEnergy() const { return this->E_el; }

    double getOrbitalEnergy() const { return this->E_orb; }
    double getKineticEnergy() const { return this->E_kin; }
    double getElectronNuclearEnergy() const { return this->E_en; }
    double getElectronElectronEnergy() const { return this->E_ee; }
    double getExchangeCorrelationEnergy() const { return this->E_xc; }
    double getExchangeEnergy() const { return this->E_x; }

    friend std::ostream& operator<<(std::ostream &o, const SCFEnergy &en) {
        double E_au = en.E_nuc + en.E_el;
        double E_kJ = E_au   * PHYSCONST::kJ;
        double E_kcal = E_au * PHYSCONST::kcal;
        double E_eV = E_au   * PHYSCONST::eV;
        int oldPrec = mrcpp::Printer::setPrecision(15);
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                         SCF Energy                         " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o <<"                                                            "<<std::endl;
        o << " Sum orbital energy:          " << std::setw(29) << en.E_orb  << std::endl;
        o << " Kinetic energy:              " << std::setw(29) << en.E_kin  << std::endl;
        o << " E-N energy:                  " << std::setw(29) << en.E_en   << std::endl;
        o << " Coulomb energy:              " << std::setw(29) << en.E_ee   << std::endl;
        o << " Exchange energy:             " << std::setw(29) << en.E_x    << std::endl;
        o << " X-C energy:                  " << std::setw(29) << en.E_xc   << std::endl;
        o <<"                                                            "<<std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o <<"                                                            "<<std::endl;
        o << " Electronic energy            " << std::setw(29) << en.E_el   << std::endl;
        o << " Nuclear energy               " << std::setw(29) << en.E_nuc  << std::endl;
        o <<"                                                            "<<std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o <<"                                                            "<<std::endl;
        o << " Total energy       (au)      " << std::setw(29) << E_au      << std::endl;
        o << "                    (kJ/mol)  " << std::setw(29) << E_kJ      << std::endl;
        o << "                    (kcal/mol)" << std::setw(29) << E_kcal    << std::endl;
        o << "                    (eV)      " << std::setw(29) << E_eV      << std::endl;
        o <<"                                                            "<<std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(oldPrec);
        return o;
    }

protected:
    double E_nuc;
    double E_el;

    double E_orb;
    double E_kin;
    double E_en;
    double E_ee;
    double E_x;
    double E_xc;
};

} //namespace mrchem
