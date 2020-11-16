#include "Intgrl.h"

#include <string>

#include <MRCPP/Printer>

#include "chemistry/Nucleus.h"
#include "AOBasis.h"

using mrcpp::GaussExp;

namespace mrchem {
namespace gto_utils {

Intgrl::Intgrl(const std::string &file) {
    std::ifstream ifs(file.c_str());
    if (not ifs) MSG_ABORT("Failed to open basis set file: " << file);
    readIntgrlFile(file, ifs);
    ifs.close();
}

Intgrl::~Intgrl() {
    for (auto &i : this->nuclei) {
        if (i != nullptr) { delete i; }
    }
    for (auto &basi : this->basis) {
        if (basi != nullptr) { delete basi; }
    }
}

void Intgrl::readIntgrlFile(const std::string &fname, std::ifstream &ifs) {
    int nTypes;
    std::string line;
    // First line is a comment
    if (not std::getline(ifs, line)) MSG_ABORT("Unexpected end of file " << fname << "  while reading basis sets!");
    // Second line is number of atom types
    if (not std::getline(ifs, line)) MSG_ABORT("Unexpected end of file " << fname << "  while reading basis sets!");
    std::istringstream iss(line);
    iss >> nTypes;
    for (int i = 0; i < nTypes; i++) { readAtomBlock(ifs); }
}

void Intgrl::readAtomBlock(std::ifstream &ifs) {
    double z;
    int nAtoms;
    int nShell;

    ifs >> z;
    ifs >> nAtoms;
    ifs >> nShell;

    int funcsPerShell[nShell];
    for (int i = 0; i < nShell; i++) { ifs >> funcsPerShell[i]; }

    readAtomData(ifs, nAtoms, z);
    AOBasis bas;
    for (int i = 0; i < nShell; i++) {
        for (int j = 0; j < funcsPerShell[i]; j++) { readContractionBlock(ifs, bas, i); }
    }
    for (int i = 0; i < nAtoms; i++) {
        auto *aoBas = new AOBasis(bas);
        this->basis.push_back(aoBas);
    }
}

void Intgrl::readAtomData(std::ifstream &ifs, int n_atoms, double z) {
    mrcpp::Coord<3> coord;
    std::string sym;
    for (int j = 0; j < n_atoms; j++) {
        ifs >> sym;
        if (sym.size() > 2) { sym = sym.erase(2); }
        for (int d = 0; d < 3; d++) { ifs >> coord[d]; }
        PeriodicTable pt;
        const Element &element = pt.getElement(sym.c_str());

        auto *nuc = new Nucleus(element, coord);
        nuc->setCharge(z);
        this->nuclei.push_back(nuc);
    }
}

void Intgrl::readContractionBlock(std::ifstream &ifs, AOBasis &bas, int l) {
    if (l > 2) MSG_ABORT("Only s, p and d orbitals are currently supported");
    int nprim, nctr;
    ifs >> nprim;
    ifs >> nctr;

    int start = bas.size();
    int nbas = 0;
    for (int i = 0; i < nctr; i++) {
        AOContraction ctr(l);
        bas.append(ctr);
        nbas += (l + 1) * (l + 2) / 2;
    }

    double expo;
    double coef;
    for (int k = 0; k < nprim; k++) {
        ifs >> expo;
        for (int i = 0; i < nctr; i++) {
            ifs >> coef;
            if (std::abs(coef) > mrcpp::MachineZero) { bas.getContraction(start + i).append(expo, coef); }
        }
    }
}

GaussExp<3> Intgrl::getAtomBasis(int i, bool norm) const {
    assert(i >= 0 and i < this->nuclei.size());
    const mrcpp::Coord<3> &coord = this->nuclei[i]->getCoord();
    if (norm) {
        return this->basis[i]->getNormBasis(coord);
    } else {
        return this->basis[i]->getBasis(coord);
    }
}

GaussExp<3> Intgrl::getMolBasis(bool norm) const {
    GaussExp<3> molexp;
    for (unsigned int i = 0; i < this->nuclei.size(); i++) {
        const mrcpp::Coord<3> &coord = this->nuclei[i]->getCoord();
        if (norm) {
            molexp.append(this->basis[i]->getNormBasis(coord));
        } else {
            molexp.append(this->basis[i]->getBasis(coord));
        }
    }
    return molexp;
}

} // namespace gto_utils
} // namespace mrchem
