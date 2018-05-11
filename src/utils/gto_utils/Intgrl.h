/** Class to parse and process basis sets in MOLECULE or INTGRL format.
 */
#pragma once

#include "MRCPP/Gaussians"

#include <vector>
#include <string>
#include <iostream>

namespace mrchem {
class Nucleus;

namespace gto_utils {
class AOBasis;

class Intgrl final {
public:
    Intgrl(const std::string &file);
    ~Intgrl();

    int getNNuclei() const { return this->nuclei.size(); }

    const Nucleus &getNucleus(int i) const { return *this->nuclei[i]; }
    const AOBasis &getAOBasis(int i) const { return *this->basis[i]; }

    Nucleus &getNucleus(int i) { return *this->nuclei[i]; }
    AOBasis &getAOBasis(int i) { return *this->basis[i]; }

    mrcpp::GaussExp<3> getAtomBasis(int i, bool norm = true) const;
    mrcpp::GaussExp<3> getMolBasis(bool norm = true) const;

protected:
    std::vector<Nucleus *> nuclei;
    std::vector<AOBasis *> basis;

    void readIntgrlFile(std::iostream &ifs);
    void readContractionBlock(std::iostream &ifs, AOBasis &basis, int l);
    void readAtomBlock(std::iostream &ifs);
    void readAtomData(std::iostream &ifs, int n_atoms, double z);
};

} //namespace gto_utils
} //namespace mrchem
