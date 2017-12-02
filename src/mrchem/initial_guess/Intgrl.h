/** Class to parse and process basis sets in MOLECULE or INTGRL format.
 */
#pragma once

#include <vector>
#include <string>
#include <iostream>

class Nucleus;
class AOBasis;
template<int> class GaussExp;

class Intgrl {
public:
    Intgrl(const std::string &file);
    virtual ~Intgrl();

    int getNNuclei() const { return this->nuclei.size(); }

    const Nucleus &getNucleus(int i) const { return *this->nuclei[i]; }
    const AOBasis &getAOBasis(int i) const { return *this->basis[i]; }

    Nucleus &getNucleus(int i) { return *this->nuclei[i]; }
    AOBasis &getAOBasis(int i) { return *this->basis[i]; }

    GaussExp<3> getAtomBasis(int i, bool norm = true) const;
    GaussExp<3> getMolBasis(bool norm = true) const;

protected:
    std::vector<Nucleus *> nuclei;
    std::vector<AOBasis *> basis;

    void readIntgrlFile(std::iostream &ifs);
    void readContractionBlock(std::iostream &ifs, AOBasis &basis, int l);
    void readAtomBlock(std::iostream &ifs);
    void readAtomData(std::iostream &ifs, int n_atoms, double z);
};

