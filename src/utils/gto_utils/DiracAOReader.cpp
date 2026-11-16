/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include "utils/gto_utils/DiracAOReader.h"

#include "chemistry/PeriodicTable.h"
#include "chemistry/PhysicalConstants.h"

#include <MRCPP/Printer>

#include <cmath>
#include <complex>
#include <filesystem>

namespace mrchem {
namespace gto_utils {

/**
 * @brief Local helper function to read a vector of values from an HDF5 group.
 * @tparam T The type of the values to read.
 * @param g The HDF5 group to read from.
 * @param name The name of the dataset to read.
 * @param type The HDF5 data type of the dataset.
 * @return A vector containing the values read from the dataset.
 */
template <typename T>
static std::vector<T> readVec(const H5::Group &g, const std::string &name, const H5::PredType &type) {
    H5::DataSet ds = g.openDataSet(name);
    std::vector<T> v(ds.getSpace().getSimpleExtentNpoints());
    ds.read(v.data(), type);
    return v;
}

/**
 * @brief Local helper function to read a vector of int64_t values from an HDF5 group.
 * @param g The HDF5 group to read from.
 * @param name The name of the dataset to read.
 * @return A vector containing the int64_t values read from the dataset.
 */
static std::vector<int64_t> readInts(const H5::Group &g, const std::string &n) { return readVec<int64_t>(g, n, H5::PredType::NATIVE_INT64); }

/**
 * @brief Local helper function to read a vector of double values from an HDF5 group.
 * @param g The HDF5 group to read from.
 * @param name The name of the dataset to read.
 * @return A vector containing the double values read from the dataset.
 */
static std::vector<double> readDoubles(const H5::Group &g, const std::string &n) { return readVec<double>(g, n, H5::PredType::NATIVE_DOUBLE); }

DiracAOReader::DiracAOReader(const std::string &basisFile) {
    if (not std::filesystem::is_regular_file(basisFile))
        MSG_ABORT("Failed to open basis set file: " << basisFile);
    this->AOFile = basisFile;
    this->cartesian = true;
    readBasisFile();
}

void DiracAOReader::readBasisFile() {
    try {
        H5::Exception::dontPrint();
        const H5::H5File file(this->AOFile, H5F_ACC_RDONLY);
        readAtomBlock(file);
        this->basis = readAOBasis(file.openGroup("/input/aobasis/1")); // large component, always present
        // Small component only exists for 4-component calculations (not for X2C / 2c Hamiltonians)
        if (file.nameExists("/input/aobasis/2"))
            this->smallBasis = readAOBasis(file.openGroup("/input/aobasis/2"));
    } catch (const H5::Exception &e) {
        MSG_ABORT("Error reading DIRAC file " << this->AOFile << ": " << e.getDetailMsg());
    }
}

void DiracAOReader::readAtomBlock(const H5::H5File &file) {
    H5::Group mol = file.openGroup("/input/molecule");
    const size_t nAtoms = readInts(mol, "n_atoms")[0];
    const std::vector<double> charge = readDoubles(mol, "nuc_charge");
    const std::vector<double> geom = readDoubles(mol, "geometry");
    if (charge.size() != nAtoms or geom.size() != 3 * nAtoms)
        MSG_ABORT("Inconsistent molecule data in DIRAC file " << this->AOFile);

    const double angToBohr = PhysicalConstants::get("angstrom2bohrs");
    PeriodicTable pt;
    for (size_t i = 0; i < nAtoms; i++) {
        mrcpp::Coord<3> coord{geom[3 * i] * angToBohr, geom[3 * i + 1] * angToBohr, geom[3 * i + 2] * angToBohr};
        const Element &element = pt.getElement(static_cast<int>(std::lround(charge[i])));
        this->nuclei.push_back(std::make_shared<Nucleus>(element, charge[i], coord));
    }
}

std::vector<std::shared_ptr<AOBasis>> DiracAOReader::readAOBasis(const H5::Group &g) const {
    if (readInts(g, "angular")[0] != 1)
        MSG_ABORT("Only Cartesian DIRAC basis sets are supported");

    const size_t nShells = readInts(g, "n_shells")[0];
    const size_t nAO = readInts(g, "n_ao")[0];
    const std::vector<int64_t> orbmom = readInts(g, "orbmom");      // l + 1, one per shell
    const std::vector<int64_t> nPrim = readInts(g, "n_prim");       // primitives per shell
    const std::vector<int64_t> nCont = readInts(g, "n_cont");       // contracted functions per shell (always 1)
    const std::vector<double> expo = readDoubles(g, "exponents");   // all primitives, shell after shell
    const std::vector<double> coef = readDoubles(g, "contractions"); // matching coefficients
    const std::vector<double> center = readDoubles(g, "center");     // x, y, z per shell, in bohr
    if (orbmom.size() != nShells or nPrim.size() != nShells or nCont.size() != nShells or center.size() != 3 * nShells)
        MSG_ABORT("Inconsistent shell data in DIRAC basis");

    std::vector<AOBasis> bas(this->nuclei.size());
    size_t offset = 0;   // position of the current shell's first primitive in expo/coef
    size_t lastAtom = 0;
    size_t aoCount = 0;
    for (size_t k = 0; k < nShells; k++) {
        // DIRAC writes general contractions as separate shells with n_cont = 1
        if (nCont[k] != 1)
            MSG_ABORT("DIRAC shells with n_cont > 1 are not supported");

        const int l = static_cast<int>(orbmom[k]) - 1;
        const size_t atom = findAtom({center[3 * k], center[3 * k + 1], center[3 * k + 2]});
        // MO coefficient rows follow the shell order of the file, OrbitalExp builds AOs atom by atom
        if (atom < lastAtom)
            MSG_ABORT("DIRAC shells are not grouped by atom, AO order would not match the MO coefficients");
        lastAtom = atom;

        // The stored coefficients already contain DIRAC's primitive normalisation: no MRChem normalisation
        AOContraction ctr(l, false);
        for (int64_t p = 0; p < nPrim[k]; p++) {
            if (coef[offset + p] != 0.0)
                ctr.append(expo[offset + p], coef[offset + p]);
        }
        offset += nPrim[k];

        bas[atom].append(ctr);
        aoCount += (l + 1) * (l + 2) / 2; // Cartesian components
    }
    if (offset != expo.size() or offset != coef.size())
        MSG_ABORT("Inconsistent number of primitives in DIRAC basis");
    if (aoCount != nAO)
        MSG_ABORT("DIRAC basis has " << aoCount << " Cartesian AOs, expected n_ao = " << nAO);

    std::vector<std::shared_ptr<AOBasis>> out;
    for (const AOBasis &b : bas)
        out.push_back(std::make_shared<AOBasis>(b));
    return out;
}

size_t DiracAOReader::findAtom(const mrcpp::Coord<3> &r) const {
    // Nuclear positions come from /input/molecule/geometry (Angstrom, 8 decimals) converted to bohr,
    // so they differ from the shell centres (bohr) by ~1e-9: compare with a loose tolerance
    const double tol = 1.0e-6;
    for (size_t i = 0; i < this->nuclei.size(); i++) {
        const mrcpp::Coord<3> &R = this->nuclei[i]->getCoord();
        if (std::abs(R[0] - r[0]) < tol and std::abs(R[1] - r[1]) < tol and std::abs(R[2] - r[2]) < tol)
            return i;
    }
    MSG_ABORT("DIRAC shell centre (" << r[0] << ", " << r[1] << ", " << r[2] << ") does not match any atom");
}

std::vector<ComplexMatrix> DiracAOReader::getMOCoefficients(const std::string &MOFile, size_t nMOs, size_t nSpinors) const {
    if (nSpinors != 2 and nSpinors != 4)
        MSG_ABORT("DIRAC MO coefficients need nSpinors = 2 or 4, got " << nSpinors);
    if (nSpinors == 4 and this->smallBasis.empty())
        MSG_ABORT("DIRAC file " << this->AOFile << " has no small component (not a 4-component calculation)");

    std::vector<double> q; // all coefficients in one flat array
    size_t nBasis = 0, nMO = 0, nPO = 0; // number of basis functions, number of stored MOs, number of negative-energy MOs
    try {
        H5::Exception::dontPrint();
        const H5::H5File file(MOFile, H5F_ACC_RDONLY);
        const H5::Group mo = file.openGroup("/result/wavefunctions/scf/mobasis");
        nBasis = readInts(mo, "n_basis")[0]; // number of AOs (large + small)
        nMO = readInts(mo, "n_mo")[0];       // number of stored MOs
        nPO = readInts(mo, "n_po")[0];       // the first nPO MOs are negative-energy: skip them
        q = readDoubles(mo, "orbitals");     // all coefficients in one flat array
    } catch (const H5::Exception &e) {
        MSG_ABORT("Error reading DIRAC file " << MOFile << ": " << e.getDetailMsg());
    }

    if (nMOs == 0)
        nMOs = 2 * (nMO - nPO); // all positive-energy spinors (for X2C)

    if (nPO + (nMOs + 1) / 2 > nMO)
        MSG_ABORT("Requested " << nMOs << " spinors, but DIRAC file " << MOFile << " has only " << 2 * (nMO - nPO));

    ComplexMatrix alpha(nBasis, nMOs);
    ComplexMatrix beta(nBasis, nMOs);
    for (size_t i = 0; i < nMOs; i++) {
        size_t m = nPO + i / 2; // spinors 0,1 come from MO nPO, spinors 2,3 from MO nPO+1, ...
        for (size_t a = 0; a < nBasis; a++) {
            ComplexDouble A(q[a + nBasis * (m + nMO * 0)], q[a + nBasis * (m + nMO * 1)]);
            ComplexDouble B(q[a + nBasis * (m + nMO * 2)], q[a + nBasis * (m + nMO * 3)]);
            if (i % 2 == 0) { // first spinor of the pair
                alpha(a, i) = A;
                beta(a, i) = -std::conj(B);
            } else {          // second spinor of the pair
                alpha(a, i) = B;
                beta(a, i) = std::conj(A);
            }
        }
    }
    size_t nLarge = 0;
    for (const auto &b : this->basis)
        nLarge += b->getNFunc();

    std::vector<ComplexMatrix> C;
    C.push_back(alpha.topRows(nLarge)); // large alpha
    C.push_back(beta.topRows(nLarge));  // large beta
    if (nSpinors == 4) {
        C.push_back(alpha.bottomRows(nBasis - nLarge)); // small alpha
        C.push_back(beta.bottomRows(nBasis - nLarge));  // small beta
    }
    return C;
}

std::shared_ptr<AOBasis> DiracAOReader::getSmallAOBasis(size_t i) const {
    if (this->smallBasis.empty())
        return nullptr;
    if (i >= this->smallBasis.size())
        MSG_ABORT("Requested small component basis for atom " << i << ", but DIRAC file " << this->AOFile << " has only " << this->smallBasis.size() << " atoms");
    return this->smallBasis[i];
}

} // namespace gto_utils
} // namespace mrchem