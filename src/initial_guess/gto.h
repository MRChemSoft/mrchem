#pragma once

#include <string>

#include "qmfunctions.h"
#include "chemistry.h"

/** @file gto.h
 *
 * @brief Module for generating initial guess of GTO molecular orbitals
 *
 * The initial_guess::gto namespace provides functionality to setup an
 * initial guess of GTO molecular orbitals. The initial guess requires
 * external input to provide basis set information and MO coefficients.
 */

namespace mrchem {

namespace gto_utils {
class OrbitalExp;
}

namespace initial_guess {
namespace gto {

OrbitalVector setup(double prec,
                    const Molecule &mol,
                    const std::string &bas_file,
                    const std::string &mo_file);
OrbitalVector setup(double prec,
                    const Molecule &mol,
                    const std::string &bas_file,
                    const std::string &moa_file,
                    const std::string &mob_file);
OrbitalVector project_ao(double prec,
                         const std::string &bas_file,
                         int spin,
                         int occ,
                         int N = -1);
OrbitalVector project_mo(double prec,
                         const std::string &bas_file,
                         const std::string &mo_file,
                         int spin,
                         int occ,
                         int N = -1);
 mrcpp::FunctionTree<3>* project_density(double prec,
                         const Nucleus &nuc,
                         const std::string &bas_file,
                         const std::string &dens_file);

} //namespace gto
} //namespace initial_guess
} //namespace mrchem
