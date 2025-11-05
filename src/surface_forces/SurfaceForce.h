#ifndef SURFACEFORCE_H
#define SURFACEFORCE_H

#include "chemistry/Molecule.h"
#include "qmfunctions/Orbital.h"
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace surface_force {

// Function declaration
Eigen::MatrixXd surface_forces(mrchem::Molecule &mol, mrchem::OrbitalVector &Phi, double prec, const json &json_fock, std::string leb_prec, double radius_factor);

} // namespace surface_force

#endif // SURFACEFORCE_H
