#pragma once

#include "mrchem.h"

namespace mrchem {

namespace SPIN { enum type { Paired, Alpha, Beta }; }
namespace NUMBER { enum type { Total, Real, Imag }; }
namespace DENSITY { enum type { Total, Spin, Alpha, Beta }; }


class Orbital;
typedef std::vector<Orbital> OrbitalVector;
typedef std::vector<std::tuple<int, Orbital> > OrbitalChunk;

namespace orbital {

ComplexDouble dot(Orbital bra, Orbital ket);

bool compare(const Orbital &orb_a, const Orbital &orb_b);
int compare_occ(const Orbital &orb_a, const Orbital &orb_b);
int compare_spin(const Orbital &orb_a, const Orbital &orb_b);

Orbital add(ComplexDouble a, Orbital inp_a, ComplexDouble b, Orbital inp_b, double prec = -1.0);
OrbitalVector add(ComplexDouble a, OrbitalVector &inp_a, ComplexDouble b, OrbitalVector &inp_b, double prec = -1.0);

Orbital multiply(Orbital inp_a, Orbital inp_b, double prec = -1.0);
Orbital multiply(const ComplexVector &c, OrbitalVector &inp, double prec = -1.0);
OrbitalVector multiply(const ComplexMatrix &U, OrbitalVector &inp, double prec = -1.0);

OrbitalVector deep_copy(OrbitalVector &inp);
OrbitalVector param_copy(const OrbitalVector &inp);

OrbitalVector adjoin(OrbitalVector &inp_a, OrbitalVector &inp_b);
OrbitalVector disjoin(OrbitalVector &inp, int spin);

void save_orbitals(OrbitalVector &Phi, const std::string &file, int n_orbs = -1);
OrbitalVector load_orbitals(const std::string &file, int n_orbs = -1);

void free(OrbitalVector &vec);
void normalize(OrbitalVector &vec);
void orthogonalize(OrbitalVector &vec);
void orthogonalize(OrbitalVector &vec, OrbitalVector &inp);

ComplexMatrix calc_overlap_matrix(OrbitalVector &braket);
ComplexMatrix calc_overlap_matrix(OrbitalVector &bra, OrbitalVector &ket);
ComplexMatrix calc_lowdin_matrix(OrbitalVector &Phi);

ComplexMatrix localize(double prec, OrbitalVector &Phi);
ComplexMatrix diagonalize(double prec, OrbitalVector &Phi, ComplexMatrix &F);
ComplexMatrix orthonormalize(double prec, OrbitalVector &Phi);

int size_empty(const OrbitalVector &vec);
int size_occupied(const OrbitalVector &vec);
int size_singly(const OrbitalVector &vec);
int size_doubly(const OrbitalVector &vec);
int size_paired(const OrbitalVector &vec);
int size_alpha(const OrbitalVector &vec);
int size_beta(const OrbitalVector &vec);
int get_multiplicity(const OrbitalVector &vec);
int get_electron_number(const OrbitalVector &vec, int spin = SPIN::Paired);

void set_spins(OrbitalVector &vec, const IntVector &spins);
void set_errors(OrbitalVector &vec, const DoubleVector &errors);
void set_occupancies(OrbitalVector &vec, const IntVector &occ);

IntVector get_spins(const OrbitalVector &vec);
IntVector get_occupancies(const OrbitalVector &vec);
DoubleVector get_errors(const OrbitalVector &vec);
DoubleVector get_norms(const OrbitalVector &vec);
DoubleVector get_squared_norms(const OrbitalVector &vec);

void print(const OrbitalVector &vec);

} //namespace orbital


typedef mrcpp::FunctionTree<3> Density;
typedef mrcpp::FunctionTreeVector<3> DensityVector;
namespace density {

void compute(double prec, Density &rho, Orbital phi, int spin);
void compute(double prec, Density &rho, OrbitalVector &Phi, int spin);

} //namespace density


} //namespace mrchem
