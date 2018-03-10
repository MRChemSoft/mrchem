#pragma once

#include "mrchem.h"

namespace mrchem {

namespace SPIN { enum type { Paired, Alpha, Beta }; }
namespace NUMBER { enum type { Total, Real, Imag }; }
namespace DENSITY { enum type { Total, Spin, Alpha, Beta }; }


class Orbital;
typedef std::vector<Orbital> OrbitalVector;

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

void free(OrbitalVector &vec);
void normalize(OrbitalVector &vec);
void orthogonalize(OrbitalVector &vec);
void orthogonalize(OrbitalVector &vec, OrbitalVector &inp);

ComplexMatrix calc_overlap_matrix(OrbitalVector &braket);
ComplexMatrix calc_overlap_matrix(OrbitalVector &bra, OrbitalVector &ket);

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

} //namespace orbital


class Density;
typedef std::vector<Density> DensityVector;
namespace density {

} //namespace density


} //namespace mrchem
