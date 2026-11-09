#include "MagneticPerturbation.h"

#include <MRCPP/MWFunctions>
#include <MRCPP/Printer>
#include <MRCPP/Timer> 

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/NablaOperator.h"
#include "qmoperators/one_electron/PositionOperator.h"
#include "utils/print_utils.h"

namespace mrchem {

MagneticPerturbation::MagneticPerturbation(
        double B3, double c,
        std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
        const mrcpp::Coord<3> &o)
    : para_coeff_(B3 / (2.0 * c))
    , origin_(o)
    , D_(D)
{
    // Diamagnetic part: W_dia = –(B3²/8c²)·(x²+y²)
    // H_BB_dia(o)[2][2] = (1/4)·(x²+y²)  →  coefficient = –(B3²/2c²)
    dia_part_ = -(B3 * B3 / (2.0 * c * c)) * H_BB_dia(o)[2][2];
    dia_part_.name() = "W_B_dia";
    name() = "W_B";
}

void MagneticPerturbation::setup(double prec) {
    prec_ = prec;
    dia_part_.setup(prec);
}

void MagneticPerturbation::clear() {
    dia_part_.clear();
    prec_ = -1.0;
}

/** Apply W to a single (complex) orbital.
 *
 *  W·φ  =  –(B3²/8c²)·(x²+y²)·φ   [dia, real × complex → complex]
 *         + (i·B3/2c)·Lz·φ          [para, i × real_function → imaginary]
 *
 *  For φ = φ_R + i·φ_I:
 *    Re(W·φ) = –(B3²/8c²)·(x²+y²)·φ_R  –  (B3/2c)·Lz·φ_I
 *    Im(W·φ) = –(B3²/8c²)·(x²+y²)·φ_I  +  (B3/2c)·Lz·φ_R
 */
Orbital MagneticPerturbation::operator()(Orbital inp) {
    if (inp.getNNodes() == 0) return inp.paramCopy(false);

    // Diamagnetic: -(B3²/8c²)·(x²+y²)·φ
    Orbital dia_phi = dia_part_(inp);

    // Paramagnetic: Lz = (x-ox)·∂y – (y-oy)·∂x
    NablaOperator nabla(D_);
    nabla.setup(prec_);
    PositionOperator r_op(origin_);
    r_op.setup(prec_);

    auto grad = nabla(inp);
    Orbital x_dy = r_op[0](grad[1]);
    Orbital y_dx = r_op[1](grad[0]);
    nabla.clear();
    r_op.clear();

    // Lz·φ = x_dy - y_dx  (via OrbitalVector API)
    OrbitalVector xdy_v, ydx_v;
    xdy_v.push_back(x_dy);
    ydx_v.push_back(y_dx);
    OrbitalVector Lz_v = orbital::add(ComplexDouble(1.0, 0.0), xdy_v,
                                       ComplexDouble(-1.0, 0.0), ydx_v, prec_);

    // out = dia_phi + i·para_coeff_·Lz_phi
    // Real inp:    out = dia·φ + i·(B3/2c)·Lz·φ ✓
    // Complex inp: orbital::add distributes ComplexDouble coeff correctly ✓
    OrbitalVector dia_v;
    dia_v.push_back(dia_phi);
    OrbitalVector out_v = orbital::add(ComplexDouble(1.0, 0.0), dia_v,
                                        ComplexDouble(0.0, para_coeff_), Lz_v, prec_);

    Orbital out = out_v[0];
    return out;
}

OrbitalVector MagneticPerturbation::operator()(OrbitalVector &inp) {
    OrbitalVector out;
    for (size_t i = 0; i < inp.size(); ++i) {
        mrcpp::Timer t1;
        Orbital out_i;
        if (mrcpp::mpi::my_func(inp[i])) {
            out_i = (*this)(inp[i]);
        } else {
            out_i = inp[i].paramCopy(false);
        }
        out.push_back(out_i);
        print_utils::qmfunction(4, "W_B|" + std::to_string(i) + ">", out_i, t1);
    }
    return out;
}

} // namespace mrchem