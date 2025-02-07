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

#include "catch2/catch_all.hpp"

#include "MRCPP/Parallel"
#include "mrchem.h"

using namespace mrchem;
using MATHCONST::pi;

namespace qmfunction_tests {

auto f = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-1.0 * R * R);
};

auto g = [](const mrcpp::Coord<3> &r) -> double {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-2.0 * R * R);
};

TEST_CASE("QMFunction", "[qmfunction]") {
    const double prec = 1.0e-3;

    SECTION("copy non-shared function") {
        mrcpp::ComplexFunction func_1(false);
        mrcpp::cplxfunc::project(func_1, f, NUMBER::Real, prec);

        SECTION("copy constructor") {
            mrcpp::ComplexFunction func_2(func_1);
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(&func_2.real() == &func_1.real());
            REQUIRE(&func_2.imag() == &func_1.imag());
        }

        SECTION("default constructor plus assignment") {
            mrcpp::ComplexFunction func_2;
            func_2 = func_1;
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(&func_2.real() == &func_1.real());
            REQUIRE(&func_2.imag() == &func_1.imag());
        }

        SECTION("assigment constructor") {
            mrcpp::ComplexFunction func_2 = func_1;
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(&func_2.real() == &func_1.real());
            REQUIRE(&func_2.imag() == &func_1.imag());
        }

        SECTION("deep copy to non-shared") {
            mrcpp::ComplexFunction func_2(false);
            mrcpp::cplxfunc::deep_copy(func_2, func_1);
            REQUIRE(not func_2.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(&func_2.real() != &func_1.real());
            REQUIRE(&func_2.imag() == &func_1.imag());
        }
#ifdef MRCHEM_HAS_MPI
        SECTION("deep copy to shared") {
            mrcpp::ComplexFunction func_2(true);
            mrcpp::cplxfunc::deep_copy(func_2, func_1);
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(&func_2.real() != &func_1.real());
            REQUIRE(&func_2.imag() == &func_1.imag());
        }
#endif
    }

#ifdef MRCHEM_HAS_MPI
    SECTION("copy shared function") {
        mrcpp::ComplexFunction func_1(true);
        mrcpp::cplxfunc::project(func_1, f, NUMBER::Real, prec);
        mrcpp::cplxfunc::project(func_1, g, NUMBER::Imag, prec);

        SECTION("copy constructor") {
            mrcpp::ComplexFunction func_2(func_1);
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.integrate().real()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(func_1.integrate().imag()));
        }

        SECTION("default constructor plus assignment") {
            mrcpp::ComplexFunction func_2;
            func_2 = func_1;
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.integrate().real()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(func_1.integrate().imag()));
        }

        SECTION("assigment constructor") {
            mrcpp::ComplexFunction func_2 = func_1;
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.integrate().real()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(func_1.integrate().imag()));
        }

        SECTION("deep copy to non-shared") {
            mrcpp::ComplexFunction func_2(false);
            mrcpp::cplxfunc::deep_copy(func_2, func_1);
            REQUIRE(not func_2.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.integrate().real()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(func_1.integrate().imag()));
        }

        SECTION("deep copy to shared") {
            mrcpp::ComplexFunction func_2(true);
            mrcpp::cplxfunc::deep_copy(func_2, func_1);
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.integrate().real()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(func_1.integrate().imag()));
        }
    }
#endif

    SECTION("rescale non-shared function") {
        mrcpp::ComplexFunction func(false);
        mrcpp::cplxfunc::project(func, f, NUMBER::Real, prec);
        mrcpp::cplxfunc::project(func, g, NUMBER::Imag, prec);

        const double ref_norm = func.norm();
        const double f_int = func.real().integrate();
        const double g_int = func.imag().integrate();
        SECTION("real scalar") {
            func.rescale(pi);
            REQUIRE(func.norm() == Catch::Approx(pi * ref_norm));
            REQUIRE(func.real().integrate() == Catch::Approx(pi * f_int));
            REQUIRE(func.imag().integrate() == Catch::Approx(pi * g_int));
        }
        SECTION("imaginary unit") {
            ComplexDouble i(0.0, 1.0);
            func.rescale(i);
            REQUIRE(func.norm() == Catch::Approx(ref_norm));
            REQUIRE(func.real().integrate() == Catch::Approx(-g_int));
            REQUIRE(func.imag().integrate() == Catch::Approx(f_int));
        }
        SECTION("unitary rotation") {
            double re = std::sin(0.5);
            double im = std::cos(0.5);
            ComplexDouble c(re, im);
            func.rescale(c);
            REQUIRE(func.norm() == Catch::Approx(ref_norm));
            REQUIRE(func.real().integrate() == Catch::Approx(re * f_int - im * g_int));
            REQUIRE(func.imag().integrate() == Catch::Approx(im * f_int + re * g_int));
        }
    }
#ifdef MRCHEM_HAS_MPI
    SECTION("rescale shared function") {
        mrcpp::ComplexFunction func(true);
        mrcpp::cplxfunc::project(func, g, NUMBER::Real, prec);
        mrcpp::cplxfunc::project(func, f, NUMBER::Imag, prec);

        const double ref_norm = func.norm();
        const double f_int = func.real().integrate();
        const double g_int = func.imag().integrate();

        SECTION("real scalar") {
            func.rescale(pi);
            REQUIRE(func.norm() == Catch::Approx(pi * ref_norm));
            REQUIRE(func.real().integrate() == Catch::Approx(pi * f_int));
            REQUIRE(func.imag().integrate() == Catch::Approx(pi * g_int));
        }
        SECTION("imaginary unit") {
            ComplexDouble i(0.0, 1.0);
            func.rescale(i);
            mrcpp::mpi::barrier(mrcpp::mpi::comm_share);
            REQUIRE(func.norm() == Catch::Approx(ref_norm));
            REQUIRE(func.real().integrate() == Catch::Approx(-g_int));
            REQUIRE(func.imag().integrate() == Catch::Approx(f_int));
        }
        SECTION("unitary rotation") {
            double re = std::sin(0.5);
            double im = std::cos(0.5);
            ComplexDouble c(re, im);
            func.rescale(c);
            REQUIRE(func.norm() == Catch::Approx(ref_norm));
            REQUIRE(func.real().integrate() == Catch::Approx(re * f_int - im * g_int));
            REQUIRE(func.imag().integrate() == Catch::Approx(im * f_int + re * g_int));
        }
    }

    SECTION("add shared function") {
        mrcpp::ComplexFunction f_re(false);
        mrcpp::ComplexFunction f_im(true);
        mrcpp::cplxfunc::project(f_re, f, NUMBER::Real, prec);
        mrcpp::cplxfunc::project(f_im, f, NUMBER::Imag, prec);

        SECTION("into non-shared function") {
            ComplexDouble c(0.5, 0.5);
            mrcpp::ComplexFunction func_h(false);
            SECTION("with complex scalar") {
                mrcpp::cplxfunc::add(func_h, c, f_re, c, f_im, -1.0);
                REQUIRE(func_h.integrate().real() == Catch::Approx(0.0));
                REQUIRE(func_h.integrate().imag() == Catch::Approx(f_im.integrate().imag()));
            }
            SECTION("with function conjugate") {
                mrcpp::cplxfunc::add(func_h, c, f_re, c, f_im.dagger(), -1.0);
                REQUIRE(func_h.integrate().real() == Catch::Approx(f_re.integrate().real()));
                REQUIRE(func_h.integrate().imag() == Catch::Approx(0.0));
            }
        }
        SECTION("into shared function") {
            ComplexDouble c(0.5, 0.5);
            mrcpp::ComplexFunction func_h(true);
            SECTION("with complex scalar") {
                mrcpp::cplxfunc::add(func_h, c, f_re, c, f_im, -1.0);
                REQUIRE(func_h.integrate().real() == Catch::Approx(0.0));
                REQUIRE(func_h.integrate().imag() == Catch::Approx(f_im.integrate().imag()));
            }
            SECTION("with function conjugate") {
                mrcpp::cplxfunc::add(func_h, c, f_re, c, f_im.dagger(), -1.0);
                REQUIRE(func_h.integrate().real() == Catch::Approx(f_re.integrate().real()));
                REQUIRE(func_h.integrate().imag() == Catch::Approx(0.0));
            }
        }
    }
#endif

    SECTION("multiply non-shared function") {
        mrcpp::ComplexFunction func_1(false);
        mrcpp::cplxfunc::project(func_1, f, NUMBER::Real, prec);
        mrcpp::cplxfunc::project(func_1, g, NUMBER::Imag, prec);

        SECTION("into non-shared function") {
            mrcpp::ComplexFunction func_2(false);
            mrcpp::cplxfunc::multiply(func_2, func_1, func_1.dagger(), -1.0);
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.squaredNorm()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(0.0));
        }
#ifdef MRCHEM_HAS_MPI
        SECTION("into shared function") {
            mrcpp::ComplexFunction func_2(true);
            mrcpp::cplxfunc::multiply(func_2, func_1, func_1.dagger(), -1.0);
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.squaredNorm()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(0.0));
        }
#endif
    }

#ifdef MRCHEM_HAS_MPI
    SECTION("multiply shared function") {
        mrcpp::ComplexFunction func_1(true);
        mrcpp::cplxfunc::project(func_1, f, NUMBER::Real, prec);
        mrcpp::cplxfunc::project(func_1, g, NUMBER::Imag, prec);

        SECTION("into non-shared function") {
            mrcpp::ComplexFunction func_2(false);
            mrcpp::cplxfunc::multiply(func_2, func_1, func_1.dagger(), -1.0);
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.squaredNorm()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(0.0));
        }
        SECTION("into shared function") {
            mrcpp::ComplexFunction func_2(true);
            mrcpp::cplxfunc::multiply(func_2, func_1, func_1.dagger(), -1.0);
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.squaredNorm()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(0.0));
        }
    }
#endif
}

} // namespace qmfunction_tests
