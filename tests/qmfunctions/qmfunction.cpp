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

ComplexDouble i1 = {0.0, 1.0};

auto f = [](const mrcpp::Coord<3> &r) -> ComplexDouble {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-1.0 * R * R);
};

auto g = [](const mrcpp::Coord<3> &r) -> ComplexDouble {
    double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    return std::exp(-2.0 * R * R) * i1;
};

TEST_CASE("QMFunction", "[qmfunction]") {
    const double prec = 1.0e-3;

    SECTION("copy non-shared function") {
        mrcpp::CompFunction func_1;
        mrcpp::project(func_1, f, prec);
        SECTION("copy constructor") {
            mrcpp::CompFunction func_2(func_1);
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(&func_2.complex() == &func_1.complex());
        }

        SECTION("default constructor plus assignment") {
            mrcpp::CompFunction func_2;
            func_2 = func_1;
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(&func_2.complex() == &func_1.complex());
        }

        SECTION("assigment constructor") {
            mrcpp::CompFunction func_2 = func_1;
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(&func_2.complex() == &func_1.complex());
        }

        SECTION("deep copy to non-shared") {
            mrcpp::CompFunction func_2(0, false);
            mrcpp::deep_copy(func_2, func_1);
            REQUIRE(!func_2.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(&func_2.complex() != &func_1.complex());
        }
#ifdef MRCHEM_HAS_MPI
        SECTION("deep copy to shared") {
            mrcpp::CompFunction func_2(0, true);
            mrcpp::deep_copy(func_2, func_1);
            REQUIRE(!func_2.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(&func_2.complex() != &func_1.complex());
        }
#endif
    }

#ifdef MRCHEM_HAS_MPI
    SECTION("copy shared function") {
        mrcpp::CompFunction func_1(0, true);
        mrcpp::project(func_1, f, prec);
        mrcpp::project(func_1, g, prec);

        SECTION("copy constructor") {
            mrcpp::CompFunction func_2(func_1);
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.integrate().real()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(func_1.integrate().imag()));
        }

        SECTION("default constructor plus assignment") {
            mrcpp::CompFunction func_2;
            func_2 = func_1;
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.integrate().real()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(func_1.integrate().imag()));
        }

        SECTION("assigment constructor") {
            mrcpp::CompFunction func_2 = func_1;
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.integrate().real()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(func_1.integrate().imag()));
        }

        SECTION("deep copy to non-shared") {
            mrcpp::CompFunction func_2(0, false);
            mrcpp::deep_copy(func_2, func_1);
            REQUIRE(func_2.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.integrate().real()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(func_1.integrate().imag()));
        }

        SECTION("deep copy to shared") {
            mrcpp::CompFunction func_2(0, true);
            mrcpp::deep_copy(func_2, func_1);
            REQUIRE(func_2.isShared() == func_1.isShared());
            REQUIRE(func_2.norm() == Catch::Approx(func_1.norm()));
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.integrate().real()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(func_1.integrate().imag()));
        }
    }
#endif

    SECTION("rescale non-shared function") {
        mrcpp::CompFunction func;
        mrcpp::project(func, f, prec);
        mrcpp::project(func, g, prec);

        const double ref_norm = func.norm();
        ComplexDouble f_int = func.complex().integrate();
        ComplexDouble g_int = func.complex().integrate();
        SECTION("real scalar") {
            func.rescale(pi);
            REQUIRE(func.norm() == Catch::Approx(pi * ref_norm));
            REQUIRE(std::real(func.complex().integrate()) == Catch::Approx(std::real(pi * f_int)));
            REQUIRE(std::real(func.complex().integrate()) == Catch::Approx(std::real(pi * g_int)));
        }
        SECTION("complexinary unit") {
            ComplexDouble i(0.0, 1.0);
            func.rescale(i);
            REQUIRE(func.norm() == Catch::Approx(ref_norm));
            REQUIRE(std::real(func.complex().integrate()) == Catch::Approx(std::real(i * f_int)));
            REQUIRE(std::imag(func.complex().integrate()) == Catch::Approx(std::imag(i * f_int)));
        }
        SECTION("unitary rotation") {
            double re = std::sin(0.5);
            double im = std::cos(0.5);
            ComplexDouble c(re, im);
            func.rescale(c);
            ComplexDouble i(0.0, 1.0);

            REQUIRE(func.norm() == Catch::Approx(ref_norm));
            REQUIRE(std::real(func.complex().integrate()) == Catch::Approx(std::real(re * f_int + i * im * f_int)));
            REQUIRE(std::imag(func.complex().integrate()) == Catch::Approx(std::imag(re * f_int + i * im * f_int)));
        }
    }
#ifdef MRCHEM_HAS_MPI
    SECTION("rescale shared function") {
        mrcpp::CompFunction func(0, true);
        mrcpp::project(func, g, prec);
        mrcpp::project(func, f, prec);

        const double ref_norm = func.norm();
        const ComplexDouble f_int = func.complex().integrate();
        const ComplexDouble g_int = func.complex().integrate();

        SECTION("real scalar") {
            func.rescale(pi);
            REQUIRE(func.norm() == Catch::Approx(pi * ref_norm));
            REQUIRE(std::real(func.complex().integrate()) == Catch::Approx(std::real(pi * f_int)));
            REQUIRE(std::real(func.complex().integrate()) == Catch::Approx(std::real(pi * g_int)));
        }
        SECTION("complexinary unit") {
            ComplexDouble i(0.0, 1.0);
            func.rescale(i);
            mrcpp::mpi::barrier(mrcpp::mpi::comm_share);
            REQUIRE(func.norm() == Catch::Approx(ref_norm));
            REQUIRE(std::real(func.complex().integrate()) == Catch::Approx(std::real(i * f_int)));
            REQUIRE(std::imag(func.complex().integrate()) == Catch::Approx(std::imag(i * f_int)));
        }
    }

    SECTION("add shared function") {
        mrcpp::CompFunction f_re;
        mrcpp::CompFunction f_im(0, true);
        mrcpp::project(f_re, f, prec);
        mrcpp::project(f_im, g, prec);

        SECTION("into non-shared function") {
            ComplexDouble c(0.5, 0.5);
            mrcpp::CompFunction func_h;
            SECTION("with complex scalar") {
                mrcpp::add(func_h, c, f_re, c, f_im, -1.0);
                REQUIRE(func_h.integrate().real() == Catch::Approx(std::real(std::real(c) * f_re.integrate() - std::imag(c) * std::imag(f_im.integrate()))));
                REQUIRE(func_h.integrate().imag() == Catch::Approx(std::real(std::imag(c) * f_re.integrate() + std::real(c) * std::imag(f_im.integrate()))));
            }
            SECTION("with conjugate flag") {
                mrcpp::add(func_h, c, f_re, c, f_im, -1.0, true);
                REQUIRE(func_h.integrate().real() == Catch::Approx(std::real(std::real(c) * f_re.integrate() + std::imag(c) * std::imag(f_im.integrate()))));
                REQUIRE(func_h.integrate().imag() == Catch::Approx(std::real(std::imag(c) * f_re.integrate() - std::real(c) * std::imag(f_im.integrate()))));
            }
            SECTION("with function conjugate") {
                f_im.dagger();
                mrcpp::add(func_h, c, f_re, c, f_im, -1.0);
                f_im.dagger();
                REQUIRE(func_h.integrate().real() == Catch::Approx(std::real(std::real(c) * f_re.integrate() + std::imag(c) * std::imag(f_im.integrate()))));
                REQUIRE(func_h.integrate().imag() == Catch::Approx(std::real(std::imag(c) * f_re.integrate() - std::real(c) * std::imag(f_im.integrate()))));
            }
        }
        SECTION("into shared function") {
            ComplexDouble c(0.5, 0.5);
            mrcpp::CompFunction func_h(0, true);
            SECTION("with complex scalar") {
                mrcpp::add(func_h, c, f_re, c, f_im, -1.0);
                REQUIRE(func_h.integrate().real() == Catch::Approx(std::real(std::real(c) * f_re.integrate() - std::imag(c) * std::imag(f_im.integrate()))));
                REQUIRE(func_h.integrate().imag() == Catch::Approx(std::real(std::imag(c) * f_re.integrate() + std::real(c) * std::imag(f_im.integrate()))));
            }
            SECTION("with function conjugate") {
                f_im.dagger();
                mrcpp::add(func_h, c, f_re, c, f_im, -1.0);
                f_im.dagger();
                REQUIRE(func_h.integrate().real() == Catch::Approx(std::real(std::real(c) * f_re.integrate() + std::imag(c) * std::imag(f_im.integrate()))));
                REQUIRE(func_h.integrate().imag() == Catch::Approx(std::real(std::imag(c) * f_re.integrate() - std::real(c) * std::imag(f_im.integrate()))));
            }
        }
    }
#endif

    SECTION("multiply non-shared function") {
        mrcpp::CompFunction func_1;
        mrcpp::project(func_1, f, prec);

        SECTION("into non-shared function") {
            mrcpp::CompFunction func_2;
            mrcpp::multiply(func_2, func_1, func_1, -1.0);
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.getSquareNorm()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(0.0));
        }
        SECTION("into non-shared function conjugate") {
            mrcpp::CompFunction func_2;
            mrcpp::multiply(func_2, func_1, func_1, -1.0, false, false, true);
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.getSquareNorm()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(0.0));
        }
#ifdef MRCHEM_HAS_MPI
        SECTION("into shared function") {
            mrcpp::CompFunction func_2(0, true);
            mrcpp::multiply(func_2, func_1, func_1, -1.0, false, false, true);
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.getSquareNorm()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(0.0));
        }
#endif
    }

#ifdef MRCHEM_HAS_MPI
    SECTION("multiply shared function") {
        mrcpp::CompFunction func_1(0, true);
        mrcpp::project(func_1, f, prec);

        SECTION("into non-shared function") {
            mrcpp::CompFunction func_2(0, false);
            ;
            mrcpp::multiply(func_2, func_1, func_1, -1.0, true);
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.getSquareNorm()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(0.0));
        }
        SECTION("into shared function") {
            mrcpp::CompFunction func_2(0, true);
            mrcpp::multiply(func_2, func_1, func_1, -1.0, true);
            REQUIRE(func_2.integrate().real() == Catch::Approx(func_1.getSquareNorm()));
            REQUIRE(func_2.integrate().imag() == Catch::Approx(0.0));
        }
    }
#endif
}

} // namespace qmfunction_tests
