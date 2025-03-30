#include "catch.hpp"
#include "pseudopotential/sphericalHarmonics.h"
#include "qmfunctions/qmfunction_utils.h"

#include "pseudopotential/projector.h"



using namespace mrchem;

namespace density_tests {

TEST_CASE("sperical_harmonics", "[spherical_harmonics]") {
    int lmax = 4;
    const double prec = 1.0e-5;
    // const double thrs = 1.0e-12;

    SECTION("calc spherical harmonics") {
        double alpha = .5;
        std::vector<ProjectorFunction> pps;
        int iOrb = 0;
        for (int l = 0; l < lmax; l++) {
            for (int m = -l; m <= l; m++) {
                for (int idim = 0; idim < 1; idim++) {
                    mrcpp::Coord<3> r;
                    r[0] = 0.0;
                    r[1] = 0.0;
                    r[2] = 0.0;
                    double rl = 2.117559e-01;
                    // ProjectorFunction pp(r, rl, idim, l, m, prec);
                    pps.push_back(ProjectorFunction(r, rl , idim, l, m, prec));
                }
            }
        }

        Eigen::MatrixXd S(pps.size(), pps.size());
        for (int i = 0; i < pps.size(); i++) {
            for (int j = 0; j < pps.size(); j++) {
                // S(i, j) = qmfunction::dot(Phi[i], Phi[j]).real();
                S(i, j) = mrcpp::cplxfunc::dot(*pps[i].projector_ptr, *pps[j].projector_ptr).real();
                // std::cout << i << " " << j << " " << S(i, j) << " " << pps[i].projector_ptr->norm() << std::endl;
            }
        }
        Eigen::MatrixXd S_ref = Eigen::MatrixXd::Zero(pps.size(), pps.size());
        for (int i = 0; i < pps.size(); i++) {
            S_ref(i, i) = 1.0;
        }
        int iorb = 0;
        int jorb = 0;
        for (int l = 0; l < lmax; l++) {
            for (int m = -l; m <= l; m++) {
                jorb = 0;
                for (int ll = 0; ll < lmax; ll++) {
                    for (int mm = -ll; mm <= ll; mm++) {
                        if (l == ll && m == mm) {
                            // std::cout << "l " << l << " m " << m << " ll " << ll << " mm " << mm << std::endl;
                            // std::cout << iorb << " " << jorb << " " << S(iorb, jorb) << " " << S_ref(iorb, jorb) << std::endl;
                            REQUIRE(std::abs(S(iorb, jorb) - (S_ref(iorb, jorb))) < prec);
                        }
                        jorb++;
                    }
                }
                iorb++;
            }
        }
    }
}

} // namespace density_tests