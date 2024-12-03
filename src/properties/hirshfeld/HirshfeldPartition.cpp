#include "HirshfeldPartition.h"
#include "utils/math_utils.h"

HirshfeldPartition::HirshfeldPartition(const mrchem::Molecule &mol, std::string data_dir) {

    this->nucs = std::make_shared<mrchem::Nuclei>(mol.getNuclei());
    this->nNucs = this->nucs->size();

    for (int i = 0; i < this->nNucs; i++) {
        std::string element = this->nucs->at(i).getElement().getSymbol();
        this->logDensities.push_back(HirshfeldRadInterpolater(element, data_dir));
        // Uncomment the following line to print the charge of the atomic density
        // std::cout << "Norm of " << element << " = " << this->logDensities[i].getNorm() << std::endl;
    }

}

double HirshfeldPartition::getHirshfeldPartitionIntegral(int index, mrcpp::ComplexFunction &rho, double prec) const {
    auto w_i_analytic = [index, this](const mrcpp::Coord<3> &r) {
        return this->evalf(r, index);
    };
    mrcpp::ComplexFunction w_i_MW;
    mrcpp::AnalyticFunction<3> w_i_analytic_func(w_i_analytic);
    mrcpp::cplxfunc::multiply(w_i_MW, rho.real(), w_i_analytic_func, prec);
    mrcpp::ComplexDouble c_complex = w_i_MW.integrate();
    double charge = c_complex.real();

    return charge;
}

double HirshfeldPartition::lseLogDens(const mrcpp::Coord<3> &r) const {
    Eigen::VectorXd lseLogDens_r(this->nNucs);
    double rr;
    mrcpp::Coord<3> nucPos;
    for (int i = 0; i < this->nNucs; i++) {
        nucPos = this->nucs->at(i).getCoord();
        rr = std::sqrt((r[0] - nucPos[0]) * (r[0] - nucPos[0])
            + (r[1] - nucPos[1]) * (r[1] - nucPos[1])
            + (r[2] - nucPos[2]) * (r[2] - nucPos[2]));
        lseLogDens_r(i) = this->logDensities[i].evalf(rr);
    }
    return mrchem::math_utils::logsumexp(lseLogDens_r);
}

double HirshfeldPartition::evalf(const mrcpp::Coord<3> &r, int iAt) const {
    mrcpp::Coord<3> nucPos = this->nucs->at(iAt).getCoord();
    double rr = std::sqrt((r[0] - nucPos[0]) * (r[0] - nucPos[0])
        + (r[1] - nucPos[1]) * (r[1] - nucPos[1])
        + (r[2] - nucPos[2]) * (r[2] - nucPos[2]));
    return std::exp(this->logDensities[iAt].evalf(rr) - this->lseLogDens(r));
}