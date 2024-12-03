#include "surface_forces/xcStress.h"

#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "mrdft/MRDFT.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/density_utils.h"
#include "qmoperators/one_electron/NablaOperator.h"

using namespace Eigen;
using namespace mrchem;
using namespace std;

namespace surface_force{

/**
 * @brief Compute the exchange-correlation stress tensor for LDA functional
 * 
 * @param mrdft_p MRDFT object
 * @param rhoGrid MatrixXd with density values, shape (nGrid, 1)
 * @return std::vector<Eigen::Matrix3d> vector of 3x3 matrices with stress tensor at each grid point
 */
std::vector<Eigen::Matrix3d> xcLDAStress(unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid){
    int nGrid = rhoGrid.rows();
    std::vector<Eigen::Matrix3d> out(nGrid);
    Eigen::MatrixXd xcOUT =  mrdft_p->functional().evaluate_transposed(rhoGrid);
    for (int i = 0; i < nGrid; i++) {
        out[i] = Matrix3d::Zero();
        for (int j = 0; j < 3; j++) {
            out[i](j, j) = xcOUT(i, 0) - xcOUT(i, 1) * rhoGrid(i);
        }
    }
    return out;
}

/**
 * @brief Compute the exchange-correlation stress tensor for LDA functional for open shell systems
 * 
 * @param mrdft_p MRDFT object
 * @param rhoGridAlpha MatrixXd with alpha density values, shape (nGrid, 1)
 * @param rhoGridBeta MatrixXd with beta density values, shape (nGrid, 1)
 * @return std::vector<Eigen::Matrix3d> vector of 3x3 matrices with stress tensor at each grid point
 */
std::vector<Matrix3d> xcLDASpinStress(unique_ptr<mrdft::MRDFT> &mrdft_p, MatrixXd &rhoGridAlpha, MatrixXd &rhoGridBeta){
    int nGrid = rhoGridAlpha.rows();
    Eigen::MatrixXd inp(rhoGridAlpha.rows(), 2);
    std::vector<Matrix3d> out = std::vector<Eigen::Matrix3d>(nGrid);
    inp.col(0) = rhoGridAlpha.col(0);
    inp.col(1) = rhoGridBeta.col(0);
    Eigen::MatrixXd xc = mrdft_p->functional().evaluate_transposed(inp);
    for (int i = 0; i < rhoGridAlpha.rows(); i++) {
        out[i] = Matrix3d::Zero();
        for (int j = 0; j < 3; j++) {
            out[i](j, j) = xc(i, 0) - xc(i, 1) * rhoGridAlpha(i) - xc(i, 2) * rhoGridBeta(i);
        }
    }
    return out;
}

/**
 * @brief Compute the exchange-correlation stress tensor for GGA functional
 * 
 * @param mrdft_p MRDFT object
 * @param rhoGrid MatrixXd with density values, shape (nGrid, 1)
 * @param nablaRhoGrid MatrixXd with gradient of density values, shape (nGrid, 3)
 * @return std::vector<Eigen::Matrix3d> vector of 3x3 matrices with stress tensor at each grid point
 */
std::vector<Matrix3d> xcGGAStress(unique_ptr<mrdft::MRDFT> &mrdft_p, mrcpp::FunctionTreeVector<3> &xc_pots, MatrixXd &rhoGrid, MatrixXd &nablaRhoGrid, Eigen::MatrixXd &gridPos){
    int nGrid = rhoGrid.rows();
    Eigen::MatrixXd inp(rhoGrid.rows(), 4);
    inp.col(0) = rhoGrid.col(0);
    inp.col(1) = nablaRhoGrid.col(0);
    inp.col(2) = nablaRhoGrid.col(1);
    inp.col(3) = nablaRhoGrid.col(2);
    Eigen::MatrixXd xcOUT =  mrdft_p->functional().evaluate_transposed(inp);
    std::vector<Matrix3d> out(nGrid);
    std::array<double, 3> pos;
    for (int i = 0; i < rhoGrid.rows(); i++) {
        out[i] = Matrix3d::Zero();
        pos[0] = gridPos(i, 0);
        pos[1] = gridPos(i, 1);
        pos[2] = gridPos(i, 2);
        for (int j = 0; j < 3; j++) {
            out[i](j, j) = xcOUT(i, 0) - rhoGrid(i) * std::get<1>(xc_pots[0])->evalf(pos);
        }
        for (int j1 = 0; j1 < 3; j1++) {
            for (int j2 = 0; j2 < 3; j2++) {
                out[i](j1, j2) = out[i](j1, j2) - xcOUT(i, 2 + j1) * nablaRhoGrid(i, j2);
            }
        }
    }
    return out;
}

/**
 * @brief Compute the exchange-correlation stress tensor for GGA functional for open shell systems
 * 
 * @param mrdft_p MRDFT object
 * @param rhoGridAlpha MatrixXd with alpha density values, shape (nGrid, 1)
 * @param rhoGridBeta MatrixXd with beta density values, shape (nGrid, 1)
 * @param nablaRhoGridAlpha MatrixXd with gradient of alpha density values, shape (nGrid, 3)
 * @param nablaRhoGridBeta MatrixXd with gradient of beta density values, shape (nGrid, 3)
 * @return std::vector<Eigen::Matrix3d> vector of 3x3 matrices with stress tensor at each grid point
 */
std::vector<Matrix3d> xcGGASpinStress(unique_ptr<mrdft::MRDFT> &mrdft_p, mrcpp::FunctionTreeVector<3> &xc_pots, MatrixXd &rhoGridAlpha, MatrixXd &rhoGridBeta, MatrixXd &nablaRhoGridAlpha, MatrixXd &nablaRhoGridBeta, Eigen::MatrixXd &gridPos){
    int nGrid = rhoGridAlpha.rows();
    Eigen::MatrixXd inp(rhoGridAlpha.rows(), 8);
    std::vector<Matrix3d> out = std::vector<Eigen::Matrix3d>(nGrid);
    inp.col(0) = rhoGridAlpha.col(0);
    inp.col(1) = rhoGridBeta.col(0);
    inp.col(2) = nablaRhoGridAlpha.col(0);
    inp.col(3) = nablaRhoGridAlpha.col(1);
    inp.col(4) = nablaRhoGridAlpha.col(2);
    inp.col(5) = nablaRhoGridBeta.col(0);
    inp.col(6) = nablaRhoGridBeta.col(1);
    inp.col(7) = nablaRhoGridBeta.col(2);
    Eigen::MatrixXd xc = mrdft_p->functional().evaluate_transposed(inp);
    std::array<double, 3> pos;
    for (int i = 0; i < rhoGridAlpha.rows(); i++) {
        out[i] = Matrix3d::Zero();
        pos[0] = gridPos(i, 0);
        pos[1] = gridPos(i, 1);
        pos[2] = gridPos(i, 2);
        for (int j = 0; j < 3; j++) {
            out[i](j, j) = xc(i, 0) - std::get<1>(xc_pots[0])->evalf(pos) * rhoGridAlpha(i) - std::get<1>(xc_pots[1])->evalf(pos) * rhoGridBeta(i);
        }
        for (int j1 = 0; j1 < 3; j1++) {
            for (int j2 = 0; j2 < 3; j2++) {
                out[i](j1, j2) = out[i](j1, j2) 
                    - xc(i, 3 + j1) * nablaRhoGridAlpha(i, j2) - xc(i, 6 + j1) * nablaRhoGridBeta(i, j2);
            }
        }
    }
    return out;
}

/**
 * @brief Compute the exchange-correlation stress tensor on a grid
 * 
 * @param mrdft_p MRDFT object
 * @param phi OrbitalVector
 * @param nabla NablaOperator (must be set up prior to calling this function)
 * @param gridPos MatrixXd with grid positions, shape (nGrid, 3)
 * @param isOpenShell bool, true if open shell calculation
 * @param prec precision to use in density representation
 */
std::vector<Eigen::Matrix3d> getXCStress(unique_ptr<mrdft::MRDFT> &mrdft_p, mrcpp::FunctionTreeVector<3> &xc_pots, std::shared_ptr<OrbitalVector> phi,
        std::shared_ptr<NablaOperator> nabla, MatrixXd &gridPos, bool isOpenShell, double prec){

    bool isGGA = mrdft_p->functional().isGGA();
    bool isHybrid = mrdft_p->functional().isHybrid();
    if (isHybrid) {
        MSG_ABORT("Exact exchange is not implemented for forces computed with surface integrals");
    }

    std::array<double, 3> pos;
    int nGrid = gridPos.rows();

    vector<Matrix3d> xcStress;

    if (isOpenShell) {
        MatrixXd rhoGridAlpha(nGrid, 1);
        MatrixXd rhoGridBeta(nGrid, 1);
        mrchem::Density rhoA(false);
        mrchem::Density rhoB(false);
        mrchem::density::compute(prec, rhoA, *phi, DensityType::Alpha);
        mrchem::density::compute(prec, rhoB, *phi, DensityType::Beta);

        for (int i = 0; i < nGrid; i++) { // compute density on grid
            pos[0] = gridPos(i, 0);
            pos[1] = gridPos(i, 1);
            pos[2] = gridPos(i, 2);
            rhoGridAlpha(i) = rhoA.real().evalf(pos);
            rhoGridBeta(i) = rhoB.real().evalf(pos);
        }

        if (isGGA) {
            mrchem::NablaOperator nablaOP = *nabla;
            std::vector<mrchem::Orbital> nablaRhoAlpha = nablaOP(rhoA);
            std::vector<mrchem::Orbital> nablaRhoBeta = nablaOP(rhoB);
            MatrixXd nablaRhoGridAlpha(nGrid, 3);
            MatrixXd nablaRhoGridBeta(nGrid, 3);
            for (int i = 0; i < nGrid; i++) {
                pos[0] = gridPos(i, 0);
                pos[1] = gridPos(i, 1);
                pos[2] = gridPos(i, 2);
                nablaRhoGridAlpha(i, 0) = nablaRhoAlpha[0].real().evalf(pos);
                nablaRhoGridAlpha(i, 1) = nablaRhoAlpha[1].real().evalf(pos);
                nablaRhoGridAlpha(i, 2) = nablaRhoAlpha[2].real().evalf(pos);
                nablaRhoGridBeta(i, 0) = nablaRhoBeta[0].real().evalf(pos);
                nablaRhoGridBeta(i, 1) = nablaRhoBeta[1].real().evalf(pos);
                nablaRhoGridBeta(i, 2) = nablaRhoBeta[2].real().evalf(pos);
            }

            xcStress = xcGGASpinStress(mrdft_p, xc_pots, rhoGridAlpha, rhoGridBeta, nablaRhoGridAlpha, nablaRhoGridBeta, gridPos);
        } else {
            xcStress = xcLDASpinStress(mrdft_p, rhoGridAlpha, rhoGridBeta);
        }

    } else { // closed shell
        MatrixXd rhoGrid(nGrid, 1);
        mrchem::Density rho(false);
        mrchem::density::compute(prec, rho, *phi, DensityType::Total);

        for (int i = 0; i < nGrid; i++) { // compute density on grid
            pos[0] = gridPos(i, 0);
            pos[1] = gridPos(i, 1);
            pos[2] = gridPos(i, 2);
            rhoGrid(i) = rho.real().evalf(pos);
        }

        if (isGGA) {
            mrchem::NablaOperator nablaOP = *nabla;
            std::vector<mrchem::Orbital> nablaRho = nablaOP(rho);
            MatrixXd nablaRhoGrid(nGrid, 3);
            for (int i = 0; i < nGrid; i++) {
                pos[0] = gridPos(i, 0);
                pos[1] = gridPos(i, 1);
                pos[2] = gridPos(i, 2);
                nablaRhoGrid(i, 0) = nablaRho[0].real().evalf(pos);
                nablaRhoGrid(i, 1) = nablaRho[1].real().evalf(pos);
                nablaRhoGrid(i, 2) = nablaRho[2].real().evalf(pos);
            }
            xcStress = xcGGAStress(mrdft_p, xc_pots, rhoGrid, nablaRhoGrid, gridPos);
        } else {
            xcStress = xcLDAStress(mrdft_p, rhoGrid);
        }
    }
    return xcStress;
}
} // namespace surface_force