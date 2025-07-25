#include "LDA.h"
#include "Factory.h"

#include <Eigen/Core>
#include <MRCPP/MWFunctions>
#include <MRCPP/MWOperators>
#include <MRCPP/trees/FunctionNode.h>
#include <MRCPP/Printer>
#include "xc_utils.h"

#include <xc_funcs.h>
#include <xc.h>

namespace mrdft {

LDA::LDA(int order, std::vector<mrlibxc::LibXCData> &data)
    : order_(order), functionals_(functionals_) {
    // Store for later use if needed
    this->order_ = order;
    this->functionals_ = functionals;

    // Set up masks (e.g., for LibXC)
    xc_mask = xc_utils::build_output_mask(true, false, this->order_);
    d_mask = xc_utils::build_density_mask(true, false, this->order_);
    xc_ = std::make_shared<LDAHandler>(data);
}



/** @brief Clear internal functions */
void LDA::clear() {
    mrcpp::clear(this->rho, false);
}

/** @brief Number of function involved in contraction step */
int LDA::getCtrInputLength() const {
    if (order_ < 2) return 0;
    if (order_ == 2) return 1;
    NOT_IMPLEMENTED_ABORT;
    return 0;
}

/** @brief Set up input for LibXC evaluation: just rho_0 */
mrcpp::FunctionTreeVector<3> LDA::setupXCInput() {
    if (rho.empty()) MSG_ERROR("Density not initialized");
    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(rho[0]);
    return out_vec;
}

/** @brief For LDA at order 2: include rho_1 in contraction input */
mrcpp::FunctionTreeVector<3> LDA::setupCtrInput() {
    mrcpp::FunctionTreeVector<3> out_vec;
    if (order_ == 2) {
        out_vec.push_back(rho[1]);
    } else if (order_ > 2) {
        NOT_IMPLEMENTED_ABORT;
    }
    return out_vec;
}

/** @brief Store rho from input */
void LDA::preprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    if ((int)inp_vec.size() != order_) MSG_ERROR("Invalid input length");
    if (!rho.empty()) MSG_ERROR("Density already initialized");

    for (auto i = 0; i < order_; i++)
        rho.push_back(inp_vec[i]);
}

/** @brief Evaluate LibXC and return energy and potential */
mrcpp::FunctionTreeVector<3> LDA::postprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    // Retrieve density function
    mrcpp::FunctionTree<3> &rho_ft = mrcpp::get_func(inp_vec, 0);
    const auto &grid = rho_ft.getGrid();
    const auto N = grid.size();

    std::vector<double> rho_vals(N), eps_vals(N), vxc_vals(N);

    Eigen::VectorXd eigen_vec = Eigen::Map<Eigen::VectorXd>(rho_vals.data(), rho_vals.size());
    rho_ft.setEndValues(eigen_vec);

    mrcpp::Coord<3> r = {x, y, z};
    double val = rho_ft.evalf(r);

    // Evaluate each functional in the list
    for (const auto &func_data : functionals_) {
        const auto &func = func_data.func;
        double weight = func_data.weight;

        std::vector<double> tmp_eps(N), tmp_vxc(N);
        xc_lda_exc_vxc(&func, N, rho_vals.data(), tmp_eps.data(), tmp_vxc.data());

        for (size_t i = 0; i < N; ++i) {
            eps_vals[i] += weight * tmp_eps[i];
            vxc_vals[i] += weight * tmp_vxc[i];
        }
    }

    // Create output FunctionTrees and load values
    mrcpp::FunctionTree<3> eps_ft(grid);
    mrcpp::FunctionTree<3> vxc_ft(grid);



    eps_ft.setEndValues(eps_vals);
    vxc_ft.setEndValues(vxc_vals);

    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(std::make_tuple(1.0, &eps_ft));
    out_vec.push_back(std::make_tuple(1.0, &vxc_ft));

    return out_vec;
}

}  // namespace mrlibxc
