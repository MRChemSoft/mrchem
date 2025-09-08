/*
 * (banner unchanged)
 */

#include "MRCPP/MWFunctions"
#include "MRCPP/Printer"

#include "LDA.h"
#include "xc_utils.h"

namespace mrdft {

LDA::LDA(int k, XC_p &f)
        : Functional(k, f) {
    xc_mask = xc_utils::build_output_mask(true, false, this->order);
    d_mask = xc_utils::build_density_mask(true, false, this->order);
}

/** @brief Clear internal functions
 *
 * Ownership of densities is outside MRDFT -> clear
 * Ownership of gradients is inside MRDFT -> free
 */
void LDA::clear() {
    mrcpp::clear(this->rho, false);
}

/** @brief Number of function involved in contraction step */
int LDA::getCtrInputLength() const {
    int length = -1;
    if (this->order < 2) length = 0;
    if (this->order == 2) length = 1;
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    return length;
}

/** @brief Collect input functions to xcfun evaluation step
 *
 * For LDA : [rho_0]
 */
mrcpp::FunctionTreeVector<3> LDA::setupXCInput() {
    if (this->rho.size() < 1) MSG_ERROR("Density not initialized");
    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(this->rho[0]);
    return out_vec;
}

/** @brief Collect input functions to contraction step
 *
 * For LDA:
 * Ground State: No contraction, empty vector
 * Linear Response: [rho_1]
 * Higher Response: NOT_IMPLEMENTED
 */
mrcpp::FunctionTreeVector<3> LDA::setupCtrInput() {
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    mrcpp::FunctionTreeVector<3> out_vec;
    if (order == 2) out_vec.push_back(this->rho[1]);
    return out_vec;
}

/** @brief Prepare input functions to xcfun
 *
 * Collects input densities and computes necessary gradients.
 *
 * Ordering of input:
 * inp_vec[0] = alpha_0
 * inp_vec[1] = beta_0
 * inp_vec[2] = alpha_1
 * inp_vec[3] = beta_1
 * ...
 */
void LDA::preprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    if (inp_vec.size() != this->order) MSG_ERROR("Invalid input length");
    if (this->rho.size() > 0) MSG_ERROR("Density not empty");

    for (auto i = 0; i < this->order; i++) this->rho.push_back(inp_vec[i]);
}

/** @brief Compute final output functions
 *
 * Combine the raw partial derivatives from xcfun into functional derivatives.
 *
 * For LDA:
 * f_xc       : out[0] = inp[0]
 * df_xc/drho : out[1] = inp[1]
 */
mrcpp::FunctionTreeVector<3> LDA::postprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    // Energy density
    mrcpp::FunctionTree<3> &f_xc = mrcpp::get_func(inp_vec, 0);
    inp_vec[0] = std::make_tuple<double, mrcpp::FunctionTree<3> *>(1.0, nullptr);

    // XC potential
    mrcpp::FunctionTree<3> &v_xc = mrcpp::get_func(inp_vec, 1);
    inp_vec[1] = std::make_tuple<double, mrcpp::FunctionTree<3> *>(1.0, nullptr);

    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(std::make_tuple(1.0, &f_xc));
    out_vec.push_back(std::make_tuple(1.0, &v_xc));
    return out_vec;
}

} // namespace mrdft
