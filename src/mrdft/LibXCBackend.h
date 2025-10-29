#pragma once

#include <vector>
#include <memory>
#include <stdexcept>
#include <Eigen/Core>

#include <xc.h>  // LibXC C header

#include "LDA.h"
#include "GGA.h"
#include "SpinLDA.h"
#include "SpinGGA.h"

namespace mrdft {

/* -------- Unpolarized LDA -------- */
class LibXCLDA : public LDA {
public:
    LibXCLDA(int k, XC_p &f, const std::vector<int>& ids)
        : LDA(k, f), ids_(ids), nspin_(XC_UNPOLARIZED) {}

protected:
    // Hook override
    Eigen::MatrixXd eval_lda_transposed(Eigen::MatrixXd &inp) const override;

private:
    std::vector<int> ids_;
    int nspin_;
};

/* -------- Unpolarized GGA -------- */
class LibXCGGA : public GGA {
public:
    LibXCGGA(int k, XC_p &f, std::unique_ptr<mrcpp::DerivativeOperator<3>> &d, const std::vector<int>& ids)
        : GGA(k, f, d), ids_(ids), nspin_(XC_UNPOLARIZED) {}

protected:
    // Hook override
    Eigen::MatrixXd eval_gga_transposed(Eigen::MatrixXd &inp) const override;

private:
    std::vector<int> ids_;
    int nspin_;
};

/* -------- Spin LDA -------- */
class LibXCSpinLDA : public SpinLDA {
public:
    LibXCSpinLDA(int k, XC_p &f, const std::vector<int>& ids)
        : SpinLDA(k, f), ids_(ids), nspin_(XC_POLARIZED) {}

protected:
    // Hook override
    Eigen::MatrixXd eval_lda_transposed(Eigen::MatrixXd &inp) const override;

private:
    std::vector<int> ids_;
    int nspin_;
};

/* -------- Spin GGA -------- */
class LibXCSpinGGA : public SpinGGA {
public:
    LibXCSpinGGA(int k, XC_p &f, std::unique_ptr<mrcpp::DerivativeOperator<3>> &d, const std::vector<int>& ids)
        : SpinGGA(k, f, d), ids_(ids), nspin_(XC_POLARIZED) {}

protected:
    // Hook override
    Eigen::MatrixXd eval_gga_transposed(Eigen::MatrixXd &inp) const override;

private:
    std::vector<int> ids_;
    int nspin_;
};

} // namespace mrdft
