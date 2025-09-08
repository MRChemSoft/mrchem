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

/* Minimal RAII wrapper for LibXC function handles */
class LibXCHandle {
public:
    LibXCHandle(const std::vector<int>& ids, int nspin)
        : funcs(ids.size()), nspin(nspin) {
        for (size_t k = 0; k < ids.size(); ++k) {
            if (xc_func_init(&funcs[k], ids[k], nspin) != 0) {
                throw std::runtime_error("LibXC: could not initialize functional id=" + std::to_string(ids[k]));
            }
        }
    }
    ~LibXCHandle() {
        for (auto &f : funcs) xc_func_end(&f);
    }
    std::vector<xc_func_type> funcs;
    int nspin{XC_UNPOLARIZED};
};

/* -------- Unpolarized LDA -------- */
class LibXCLDA : public LDA {
public:
    LibXCLDA(int k, XC_p &f, const std::vector<int>& ids)
        : LDA(k, f), handle(ids, XC_UNPOLARIZED) {}

    Eigen::MatrixXd evaluate_transposed(Eigen::MatrixXd &inp) const override;

private:
    LibXCHandle handle;
};

/* -------- Unpolarized GGA -------- */
class LibXCGGA : public GGA {
public:
    LibXCGGA(int k, XC_p &f, std::unique_ptr<mrcpp::DerivativeOperator<3>> &d, const std::vector<int>& ids)
        : GGA(k, f, d), handle(ids, XC_UNPOLARIZED) {}

    Eigen::MatrixXd evaluate_transposed(Eigen::MatrixXd &inp) const override;

private:
    LibXCHandle handle;
};

/* -------- Spin LDA -------- */
class LibXCSpinLDA : public SpinLDA {
public:
    LibXCSpinLDA(int k, XC_p &f, const std::vector<int>& ids)
        : SpinLDA(k, f), handle(ids, XC_POLARIZED) {}

    Eigen::MatrixXd evaluate_transposed(Eigen::MatrixXd &inp) const override;

private:
    LibXCHandle handle;
};

/* -------- Spin GGA -------- */
class LibXCSpinGGA : public SpinGGA {
public:
    LibXCSpinGGA(int k, XC_p &f, std::unique_ptr<mrcpp::DerivativeOperator<3>> &d, const std::vector<int>& ids)
        : SpinGGA(k, f, d), handle(ids, XC_POLARIZED) {}

    Eigen::MatrixXd evaluate_transposed(Eigen::MatrixXd &inp) const override;

private:
    LibXCHandle handle;
};

} // namespace mrdft
