#include "xc_backend.h"

#include <XCFun/xcfun.h>
#include <memory>
#include <stdexcept>
#include <utility>

namespace mrdft {

namespace {
struct XCFunDeleter {
    void operator()(xcfun_t *p) const { if (p) xcfun_delete(p); }
};
} // namespace

class XCFunBackend final : public XCBackend {
public:
    XCFunBackend()
        : xcfun_(xcfun_new(), XCFunDeleter{}) {}

    void set_functional(const std::string &name, double coeff) override {
        xcfun_set(xcfun_.get(), name.c_str(), coeff);
    }

    void configure(int order, bool spin, bool use_gamma) override {
        order_ = order;
        spin_  = spin;

        const bool gga = xcfun_is_gga(xcfun_.get());
        unsigned int mode         = 1;                 // partial derivatives
        unsigned int func_type    = gga ? 1u : 0u;     // 0=LDA, 1=GGA
        unsigned int dens_type    = static_cast<unsigned>(1 + (spin ? 1 : 0)); // 1 or 2
        unsigned int laplacian    = 0;
        unsigned int kinetic      = 0;
        unsigned int current      = 0;
        unsigned int exp_deriv    = static_cast<unsigned>(!use_gamma);
        if (!gga) exp_deriv = 0; // fall back to gamma-type for LDA

        xcfun_user_eval_setup(xcfun_.get(), order_, func_type, dens_type,
                              mode, laplacian, kinetic, current, exp_deriv);
    }

    bool is_gga() const override     { return xcfun_is_gga(xcfun_.get()); }
    bool is_metagga() const override { return xcfun_is_metagga(xcfun_.get()); }

    double amount_exx() const override {
        double exx = 0.0;
        xcfun_get(xcfun_.get(), "exx", &exx);
        return exx;
    }

    int input_length()  const override { return xcfun_input_length(xcfun_.get()); }
    int output_length() const override { return xcfun_output_length(xcfun_.get()); }

    Eigen::MatrixXd eval_transposed(const Eigen::MatrixXd &inp,
                                    double cutoff,
                                    bool is_spin_sep) const override {
        const int nInp = input_length();
        const int nOut = output_length();
        const int nPts = static_cast<int>(inp.rows());
        if (inp.cols() != nInp) throw std::runtime_error("XCFunBackend: invalid input shape");

        Eigen::MatrixXd out(nPts, nOut);
        out.setZero();

        Eigen::VectorXd in_row(nInp);
        Eigen::VectorXd out_row(nOut);
        for (int i = 0; i < nPts; ++i) {
            bool calc = true;
            if (is_spin_sep) {
                if (inp(i, 0) < cutoff && inp(i, 1) < cutoff) calc = false;
            } else {
                if (inp(i, 0) < cutoff) calc = false;
            }
            for (int j = 0; j < nInp; ++j) in_row(j) = inp(i, j);
            if (calc) xcfun_eval(xcfun_.get(), in_row.data(), out_row.data());
            for (int j = 0; j < nOut; ++j) out(i, j) = out_row(j);
        }
        return out;
    }

private:
    std::unique_ptr<xcfun_t, XCFunDeleter> xcfun_;
    int  order_ = 1;
    bool spin_  = false;
};

XCBackend_p make_xcfun_backend() {
    return std::make_shared<XCFunBackend>();
}

} // namespace mrdft
