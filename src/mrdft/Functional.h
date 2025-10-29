#pragma once

#include <memory>
#include <Eigen/Core>
#include <MRCPP/MWFunctions>
#include <MRCPP/MWOperators>
#include <MRCPP/trees/FunctionNode.h>
#include <XCFun/xcfun.h>

namespace mrdft {

using XC_p = std::unique_ptr<xcfun_t, decltype(&xcfun_delete)>;

class Functional {
public:
    Functional(int k, XC_p &f)
        : order(k)
        , xcfun(std::move(f)) {}
    virtual ~Functional() = default;

    void makepot(mrcpp::FunctionTreeVector<3> &inp,
                 std::vector<mrcpp::FunctionNode<3> *> xcNodes) const;

    void setLogGradient(bool log) { log_grad = log; }
    void setDensityCutoff(double cut) { cutoff = cut; }
    void setDerivOp(std::unique_ptr<mrcpp::DerivativeOperator<3>> &d) { derivOp = std::move(d); }

    // What kind of functional?
    virtual bool isSpin() const = 0;

    // Virtual + safe default implementations (so LibXC backends don’t poke XCFun).
    virtual bool isGGA() const {
        if (!xcfun) return false;
        return xcfun_is_gga(xcfun.get());
    }
    virtual bool isMetaGGA() const {
        if (!xcfun) return false;
        return xcfun_is_metagga(xcfun.get());
    }
    virtual bool isLDA() const { return !(isGGA() || isMetaGGA()); }
    virtual bool isHybrid() const { return (std::abs(amountEXX()) > 1.0e-10); }

    // Virtual so LibXC subclasses can override (e.g., return 0.0 for LDA/GGA).
    virtual double amountEXX() const {
        if (!xcfun) return 0.0;
        double exx = 0.0;
        xcfun_get(xcfun.get(), "exx", &exx);
        return exx;
    }

    double XCenergy = 0.0;

    Eigen::MatrixXd evaluate(Eigen::MatrixXd &inp) const;
    Eigen::MatrixXd evaluate_transposed(Eigen::MatrixXd &inp) const;

    friend class MRDFT;

protected:
    const int order;
    bool log_grad{false};
    double cutoff{-1.0};
    Eigen::VectorXi d_mask;
    Eigen::MatrixXi xc_mask;
    XC_p xcfun;
    std::unique_ptr<mrcpp::DerivativeOperator<3>> derivOp{nullptr};

    int getXCInputLength() const { return xcfun ? xcfun_input_length(xcfun.get()) : 0; }
    int getXCOutputLength() const { return xcfun ? xcfun_output_length(xcfun.get()) : 0; }
    virtual int getCtrInputLength() const = 0;
    virtual int getCtrOutputLength() const = 0;

    Eigen::MatrixXd contract(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const;
    Eigen::MatrixXd contract_transposed(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const;

    virtual void clear() = 0;
    virtual mrcpp::FunctionTreeVector<3> setupXCInput() = 0;
    virtual mrcpp::FunctionTreeVector<3> setupCtrInput() = 0;

    virtual void preprocess(mrcpp::FunctionTreeVector<3> &inp) = 0;
    virtual mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) = 0;

    // ---- Hook pattern: LibXC backends override these ----
    virtual Eigen::MatrixXd eval_lda_transposed(Eigen::MatrixXd &inp) const;
    virtual Eigen::MatrixXd eval_gga_transposed(Eigen::MatrixXd &inp) const;

    // Helper to put XCFun in a consistent user-eval mode
    static void ensure_xcfun_user_setup(xcfun_t* xf, int order, bool spin, bool is_gga);
};

} // namespace mrdft
