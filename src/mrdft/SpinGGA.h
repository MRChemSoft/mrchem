#pragma once

#include <XCFun/xcfun.h>
#include "Functional.h"

namespace mrdft {

class SpinGGA : public Functional {
public:
    SpinGGA(int k, XC_p &f, std::unique_ptr<mrcpp::DerivativeOperator<3>> &d);
    ~SpinGGA() override = default;

    bool isSpin() const override { return true; }

private:
    std::unique_ptr<mrcpp::DerivativeOperator<3>> derivative{nullptr};
    mrcpp::FunctionTreeVector<3> rho_a;
    mrcpp::FunctionTreeVector<3> rho_b;
    mrcpp::FunctionTreeVector<3> grad_a;
    mrcpp::FunctionTreeVector<3> grad_b;

    int getCtrInputLength() const override;
    int getCtrOutputLength() const override { return 9; }

    void clear() override;
    mrcpp::FunctionTreeVector<3> setupXCInput() override;
    mrcpp::FunctionTreeVector<3> setupCtrInput() override;

    void preprocess(mrcpp::FunctionTreeVector<3> &inp) override;
    mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) override;
};

} // namespace mrdft
