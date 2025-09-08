#pragma once

#include <XCFun/xcfun.h>
#include "Functional.h"

namespace mrdft {

class GGA : public Functional {
public:
    GGA(int k, XC_p &f, std::unique_ptr<mrcpp::DerivativeOperator<3>> &d);
    ~GGA() override = default;

    bool isSpin() const override { return false; }

private:
    std::unique_ptr<mrcpp::DerivativeOperator<3>> derivative{nullptr};
    mrcpp::FunctionTreeVector<3> rho;
    mrcpp::FunctionTreeVector<3> grad;

    int getCtrInputLength() const override;
    int getCtrOutputLength() const override { return 5; }

    void clear() override;
    mrcpp::FunctionTreeVector<3> setupXCInput() override;
    mrcpp::FunctionTreeVector<3> setupCtrInput() override;

    void preprocess(mrcpp::FunctionTreeVector<3> &inp) override;
    mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) override;
};

} // namespace mrdft
