#pragma once

#include <XCFun/xcfun.h>
#include "Functional.h"

namespace mrdft {

class SpinLDA : public Functional {
public:
    SpinLDA(int k, XC_p &f);
    ~SpinLDA() override = default;

    bool isSpin() const override { return true; }

private:
    mrcpp::FunctionTreeVector<3> rho_a;
    mrcpp::FunctionTreeVector<3> rho_b;

    int getCtrInputLength() const override;
    int getCtrOutputLength() const override { return 3; }

    void clear() override;
    mrcpp::FunctionTreeVector<3> setupXCInput() override;
    mrcpp::FunctionTreeVector<3> setupCtrInput() override;

    void preprocess(mrcpp::FunctionTreeVector<3> &inp) override;
    mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) override;
};

} // namespace mrdft
