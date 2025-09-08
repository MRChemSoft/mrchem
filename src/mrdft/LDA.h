#pragma once

#include <XCFun/xcfun.h>
#include "Functional.h"

namespace mrdft {

class LDA : public Functional {
public:
    LDA(int k, XC_p &f);
    ~LDA() override = default;

    bool isSpin() const override { return false; }

private:
    mrcpp::FunctionTreeVector<3> rho;

    int getCtrInputLength() const override;
    int getCtrOutputLength() const override { return 2; }

    void clear() override;
    mrcpp::FunctionTreeVector<3> setupXCInput() override;
    mrcpp::FunctionTreeVector<3> setupCtrInput() override;

    void preprocess(mrcpp::FunctionTreeVector<3> &inp) override;
    mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) override;
};

} // namespace mrdft
