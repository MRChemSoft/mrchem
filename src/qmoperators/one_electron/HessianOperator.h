#pragma once

#include "tensor/RankOneOperator.h"

#include "qmoperators/QMDerivative.h"
#include "qmoperators/one_electron/NablaOperator.h"

namespace mrchem {

class HessianOperator final : public RankOneOperator<6> {
public:
    HessianOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D1, std::shared_ptr<mrcpp::DerivativeOperator<3>> D2, double prec) {
        bool imag = false;
        NablaOperator N1 = NablaOperator(D1, imag);
        NablaOperator N2 = NablaOperator(D2, imag);
        N1.setup(prec);
        N2.setup(prec);

        // Invoke operator= to assign *this operator
        RankOneOperator<6> &d = (*this);
        d[0] = N2[0];
        d[1] = N2[1];
        d[2] = N2[2];
        d[3] = N1[0] * N1[1];
        d[4] = N1[0] * N1[2];
        d[5] = N1[1] * N1[2];

        d[0].name() = "del[x]del[x]";
        d[1].name() = "del[y]del[y]";
        d[2].name() = "del[z]del[z]";
        d[3].name() = "del[x]del[y]";
        d[4].name() = "del[x]del[z]";
        d[5].name() = "del[y]del[z]";
    }
};

}