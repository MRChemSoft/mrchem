#pragma once

#include "WaveletAdaptor.h"

class OperatorAdaptor : public WaveletAdaptor<2> {
public:
    OperatorAdaptor(double pr, int ms, bool ap = false)
            : WaveletAdaptor<2>(pr, ms, ap) { }
    virtual ~OperatorAdaptor() { }

protected:
    virtual bool splitNode(const MWNode<2> &node) const {
        int chkCompNorm = 0;
        for (int i = 1; i < 4; i++) {
            if (node.getComponentNorm(i) > 0.0) {
                chkCompNorm = 1;
            }
        }

        const int *l = node.getTranslation();
        int chkTransl = (l[0] == 0 or l[1] == 0);

        int split = 1;
        split *= chkTransl;
        split *= chkCompNorm;

        return split;
    }
};

