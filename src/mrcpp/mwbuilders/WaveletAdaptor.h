#pragma once

#include "TreeAdaptor.h"
#include "constants.h"

template<int D>
class WaveletAdaptor : public TreeAdaptor<D> {
public:
    WaveletAdaptor(double pr,
		   int ms,
		   bool ap = false,
		   double sf = 1.0)
            : TreeAdaptor<D>(ms),
              absPrec(ap),
              prec(pr),
	      splitFac(sf) { }
    virtual ~WaveletAdaptor() { }

    void setPrecision(double pr) { this->prec = pr; }
    void setAbsPrec(bool ap) { this->absPrec = ap; }

protected:
    bool absPrec;
    double prec;
    double splitFac;

    virtual bool splitNode(const MWNode<D> &node) const {
        return node.splitCheck(this->prec, this->splitFac, this->absPrec);
    }
};

