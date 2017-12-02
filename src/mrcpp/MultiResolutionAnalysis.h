#pragma once

#include "BoundingBox.h"
#include "InterpolatingBasis.h"
#include "LegendreBasis.h"
#include "FilterCache.h"

template<int D>
class MultiResolutionAnalysis {
public:
    MultiResolutionAnalysis(const MultiResolutionAnalysis<D> &mra)
            : maxDepth(mra.maxDepth),
              basis(mra.basis),
              world(mra.world) {
        if (getMaxDepth() > MaxDepth) MSG_FATAL("Beyond MaxDepth");
        if (getMaxScale() > MaxScale) MSG_FATAL("Beyond MaxScale");
        setupFilter();
    }
    MultiResolutionAnalysis(const BoundingBox<D> &bb,
                            const ScalingBasis &sb,
                            int depth = MaxDepth)
            : maxDepth(depth),
              basis(sb),
              world(bb) {
        if (getMaxDepth() > MaxDepth) MSG_FATAL("Beyond MaxDepth");
        if (getMaxScale() > MaxScale) MSG_FATAL("Beyond MaxScale");
        setupFilter();
    }
    virtual ~MultiResolutionAnalysis() { }

    int getOrder() const { return this->basis.getScalingOrder(); }
    int getMaxDepth() const { return this->maxDepth; }
    int getMaxScale() const { return this->world.getScale() + this->maxDepth; }

    const MWFilter &getFilter() const { return *this->filter; }
    const ScalingBasis &getScalingBasis() const { return this->basis; }
    const BoundingBox<D> &getWorldBox() const { return this->world; }

    MultiResolutionAnalysis<1> getKernelMRA() const {
        const BoundingBox<D> &box = getWorldBox();
        const ScalingBasis &basis = getScalingBasis();

        int type = basis.getScalingType();
        int kern_order = 2*basis.getScalingOrder() + 1;

        ScalingBasis *kern_basis = 0;
        if (type == Interpol) {
            kern_basis = new InterpolatingBasis(kern_order);
        } else if (type == Legendre) {
            kern_basis = new LegendreBasis(kern_order);
        } else {
            MSG_FATAL("Invalid scaling type");
        }

        int max_l = 0;
        for (int i = 0; i < D; i++) {
            if (box.size(i) > max_l) {
                max_l = box.size(i);
            }
        }
        int start_l = -max_l;
        int tot_l = 2*max_l;
        NodeIndex<1> idx(box.getScale(), &start_l);
        BoundingBox<1> kern_box(idx, &tot_l);
        MultiResolutionAnalysis<1> mra(kern_box, *kern_basis);
        delete kern_basis;
        return mra;
    }

    MultiResolutionAnalysis<2> getOperatorMRA() const {
        const BoundingBox<D> &box = getWorldBox();
        const ScalingBasis &basis = getScalingBasis();

        int maxn = 0;
        for (int i = 0; i < D; i++) {
            if (box.size(i) > maxn) {
                maxn = box.size(i);
            }
        }
        int nbox[2] = { maxn, maxn};
        NodeIndex<2> idx(box.getScale());
        BoundingBox<2> oper_box(idx, nbox);
        return MultiResolutionAnalysis<2>(oper_box, basis);
    }

    bool operator==(const MultiResolutionAnalysis<D> &mra) const {
        if (this->basis != mra.basis) return false;
        if (this->world != mra.world) return false;
        if (this->maxDepth != mra.maxDepth) return false;
        return true;
    }
    bool operator!=(const MultiResolutionAnalysis<D> &mra) const {
        if (this->basis != mra.basis) return true;
        if (this->world != mra.world) return true;
        if (this->maxDepth != mra.maxDepth) return true;
        return false;
    }

    void print() const {
        println(0, std::endl);
        println(0, "============================================================");
        println(0, "                  MultiResolution Analysis                  ");
        println(0, "------------------------------------------------------------");
        println(0, this->basis);
        println(0, "------------------------------------------------------------");
        println(0, this->world);
        println(0, "============================================================");
        println(0, std::endl);
    }

protected:
    const int maxDepth;
    const ScalingBasis basis;
    const BoundingBox<D> world;
    MWFilter *filter;

    void setupFilter() {
        getLegendreFilterCache(lfilters);
        getInterpolatingFilterCache(ifilters);
        int k = this->basis.getScalingOrder();
        int type = this->basis.getScalingType();
        switch (type) {
        case Legendre:
            this->filter = &lfilters.get(k);
            break;
        case Interpol:
            this->filter = &ifilters.get(k);
            break;
        default:
            MSG_ERROR("Invalid scaling basis selected.")
        }
    }
};

