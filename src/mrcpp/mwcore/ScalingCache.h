#pragma once

#include "ObjectCache.h"

#define getLegendreScalingCache(X)\
    ScalingCache<LegendreBasis> &X = \
    ScalingCache<LegendreBasis>::getInstance()
#define getInterpolatingScalingCache(X)\
    ScalingCache<InterpolatingBasis> &X = \
    ScalingCache<InterpolatingBasis>::getInstance()

template<class P>
class ScalingCache: public ObjectCache<P> {
public:
    static ScalingCache &getInstance() {
        static ScalingCache theScalingCache;
        return theScalingCache;
    }
    virtual void load(int order) {
        SET_CACHE_LOCK();
        if (not this->hasId(order)) {
            P *f = new P(order);
            int memo = 2 * SQUARE(order+1) * sizeof(double); //approx
            ObjectCache<P>::load(order, f, memo);
        }
        UNSET_CACHE_LOCK();
    }

    P &get(int order) {
        if (not this->hasId(order)) {
            load(order);
        }
        return ObjectCache<P>::get(order);
    }
private:
    ScalingCache() { }
    virtual ~ScalingCache() { }
    ScalingCache(const ScalingCache<P> &sc);
    ScalingCache<P> &operator=(const ScalingCache<P> &sc) { return *this;	}
};

