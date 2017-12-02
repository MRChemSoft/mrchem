#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "TreeCalculator.h"

template<int D>
class ProjectionCalculator : public TreeCalculator<D> {
public:
    ProjectionCalculator(const RepresentableFunction<D> &inp_func) : func(&inp_func) { }
    virtual ~ProjectionCalculator() { }

protected:
    const RepresentableFunction<D> *func;

    virtual void calcNode(MWNode<D> &node);
};

