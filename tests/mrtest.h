#pragma once

#include "MRCPP/MWFunctions"

using mrcpp::BoundingBox;
using mrcpp::InterpolatingBasis;
using mrcpp::MultiResolutionAnalysis;

namespace mrtest {

inline MultiResolutionAnalysis<3>* initialize_mra() {
    // Defining resolution levels 2^{-n}
    int min_scale = -4;
    int max_depth = 25;

    // Constructing world box
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    BoundingBox<3> world(min_scale, corner, boxes);

    // Constructing basis
    int order = 5;
    InterpolatingBasis basis(order);

    // Constructing global MRA
    return new MultiResolutionAnalysis<3>(world, basis, max_depth);
}

} // namespace mrtest
