//#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "mrchem.h"

#include "MRCPP/MWFunctions"

mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

void initialize_mra();
void finalize_mra();

int main(int argc, char *argv[]) {
    // global setup
    initialize_mra();

    // run tests
    int result = Catch::Session().run(argc, argv);

    // global cleanup
    finalize_mra();

    return (result < 0xff ? result : 0xff);
}

void initialize_mra() {
    // Defining resolution levels 2^{-n}
    int min_scale = -4;
    int max_depth = 25;

    // Constructing world box
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    mrcpp::BoundingBox<3> world(min_scale, corner, boxes);

    // Constructing basis
    int order = 5;
    mrcpp::InterpolatingBasis basis(order);

    // Constructing global MRA
    mrchem::MRA = new mrcpp::MultiResolutionAnalysis<3>(world, basis, max_depth);
}

void finalize_mra() {
    if (mrchem::MRA != 0) delete mrchem::MRA;
    mrchem::MRA = 0;
}

