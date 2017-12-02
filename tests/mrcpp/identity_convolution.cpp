#include "catch.hpp"

#include "factory_functions.h"
#include "MWProjector.h"
#include "IdentityConvolution.h"
#include "IdentityKernel.h"
#include "OperatorTree.h"
#include "MWConvolution.h"
#include "OperatorAdaptor.h"
#include "CrossCorrelationCalculator.h"
#include "BandWidth.h"

namespace identity_convolution {

template<int D> void applyIdentity();

TEST_CASE("Initialize identity convolution operator", "[init_identity], [identity_convolution], [mw_operator]") {
    double exp_prec  = 1.0e-6;
    double proj_prec = 1.0e-6;
    double ccc_prec  = 1.0e-4;

    const int n = -6;
    const int k = 5;

    SECTION("Initialize identity kernel") {
        IdentityKernel id_kern(exp_prec);
        REQUIRE( (id_kern.size() == 1) );

        SECTION("Project identity kernel") {
            int l = -1;
            int nbox = 2;
            NodeIndex<1> idx(n, &l);
            BoundingBox<1> box(idx, &nbox);

            InterpolatingBasis basis(2*k+1);
            MultiResolutionAnalysis<1> kern_mra(box, basis);
            GridGenerator<1> G;
            MWProjector<1> Q(proj_prec);

            FunctionTree<1> kern_tree(kern_mra);
            G(kern_tree, id_kern);
            Q(kern_tree, id_kern);
            REQUIRE( (kern_tree.integrate() == Approx(1.0).epsilon(proj_prec)) );

            SECTION("Build operator tree by cross correlation") {
                NodeIndex<2> idx(n);
                BoundingBox<2> box(idx);

                InterpolatingBasis basis(k);
                MultiResolutionAnalysis<2> oper_mra(box, basis);

                TreeBuilder<2> builder;
                OperatorAdaptor adaptor(ccc_prec, oper_mra.getMaxScale());
                CrossCorrelationCalculator calculator(kern_tree);

                OperatorTree oper_tree(oper_mra, ccc_prec);
                builder.build(oper_tree, calculator, adaptor, -1);
                oper_tree.setupOperNodeCache();

                oper_tree.calcBandWidth(1.0);
                BandWidth bw_1 = oper_tree.getBandWidth();
                oper_tree.clearBandWidth();

                oper_tree.calcBandWidth(0.001);
                BandWidth bw_2 = oper_tree.getBandWidth();
                oper_tree.clearBandWidth();

                oper_tree.calcBandWidth(-1.0);
                BandWidth bw_3 = oper_tree.getBandWidth();
                oper_tree.clearBandWidth();

                for (int i = 0; i < oper_tree.getDepth(); i++) {
                    REQUIRE( (bw_1.getMaxWidth(i) <= bw_2.getMaxWidth(i)) );
                    REQUIRE( (bw_2.getMaxWidth(i) <= bw_3.getMaxWidth(i)) );
                }
            }
        }
    }
}

TEST_CASE("Apply identity convolution operator", "[apply_identity], [identity_convolution], [mw_operator]") {
    SECTION("1D") {
        applyIdentity<1>();
    }
    SECTION("2D") {
        applyIdentity<2>();
    }
    SECTION("3D") {
        applyIdentity<3>();
    }
}

template<int D> void applyIdentity() {
    double proj_prec = 1.0e-3;
    double apply_prec = 1.0e-3;
    double build_prec = 1.0e-4;

    MultiResolutionAnalysis<D> *mra = 0;
    GaussFunc<D> *fFunc = 0;
    initialize(&fFunc);
    initialize(&mra);

    MWProjector<D> Q(proj_prec);
    IdentityConvolution<D> I(*mra, build_prec);
    MWConvolution<D> apply(apply_prec);

    FunctionTree<D> fTree(*mra);
    FunctionTree<D> gTree(*mra);

    Q(fTree, *fFunc);
    apply(gTree, I, fTree);

    REQUIRE( (gTree.getDepth() <= fTree.getDepth()) );
    REQUIRE( (gTree.getNNodes() <= fTree.getNNodes()) );
    REQUIRE( (gTree.integrate() == Approx(fTree.integrate()).epsilon(apply_prec)) );

    finalize(&fFunc);
    finalize(&mra);
}

} // namespace
