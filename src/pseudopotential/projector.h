#pragma once
#include <mrchem.h>
#include "pseudopotential/sphericalHarmonics.h"

class ProjectorFunction{

public:
    ProjectorFunction(mrcpp::Coord<3> pos, double rl, int i, int l, int m, double prec);

    mrcpp::Coord<3> pos;
    int i;
    int l;
    int m;
    double rl;
    double prec;
    // mrcpp::ComplexFunction projector;
    std::shared_ptr<mrcpp::CompFunction<3>> projector_ptr;



    // destructor
    // ~ProjectorFunction(){
    //     std::cout << "Whya asdf "<< std::endl;
    //     projector.free(mrcpp::NUMBER::Total);
    // }

private:
    /**
     * @brief Contains analytic form of projector.
    */
    double (*s)(const std::array<double, 3> &r, const double &normr);

    void switch_sperics(int l, int m);

};
