/**
 * @file sphericalHarmonics.h
 * @brief Contains the relevant, real spherical harmonics functions for the non local pseudopotential.
 * The naming of the functions is as follows: slm where l is the angular momentum and m is the magnetic quantum number. m stands for minus.
 * An example of the naming is s2m1 which is the spherical harmonic function with l=2 and m=-1.
 */

#include "pseudopotential/sphericalHarmonics.h"
#include <cmath>
#include <array>
#include <iostream>

// double (*)(const std::array<double, 3> &r, const double &normr) get_spherical_harmonics(const int &l, const int &m) {
double (*get_spherical_harmonics(const int &l, const int &m))(const std::array<double, 3> &r, const double &normr){
    if (l == 0 && m == 0) return s0;
    else if (l == 1 && m == -1) return s1m1;
    else if (l == 1 && m == 0) return s10;
    else if (l == 1 && m == 1) return s11;
    else if (l == 2 && m == -2) return s2m2;
    else if (l == 2 && m == -1) return s2m1;
    else if (l == 2 && m == 0) return s20;
    else if (l == 2 && m == 1) return s21;
    else if (l == 2 && m == 2) return s22;
    else if (l == 3 && m == -3) return s3m3;
    else if (l == 3 && m == -2) return s3m2;
    else if (l == 3 && m == -1) return s3m1;
    else if (l == 3 && m == 0) return s30;
    else if (l == 3 && m == 1) return s31;
    else if (l == 3 && m == 2) return s32;
    else if (l == 3 && m == 3) return s33;
    else if (l == 4 && m == -4) return s4m4;
    else if (l == 4 && m == -3) return s4m3;
    else if (l == 4 && m == -2) return s4m2;
    else if (l == 4 && m == -1) return s4m1;
    else if (l == 4 && m == 0) return s40;
    else if (l == 4 && m == 1) return s41;
    else if (l == 4 && m == 2) return s42;
    else if (l == 4 && m == 3) return s43;
    else if (l == 4 && m == 4) return s44;
    else {
        std::cerr << "Spherical harmonic function not found for l = " << l << " and m = " << m << std::endl;
        std::cerr << "Returning nullptr" << std::endl;
        return nullptr;
    }
}

double s0(const std::array<double, 3> &r, const double &normr) {
    return 0.5 * std::sqrt(1.0 / M_PI);
}

double s1m1(const std::array<double, 3> &r, const double &normr) {
    return std::sqrt(3.0 / (4.0 * M_PI)) * r[1];
}

double s10(const std::array<double, 3> &r, const double &normr) {
    return std::sqrt(3.0 / (4.0 * M_PI)) * r[2];
}

double s11(const std::array<double, 3> &r, const double &normr) {
    return std::sqrt(3.0 / (4.0 * M_PI)) * r[0];
}

double s2m2(const std::array<double, 3> &r, const double &normr) {
    return 0.5 * std::sqrt(15.0 / ( M_PI)) * r[0] * r[1];
}

double s2m1(const std::array<double, 3> &r, const double &normr) {
    return 0.5 * std::sqrt(15.0 / ( M_PI)) * r[1] * r[2];
}

double s20(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(5.0 / ( M_PI)) * (3.0 * r[2] * r[2] - normr * normr);
}

double s21(const std::array<double, 3> &r, const double &normr) {
    return 0.5 * std::sqrt(15.0 / ( M_PI)) * r[0] * r[2];
}

double s22(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(15.0 / ( M_PI)) * (r[0] * r[0] - r[1] * r[1]);
}

double s3m3(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(35.0 / ( 2 * M_PI)) * (r[1] * (3 * r[0] * r[0] - r[1] * r[1]));
}

double s3m2(const std::array<double, 3> &r, const double &normr) {
    return 0.5 * std::sqrt(105.0 / ( M_PI)) * r[0] * r[1] * r[2];
}

double s3m1(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(21.0 / ( 2 * M_PI)) * r[1] * (5 * r[2] * r[2] - normr * normr);
}

double s30(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(7.0 / ( M_PI)) * (5 * r[2] * r[2] * r[2] - 3 * r[2] * normr * normr );
}

double s31(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(21.0 / ( 2 * M_PI)) * r[0] * (5 * r[2] * r[2] - normr * normr);
}

double s32(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(105.0 / ( M_PI)) * (r[0] * r[0] - r[1] * r[1]) * r[2];
}

double s33(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(35.0 / ( 2 * M_PI)) * (r[0] * r[0] * r[0] - 3 * r[1] * r[1] * r[0]);
}

double s4m4(const std::array<double, 3> &r, const double &normr) {
    return .75 * std::sqrt(35.0 / ( M_PI)) * (r[0] * r[1]) * (r[0] * r[0] - r[1] * r[1]);
}

double s4m3(const std::array<double, 3> &r, const double &normr) {
    return .75 * std::sqrt(35.0 / ( 2 * M_PI)) * r[1] * (3 * r[0] * r[0] - r[1] * r[1]) * r[2];
}

double s4m2(const std::array<double, 3> &r, const double &normr) {
    return .75 * std::sqrt(5.0 / ( M_PI)) * r[0] * r[1] * (7 * r[2] * r[2] - normr * normr);
}

double s4m1(const std::array<double, 3> &r, const double &normr) {
    return .75 * std::sqrt(5.0 / ( 2 * M_PI)) * r[1] * (7 * r[2] * r[2] * r[2] - 3 * r[2] * normr * normr);
}

double s40(const std::array<double, 3> &r, const double &normr) {
    return (3. / 16.) * std::sqrt(1.0 / ( M_PI)) * (35 * r[2] * r[2] * r[2] * r[2] - 30 * r[2] * r[2] * normr * normr + 3 * normr * normr * normr * normr);
}

double s41(const std::array<double, 3> &r, const double &normr) {
    return .75 * std::sqrt(5.0 / (2.0 * M_PI)) * r[0] * (7 * r[2] * r[2] * r[2] - 3 * r[2] * normr * normr);
}

double s42(const std::array<double, 3> &r, const double &normr) {
    return (3. / 8.) * std::sqrt(5.0 / ( M_PI)) * (r[0] * r[0] - r[1] * r[1]) * (7 * r[2] * r[2] - normr * normr);
}

double s43(const std::array<double, 3> &r, const double &normr) {
    return .75 * std::sqrt(35.0 / ( 2 * M_PI)) * r[0] * (r[0] * r[0] - 3 * r[1] * r[1]) * r[2];
}

double s44(const std::array<double, 3> &r, const double &normr) {
    return (3. / 16.) * std::sqrt(35.0 / ( M_PI)) * ( 
        r[0] * r[0] * ( r[0] * r[0] - 3 * r[1] * r[1] ) -
        r[1] * r[1] * ( 3 * r[0] * r[0] - r[1] * r[1] ));
}