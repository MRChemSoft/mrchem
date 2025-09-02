#pragma once
/**
 * @file sphericalHarmonics.h
 * @brief Contains the relevant, real spherical harmonics functions for the non local pseudopotential.
 * The naming of the functions is as follows: slm where l is the angular momentum and m is the magnetic quantum number. m stands for minus.
 * An example of the naming is s2m1 which is the spherical harmonic function with l=2 and m=-1.
 * Important: functions are not normalized, the factor r^-l is not included for numerical stability.
 */
#include <array>

// double (*)(const std::array<double, 3> &r, const double &normr) get_spherical_harmonics(const int &l, const int &m);
double (*get_spherical_harmonics(const int &l, const int &m))(const std::array<double, 3> &r, const double &normr);


double s0(const std::array<double, 3> &r, const double &normr);

double s1m1(const std::array<double, 3> &r, const double &normr);

double s10(const std::array<double, 3> &r, const double &normr);

double s11(const std::array<double, 3> &r, const double &normr);

double s2m2(const std::array<double, 3> &r, const double &normr);

double s2m1(const std::array<double, 3> &r, const double &normr);

double s20(const std::array<double, 3> &r, const double &normr);

double s21(const std::array<double, 3> &r, const double &normr);

double s22(const std::array<double, 3> &r, const double &normr);

double s3m3(const std::array<double, 3> &r, const double &normr);

double s3m2(const std::array<double, 3> &r, const double &normr);

double s3m1(const std::array<double, 3> &r, const double &normr);

double s30(const std::array<double, 3> &r, const double &normr);

double s31(const std::array<double, 3> &r, const double &normr);

double s32(const std::array<double, 3> &r, const double &normr);

double s33(const std::array<double, 3> &r, const double &normr);

double s4m4(const std::array<double, 3> &r, const double &normr);

double s4m3(const std::array<double, 3> &r, const double &normr);

double s4m2(const std::array<double, 3> &r, const double &normr);

double s4m1(const std::array<double, 3> &r, const double &normr);

double s40(const std::array<double, 3> &r, const double &normr);

double s41(const std::array<double, 3> &r, const double &normr);

double s42(const std::array<double, 3> &r, const double &normr);

double s43(const std::array<double, 3> &r, const double &normr);

double s44(const std::array<double, 3> &r, const double &normr);