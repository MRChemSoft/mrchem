#pragma once
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <filesystem>
#include <iostream>

namespace interpolation_utils{

/**
 * @brief Interpolate a 5th order polynomial through 5 points and evaluate at x.
 * @param x_in x values of the 5 points
 * @param y_in y values of the 5 points
 * @param x x value at which to evaluate the polynomial
 */
double polynomialInterpolate5(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in, double x);

/**
 * @brief Evaluate the derivative of a 5th order polynomial through 5 points at x.
 * @param x_in x values of the 5 points
 * @param y_in y values of the 5 points
 * @param x x value at which to evaluate the derivative of the polynomial
 */
double polynomialInterpolate5_deriv(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in, double x);

/**
 * @brief Binary search for the index of the largest element in x that is smaller than x0.
 * @param x Vector to search in
 * @param x0 Value to search for
 */
int binarySearch(const Eigen::VectorXd &x, const double &x0);

/**
 * @brief Class to interpolate a functions using a 5th order polynomials.
 */
class PolyInterpolator {
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    int n;
    double xmin;
    double xmax;
    /**
     * @brief Derivative of the interpolated function at xmax. Used for linear extrapolation.
     */
    double ypxmax;

    /**
     * @brief Derivative of the interpolated function at xmin. Used for linear extrapolation.
     */
    double ypxmin;

    public:
    /**
     * @brief Constructor
     * @param x_in x values of the points to interpolate
     * @param y_in y values of the points to interpolate
     */
    PolyInterpolator(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in){
        x = x_in;
        y = y_in;
        n = x_in.size();
        xmin = x_in(0);
        xmax = x_in(x_in.size() - 1);
        Eigen::VectorXd x_in_poly(5), y_in_poly(5);

        int j = 0;
        for (int i = n - 5; i < n; i++) {
            x_in_poly(j) = x(i);
            y_in_poly(j) = y(i);
            j++;
        }
        ypxmax = polynomialInterpolate5_deriv(x_in_poly, y_in_poly, xmax);
        for (int i = 0; i < 5; i++) {
            x_in_poly(i) = x(i);
            y_in_poly(i) = y(i);
        }
        ypxmin = polynomialInterpolate5_deriv(x_in_poly, y_in_poly, xmin);
    }
    
    /**
     * @brief Evaluate the interpolated function at x. 
     * No extrapolation for x < xmin (meaning that the polynomial is evaluated at x < xmin), linear extrapolation for x > xmax.
     * Useful when the logarithm of the density is interpolated.
     * @param x x value at which to evaluate the function
     */
    double evalfLeftNoRightLinear(const double &xval) const{
        double y;
        if (xval > xmax) { // linear extrapolation
                y = this->y(n - 1) + ypxmax*(xval - xmax);
            return y;
        }

        Eigen::VectorXd x_in_poly(5), y_in_poly(5);
        int i = adjustIndexToBoundaries(binarySearch(this->x, xval));
        for (int k = 0; k < 5; k++) {
            x_in_poly(k) = this->x(i - 2 + k);
            y_in_poly(k) = this->y(i - 2 + k);
        }
        y = polynomialInterpolate5(x_in_poly, y_in_poly, xval);
        return y;
    }

    /**
     * @brief Evaluate the interpolated function at x.
     * No extrapolation for x < xmin (meaning that the polynomial is evaluated at x < xmin), zero is returned for x > xmax.
     * Useful when the density is interpolated.
     * @param x x value at which to evaluate the function
     */
    double evalfLeftNoRightZero(const double &xval) const{
        double y;
        if (xval > xmax) { // constant
            y = 0.0;
            return y;
        }

        Eigen::VectorXd x_in_poly(5), y_in_poly(5);
        int j = adjustIndexToBoundaries(binarySearch(this->x, xval));
        int i = adjustIndexToBoundaries(binarySearch(this->x, xval));
        for (int k = 0; k < 5; k++) {
            x_in_poly(k) = this->x(i - 2 + k);
            y_in_poly(k) = this->y(i - 2 + k);
        }
        y = polynomialInterpolate5(x_in_poly, y_in_poly, xval);
        return y;
    }

    /**
     * @brief Evaluate the interpolated function at x.
     * No extrapolation for x < xmin (meaning that the polynomial is evaluated at x < xmin), the last value is returned for x > xmax.
     */
    double evalfLeftNoRightConstant(const double &xval) const {
        double y;
        if (xval > xmax) { // constant
            y = this->y(n - 1);
            return y;
        }

        Eigen::VectorXd x_in_poly(5), y_in_poly(5);
        int i = adjustIndexToBoundaries(binarySearch(this->x, xval));
        for (int k = 0; k < 5; k++) {
            x_in_poly(k) = this->x(i - 2 + k);
            y_in_poly(k) = this->y(i - 2 + k);
        }
        y = polynomialInterpolate5(x_in_poly, y_in_poly, xval);
        return y;
    }

    private: 
    int adjustIndexToBoundaries(int i) const {
        int j = i;
        if (i == 0) j = 2;
        if (i == 1) j = 2;
        if (i == n - 1) j = n - 3;
        if (i == n - 2) j = n - 3;
        return j;
    }
};    
} // namespace interpolation_utils