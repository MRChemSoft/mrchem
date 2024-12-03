/**
 * @file lebedev.cpp
 * @brief Implementation of the LebedevIntegrator class for surface integration calculations.
 */

#include "lebedev.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "surface_forces/LebedevData.h"


    /**
     * @brief Constructor for the LebedevIntegrator class.
     * @param filename The name of the file containing the Lebedev grid data.
     * @param radius The radius of the sphere.
     * @param center Center of the sphere.
     */
    LebedevIntegrator::LebedevIntegrator(int nPoints, double radius, const Eigen::Vector3d& center) {
        getLebedevData(nPoints);
        calculateCartesianPoints(radius, center);
    }

    /**
     * @brief Get the Cartesian coordinates of the Lebedev grid points.
     * @return A matrix containing the Cartesian coordinates of the Lebedev grid points.
     */
    Eigen::MatrixXd LebedevIntegrator::getPoints() const { return points; }

    /**
     * @brief Get the weights associated with each Lebedev grid point.
     * @return A vector containing the weights associated with each Lebedev grid point.
     */
    Eigen::VectorXd LebedevIntegrator::getWeights() const { return weights; }

    /**
     * @brief Get the normal vectors at each Lebedev grid point.
     * @return A matrix containing the normal vectors at each Lebedev grid point.
     */
    Eigen::MatrixXd LebedevIntegrator::getNormals() const { return normals; }

    /**
     * @brief Read the Lebedev grid data from a file.
     * @param filename The name of the file containing the Lebedev grid data.
     * @throw std::runtime_error if the file cannot be opened.
     */
    void LebedevIntegrator::readLebedevFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open the lebedev data file.");
        }
        std::string line;
        std::vector<double> thetaVec, phiVec, weightVec;

        while (getline(file, line)) {
            double theta, phi, weight;
            std::istringstream iss(line);
            if (!(iss >> theta >> phi >> weight)) { break; } // error

            thetaVec.push_back(theta * M_PI / 180.0);
            phiVec.push_back(phi * M_PI / 180.0);
            weightVec.push_back(weight);
        }

        n = thetaVec.size();
        points = Eigen::MatrixXd(n, 3);
        weights = Eigen::VectorXd(n);
        normals = Eigen::MatrixXd(n, 3);

        for (int i = 0; i < n; ++i) {
            weights[i] = weightVec[i];
            normals.row(i) << cos(thetaVec[i]) * sin(phiVec[i]), sin(thetaVec[i]) * sin(phiVec[i]), cos(phiVec[i]);
        }
    }

    /**
     * @brief Get the Lebedev data for a given number of grid points.
     * @param n The number of Lebedev grid points.
     */
    void LebedevIntegrator::getLebedevData(int nPoints) {
        std::tuple<Eigen::VectorXd, Eigen::Matrix3Xd> dat = lebedev(nPoints);
        weights = std::get<0>(dat);
        n = weights.size();
        normals = std::get<1>(dat).transpose();
        points = Eigen::MatrixXd(n, 3);
    }
        

    /**
     * @brief Calculate the Cartesian coordinates of the Lebedev grid points.
     * @param radius The radius of the sphere.
     * @param center Center of the sphere.
     */
    void LebedevIntegrator::calculateCartesianPoints(double radius, const Eigen::Vector3d &center) {
        for (int i = 0; i < n; ++i) {
            points.row(i) = radius * normals.row(i) + center.transpose();
        }
        weights *= 4.0 * M_PI * radius * radius;
    }


/*
This can be used for testing the LebedevIntegrator class.
Eigen::Vector3d pointEField(const Eigen::Vector3d& r) {
    double rNorm = r.norm();
    return r / (rNorm * rNorm * rNorm);
}
int main() {
    try {
        double radius = 4.0;
        Eigen::Vector3d shift(5.0, 2.0, 1.0);
        LebedevIntegrator integrator("lebedev.txt", radius, shift);

        double sum = integrator.getWeights().sum();
        std::cout << "Sum of weights: " << sum << " Refererence 4 * pi * r**2 " << 4 * M_PI * radius * radius << std::endl;

        Eigen::MatrixXd fields = integrator.getPoints();
        for (int i = 0; i < integrator.n; i++)
        {
            fields.row(i) = fields.row(i) / (fields.row(i).norm() * fields.row(i).norm() * fields.row(i).norm());
        }
        
        double integral = 0.0;
        for (int i = 0; i < integrator.n; i++)
        {
            integral += integrator.getWeights()[i] * fields.row(i).dot(integrator.getNormals().row(i));
        }

        std::cout << "Surface Integral: " << integral << " Reference 4 * pi " << 4 * M_PI << std::endl;

        Eigen::VectorXd pots = integrator.getPoints().rowwise().norm().cwiseInverse();
        pots = pots.cwiseProduct(pots);
        double potIntegral = integrator.getWeights().dot(pots);
        std::cout << "Potential Integral: " << potIntegral << " Reference 1 " << 1 << std::endl;
        

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    return 0;
}
*/