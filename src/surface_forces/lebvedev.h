#ifndef LEBVEDEV_H
#define LEBVEDEV_H

#include <Eigen/Dense>
#include <string>

class LebedevIntegrator {
public:
    Eigen::MatrixXd points; /**< Matrix storing the Cartesian coordinates of the Lebedev grid points. */
    Eigen::VectorXd weights; /**< Vector storing the weights associated with each Lebedev grid point. */
    Eigen::MatrixXd normals; /**< Matrix storing the normal vectors (to the integration sphere) at each Lebedev grid point. */
    int n; /**< Number of Lebedev grid points. */


    LebedevIntegrator(const std::string& filename, double radius, const Eigen::Vector3d& center);
    Eigen::MatrixXd getPoints() const;
    Eigen::VectorXd getWeights() const;
    Eigen::MatrixXd getNormals() const;
private:
    void readLebedevFile(const std::string& filename);
    void calculateCartesianPoints(double radius, const Eigen::Vector3d &center);
};

#endif // LEBVEDEV_H