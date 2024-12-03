#include "HirshfeldInterpolator.h"
#include "qmfunctions/density_utils.h"

double HirshfeldRadInterpolater::getNorm() const {
    double norm = 0.0;
    double xmax = 20;
    int n = 1000;
    double dx = xmax / n;
    for (int i = 0; i < n; i++) {
        double x = i * dx;
        double y = std::exp(evalf(x));
        norm += y * x * x;
    }
    norm *= dx * 4 * M_PI;
    return norm;
}

// Constructor
HirshfeldRadInterpolater::HirshfeldRadInterpolater(const std::string element, std::string data_dir, bool writeToFile) {
    Eigen::VectorXd rGrid;
    Eigen::VectorXd rhoGrid;

    std::string filename = data_dir + '/' + element + ".density";

    mrchem::density::readAtomicDensity(filename, rGrid, rhoGrid);

    rhoGrid = rhoGrid.array().log();

    lnRho = std::make_shared<interpolation_utils::PolyInterpolator>(rGrid, rhoGrid);
    if (writeToFile) {
        writeInterpolatedDensity(element + ".interpolated");
    }
}

// Function to evaluate the interpolated function
double HirshfeldRadInterpolater::evalf(const double &r) const {
    double y;
    return lnRho->evalfLeftNoRightLinear(r);
}

void HirshfeldRadInterpolater::writeInterpolatedDensity(const std::string path) {
    std::ofstream file;
    file.open(path);
    for (double r = 0.0; r < 50.0; r += 0.01) {
        double rho = this->evalf(r);
        file << r << " " << rho << std::endl;
    }
    file.close();
}