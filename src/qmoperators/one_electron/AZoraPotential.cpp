#include "AZoraPotential.h"

#include "mrchem.h"
#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"
#include <vector>
#include <fstream>
#include <Eigen/Core>
#include <string>
#include "utils/PolyInterpolator.h"
#include "qmoperators/QMOperator.h"
#include "qmoperators/QMPotential.h"

/**
 * @brief Read ZORA potential from file. Check if file exists and abort if it does not.
 * @param path Path to the file containing the ZORA potential
 * @param rGrid Vector containing the radial grid
 * @param vZora Vector containing the ZORA potential
 * @param kappa Vector containing the kappa parameter
*/
void readZoraPotential(const std::string path, Eigen::VectorXd &rGrid, Eigen::VectorXd &vZora, Eigen::VectorXd &kappa){
    std::vector<double> r, v, k;
    bool file_exists = std::filesystem::exists(path);
    if (!file_exists) {
        std::cerr << "File " << path << " does not exist." << std::endl;
        std::cout << "File " << path << " does not exist." << std::endl;
        exit(1);
    }
    std::ifstream file(path);
    std::string line;
    double r_, v_, k_;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        iss >> k_ >> r_ >> v_;
        r.push_back(r_);
        v.push_back(v_);
        k.push_back(k_);
    }
    file.close();
    rGrid = Eigen::Map<Eigen::VectorXd>(r.data(), r.size());
    vZora = Eigen::Map<Eigen::VectorXd>(v.data(), v.size());
    kappa = Eigen::Map<Eigen::VectorXd>(k.data(), k.size());
    // The kappa function is half of what is defined in the paper Scalar 
    // Relativistic Effects with Multiwavelets: Implementation and Benchmark
    // it is not used in the code, only the potential is used
}

namespace mrchem {

    /**
     * Initialize the azora potential based on the molecule.
     * This method would typically setup the real and imaginary function trees
     * representing the potential.
     */
    void AZoraPotential::initAzoraPotential() {

        int n = nucs.size();
        Eigen::VectorXd rGrid, vZora, kappa;
        atomicPotentials.clear();

        // std::vector<RadInterpolater> atomicPotentials;

        std::string element;
        std::string filename;
        for (int i = 0; i < n; i++) {
            element = nucs[i].getElement().getSymbol();
            filename = this->azora_dir + '/' + element + ".txt";
            readZoraPotential(filename, rGrid, vZora, kappa);
            interpolation_utils::PolyInterpolator potential_interpolator(rGrid, vZora);

            atomicPotentials.push_back(potential_interpolator);
        }
        // project(this->prec);
    }

    void AZoraPotential::project(double proj_prec){
        if (isProjected) free(mrchem::NUMBER::Total);
        mrcpp::ComplexFunction vtot;
        this->prec = proj_prec;
        auto chi_analytic = [this](const mrcpp::Coord<3>& r) {
            return this->evalf_analytic(r);
        };
        mrcpp::cplxfunc::project(vtot, chi_analytic, mrcpp::NUMBER::Real, proj_prec);
        this->add(1.0, vtot);
        isProjected = true;
    }
    
    double AZoraPotential::evalf_analytic(const mrcpp::Coord<3>& r){
        double V = 0.0;
        // Loop over all atoms:
        for (int i = 0; i < this->atomicPotentials.size(); i++) {
            mrcpp::Coord<3> r_i = nucs[i].getCoord();
            double rr = std::sqrt((r[0] - r_i[0]) * (r[0] - r_i[0]) + (r[1] - r_i[1]) * (r[1] - r_i[1]) + (r[2] - r_i[2]) * (r[2] - r_i[2]));
            V += this->atomicPotentials[i].evalfLeftNoRightConstant(rr);
        }
        return 1 / (1 - V / (2.0 * c * c)) - 1;
    }

} // namespace mrchem
