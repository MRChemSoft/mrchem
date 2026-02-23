#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include "MRCPP/Printer"
#include "utils/print_utils.h"

/**
 * Splits a given string into a vector of words.
 *
 * @param str The string to be split.
 * @return A vector of words obtained from the input string.
 */
inline std::vector<std::string> splitStringToWords(const std::string& str) {
    std::istringstream iss(str);
    std::vector<std::string> words;
    std::string word;

    while (iss >> word) {
        words.push_back(word);
    }

    return words;
}


/**
 * @class PseudopotentialData
 * @brief Class representing pseudopotential data.
 * 
 * This class stores the data related to a pseudopotential, including the effective charge of the nucleus,
 * the atomic number of the nucleus, the radius of the local part of the pseudopotential, the coefficients
 * of the local functions, the radii of the projectors, the maximum angular momentum, the number of projectors,
 * the projector matrices, and the number of separable components.
 */
class PseudopotentialData {

public:




    PseudopotentialData(nlohmann::json pp_json_in) {

        nlohmann::json pp_json = pp_json_in["pseuodopotential"];

        zeff = pp_json["zeff"];
        zion = pp_json["zion"];
        rloc = pp_json["local"]["rloc"];
        nloc = pp_json["local"]["nloc"];
        alpha_pp = pp_json["local"]["alpha_pp"];
        std::vector<double> c_vec = pp_json["local"]["c"];
        c = Eigen::Map<Eigen::VectorXd>(c_vec.data(), c_vec.size());
        nsep = pp_json["nonlocal"]["nsep"];
        std::vector<double> rl_vec = pp_json["nonlocal"]["rl"];
        rl = rl_vec;
        std::vector<int> dim_h_vec = pp_json["nonlocal"]["dim_h"];
        dim_h = dim_h_vec;
        for (int l = 0; l < nsep; l++) {
            std::vector<std::vector<double>> h_l_vec = pp_json["nonlocal"]["h"][l];
            Eigen::MatrixXd h_l_mat(dim_h[l], dim_h[l]);
            for (int i = 0; i < dim_h[l]; i++) {
                for (int j = 0; j < dim_h[l]; j++) {
                    h_l_mat(i, j) = h_l_vec[i][j];
                }
            }
            h.push_back(h_l_mat);
        }

        if (pp_json.contains("nlcc")) {
            this->has_nlcc = true;
            this->qnlcc = pp_json["nlcc"]["qcore"];
            this->rnlcc = pp_json["nlcc"]["rcore"];
        }
    }

    /**
     * Prints the pseudopotential.
     */
    void print() {
        mrcpp::print::separator(0, '-');
        mrcpp::print::header(0, "Pseudopotential data");
        mrchem::print_utils::scalar(0, "Zeff", zeff);
        mrchem::print_utils::scalar(0, "Zion", zion);
        mrchem::print_utils::scalar(0, "rloc", rloc);
        mrchem::print_utils::scalar(0, "alpha_pp", alpha_pp);
        mrchem::print_utils::scalar(0, "nloc", nloc);
        mrchem::print_utils::vector(0, "c", c);
        mrchem::print_utils::scalar(0, "nsep", nsep);

        for (int l = 0; l < nsep; l++) {
            mrcpp::print::separator(0, '+');
            mrchem::print_utils::scalar(0, "l", l);
            mrchem::print_utils::scalar(0, "rl", rl[l]);
            mrchem::print_utils::matrix(0, "h", h[l]);
        }
        mrcpp::print::separator(0, '+');
        if (has_nlcc) {
            mrchem::print_utils::text(0, "NLCC", "enabled");
            mrchem::print_utils::scalar(0, "qnlcc", qnlcc);
            mrchem::print_utils::scalar(0, "rnlcc", rnlcc);
        }
        mrcpp::print::separator(0, '-');
    }

    /**
     * Returns the effective charge of the nucleus.
     *
     * @return The effective charge of the nucleus.
     */
    int getZeff() const {
        return zeff;
    }

    /**
     * Returns the atomic number of the nucleus.
     *
     * @return The atomic number of the nucleus.
     */
    int getZion() const {
        return zion;
    }

    /**
     * Returns the radius of the local part of the pseudopotential.
     *
     * @return The radius of the local part of the pseudopotential.
     */
    double getRloc() const {
        return rloc;
    }

    /**
     * Returns the number of local functions.
     *
     * @return The number of local functions.
     */
    int getNloc() const {
        return nloc;
    }

    /**
     * Returns the coefficients of the local functions.
     *
     * @return The coefficients of the local functions.
     */
    Eigen::VectorXd getC() const {
        return c;
    }

    /**
     * Returns the radii of the projectors.
     *
     * @return The radii of the projectors.
     */
    std::vector<double> getRl() const {
        return rl;
    }

    /**
     * Returns the projector matrices.
     *
     * @return The projector matrices.
     */
    std::vector<Eigen::MatrixXd> getH() const {
        return h;
    }

    /**
     * Returns the dimension of the projector matrices.
     *
     * @return The dimension of the projector matrices.
     */
    std::vector<int> getDimH() const {
        return dim_h;
    }

    /**
     * Returns the number of separable components.
     *
     * @return The number of separable components.
     */
    int getNsep() const {
        return nsep;
    }

    /**
     * Returns whether the pseudopotential has a non-local core correction.
     *
     * @return Whether the pseudopotential has a non-local core correction.
     */
    bool getHasNlcc() const {
        return has_nlcc;
    }

    /**
     * Returns the radius of the non-local core correction.
     */
    double getRnlcc() const {
        if (!has_nlcc) {
            MSG_WARN("Pseudopotential does not have a non-local core correction but getRnlcc() is called");
        }
        return rnlcc;
    }

    /**
     * Returns the charge of the non-local core correction. (c_core in the paper)
     */
    double getQnlcc() const {
        if (!has_nlcc) {
            MSG_WARN("Pseudopotential does not have a non-local core correction but getQnlcc() is called");
        }
        return qnlcc;
    }

    /**
     * Returns the maximum radius of the pseudopotential.
     */
    double getMaxRpp() const {
        double max_r_pp = rloc;
        for (int l = 0; l < nsep; l++) {
            if (rl[l] > max_r_pp) {
                max_r_pp = rl[l];
            }
        }
        if (getHasNlcc()) {
            if (rnlcc > max_r_pp) {
                max_r_pp = rnlcc;
            }
        }
        return max_r_pp;
    }

    int zeff; /** Effective charge of nucleus */
    int zion; /** Atomic number of nucleus */
    double rloc; /** Radius of local part of pseudopotential */
    double alpha_pp;
    int nloc; /** Number of local functions */
    Eigen::VectorXd c; /** Coefficienst of local functions */
    std::vector<double> rl; /** Radii of projectors (one for each angular momentum) */
    std::vector<Eigen::MatrixXd> h; /** Projector matrices */
    std::vector<int> dim_h; /** Dimension of projector matrices */
    int nsep; /** Number of different angular momenta from 0 to nsep - 1*/

    private:
    
    bool has_nlcc = false; /** Whether the pseudopotential has a non-local core correction */
    double rnlcc = 0.0; /** Radius of non-local core correction */
    double qnlcc = 0.0; /** Charge of non-local core correction */

};

