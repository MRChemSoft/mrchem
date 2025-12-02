#pragma once

#include <nlohmann/json.hpp>

#include "mrchem.h"

#include "utils/print_utils.h"
#include <string>

namespace mrchem {

class PopulationAnalysis final {
public:
    DoubleMatrix getMatrix() const { return this->populations; }

    void setMatrix(const DoubleMatrix &v) { this->populations = v; }

    void print(const std::string &id) const {
        mrcpp::print::header(0, "Orbital Populations (" + id + ")");
        mrcpp::print::separator(0, '-');
        if (populations.cols() == 1) {
            for (int i = 0; i < populations.rows(); i++) {
                std::string text = "Orbital " + std::to_string(i);
                print_utils::scalar(0, text, populations(i, 0));
            }
        }
        else {
            for (int i = 0; i < populations.rows(); i++) {
                std::string text = "Orbital " + std::to_string(i);
                mrcpp::Coord<3> vec = {populations(i, 0), populations(i, 1), populations(i, 2)};
                print_utils::coord(0, text, vec);
            }
        }
        mrcpp::print::separator(0, '=');
    }

    nlohmann::json json() const { return {{"total", print_utils::eigen_to_vector(getMatrix(), 1.0e-12)}}; }

protected:
    DoubleMatrix populations;
};

} // namespace mrchem