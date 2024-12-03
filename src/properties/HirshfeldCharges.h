#pragma once

#include <nlohmann/json.hpp>

#include "mrchem.h"

#include "utils/math_utils.h"
#include "utils/print_utils.h"
#include <string>

namespace mrchem {

class HirshfeldCharges final {
public:

    DoubleVector getVector() const { return this->hirshfeld_charges; }

    void setVector(const DoubleVector &v) { this->hirshfeld_charges = v; }

    void print(const std::string &id) const {
        mrcpp::print::header(0, "Hirshfeld Charges (" + id + ")");
        mrcpp::print::separator(0, '-');    
        for (int i = 0; i < hirshfeld_charges.size(); i++) {
            std::string text = "Charge of atom " + std::to_string(i);
            print_utils::scalar(0, text, hirshfeld_charges(i));
        }
        mrcpp::print::separator(0, '-');
        print_utils::scalar(0, "Sum of Hirshfeld charges", getVector().sum(), "(au)", -1, true);
        mrcpp::print::separator(0, '=');
    }

    nlohmann::json json() const {
        return {{"total", print_utils::eigen_to_vector(getVector(), 1.0e-12)}};
    }

protected:
    DoubleVector hirshfeld_charges;
};

} // namespace mrchem