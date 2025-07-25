#include "Factory.h"


#include "LibXC.h"
#include "LDA.h"



namespace mrlibxc {

Factory::Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA)
        : mra(MRA) {}

Factory::~Factory() {
    cleanupFunctionals();
}

void Factory::cleanupFunctionals() {
    for (auto &func_data : functionals) {
        if (func_data.initialized) {
            xc_func_end(&func_data.func);
            func_data.initialized = false;
        }
    }
    functionals.clear();
}

int Factory::mapFunctionalName(const std::string &name) const {
    // Map common functional names to LibXC IDs
    if (name == "LDA" || name == "LDA_X" || name == "svwn5") return XC_LDA_X;

    
    // If not a common name, try to get it directly from LibXC
    int func_id = xc_functional_get_number(name.c_str());
    if (func_id <= 0) {
        std::string msg = "Unknown functional: " + name;
        MSG_ABORT(msg.c_str());
    }
    return func_id;
}

void Factory::setFunctional(const std::string &name, double weight) {
    int func_id = mapFunctionalName(name);
    
    LibXCData func_data;
    func_data.func_id = func_id;
    func_data.weight = weight;
    func_data.initialized = false;
    
    functionals.push_back(func_data);
}



/** @brief Build a MRDFT object from the currently defined parameters */
std::unique_ptr<mrdft::MRDFT> Factory::build() {
    // Init DFT grid
    auto grid_p = std::make_unique<mrdft::Grid>(mra);
    
    // Initialize LibXC functionals
    for (auto &func_data : functionals) {
        int polarization = spin ? XC_POLARIZED : XC_UNPOLARIZED;
        if (xc_func_init(&func_data.func, func_data.func_id, polarization) != 0) {
            std::string msg = "Error initializing LibXC functional";
            MSG_ABORT(msg.c_str());
        }
        func_data.initialized = true;
    }
    
    // Check if we have any functionals
    if (functionals.empty()) {
        MSG_ABORT("No functionals defined");
    }
    
    
    // Init XC functional
    std::unique_ptr<mrdft::Functional> func_p{nullptr};

    func_p = std::make_unique<mrdft::LDA>(order, functionals);
    
    if (func_p == nullptr) MSG_ABORT("Invalid functional type");
    
    diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
    func_p->setDerivOp(diff_p);
    func_p->setLogGradient(log_grad);
    func_p->setDensityCutoff(cutoff);
    
    auto mrdft_p = std::make_unique<mrdft::MRDFT>(grid_p, func_p);
    return mrdft_p;
}

} // namespace mrdft