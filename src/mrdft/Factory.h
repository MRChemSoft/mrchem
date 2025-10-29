#pragma once

#include <MRCPP/MWOperators>
#include <XCFun/xcfun.h>
#include "MRDFT.h"

namespace mrdft {

using XC_p = std::unique_ptr<xcfun_t, decltype(&xcfun_delete)>;

class Factory final {
public:
    Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA);
    ~Factory() = default;

    void setSpin(bool s) { spin = s; }
    void setOrder(int k) { order = k; }
    void setUseGamma(bool g) { gamma = g; }
    void setLogGradient(bool lg) { log_grad = lg; }
    void setDensityCutoff(double c) { cutoff = c; }
    void setDerivative(const std::string &n) { diff_s = n; }
    void setFunctional(const std::string &n, double c = 1.0) { xcfun_set(xcfun_p.get(), n.c_str(), c); }

    // Optional backend controls (env can override):
    //   MRCHEM_XC_BACKEND=libxc|xcfun
    //   MRCHEM_LIBXC_IDS="ID[,ID,...]"  (or family-specific envs)
    void setBackend(const std::string &b) { backend = b; }
    void setLibXCIDs(const std::vector<int> &ids_in) { libxc_ids = ids_in; }

    std::unique_ptr<MRDFT> build();

private:
    int order{1};
    bool spin{false};
    bool gamma{false};
    bool log_grad{false};
    double cutoff{-1.0};
    std::string diff_s{"abgv_00"};
    std::string backend{"xcfun"};           // "xcfun" (default) or "libxc"
    std::vector<int> libxc_ids;             // used only if backend == "libxc"
    const mrcpp::MultiResolutionAnalysis<3> mra;

    XC_p xcfun_p;
    std::unique_ptr<mrcpp::DerivativeOperator<3>> diff_p;
};

} // namespace mrdft
