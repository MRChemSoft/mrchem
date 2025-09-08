/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 *
 * This file is part of MRChem.
 */

#include "Factory.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <XCFun/xcfun.h>

#include "GGA.h"
#include "Grid.h"
#include "LDA.h"
#include "MRDFT.h"
#include "SpinGGA.h"
#include "SpinLDA.h"

// LibXC adapters (compiled only if enabled at configure time)
#ifdef MRCHEM_ENABLE_LIBXC
#  include "LibXCBackend.h"
#endif

namespace mrdft {

// ------------------------------- helpers ---------------------------------

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    return s;
}

static std::vector<int> parse_ids_csv(const char* csv_env) {
    std::vector<int> out;
    if (!csv_env) return out;
    std::stringstream ss(csv_env);
    for (std::string tok; std::getline(ss, tok, ','); ) {
        if (!tok.empty()) out.push_back(std::stoi(tok));
    }
    return out;
}

// ------------------------------ Factory ----------------------------------

Factory::Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA)
    : mra(MRA)
    , xcfun_p(xcfun_new(), xcfun_delete) {}

/** @brief Build a MRDFT object from the currently defined parameters
 *
 * Backend selection policy:
 *  - Default: XCFun
 *  - If env MRCHEM_XC_BACKEND=libxc (or MRCHEM_FORCE_LIBXC=1), choose LibXC (if compiled)
 *  - LibXC IDs are taken from Factory::libxc_ids if set; otherwise from env MRCHEM_LIBXC_IDS
 *  - If LibXC chosen but no IDs available, abort with a clear message
 */
std::unique_ptr<MRDFT> Factory::build() {
    // Init DFT grid
    auto grid_p = std::make_unique<Grid>(mra);

    // Decide backend from member + env (env wins)
    std::string backend_eff = to_lower(backend);
    if (const char* be = std::getenv("MRCHEM_XC_BACKEND")) {
        backend_eff = to_lower(std::string(be));
    }
    const bool force_libxc = (std::getenv("MRCHEM_FORCE_LIBXC") != nullptr);
    const bool want_libxc  = force_libxc || (backend_eff == "libxc");

    // Determine LibXC IDs if needed
    std::vector<int> ids = libxc_ids; // may be pre-filled by caller
    if (ids.empty()) {
        // env override
        auto env_ids = parse_ids_csv(std::getenv("MRCHEM_LIBXC_IDS"));
        if (!env_ids.empty()) ids = std::move(env_ids);
    }

    // Flag whether LibXC adapters are compiled in
#ifdef MRCHEM_ENABLE_LIBXC
    const bool libxc_available = true;
#else
    const bool libxc_available = false;
#endif
    const bool use_libxc = (want_libxc && libxc_available);

    // ---------------------------------------------------------------------
    // XCFun user-eval setup (parity with upstream)
    // ---------------------------------------------------------------------
    bool gga = xcfun_is_gga(xcfun_p.get());
    bool lda = !gga;

    unsigned int mode          = 1;                    // partial derivative mode
    unsigned int func_type     = (gga ? 1u : 0u);      // LDA(0) or GGA(1)
    unsigned int dens_type     = 1 + spin;             // n (1) or alpha&beta (2)
    unsigned int laplacian     = 0;                    // no laplacian
    unsigned int kinetic       = 0;                    // no kinetic energy density
    unsigned int current       = 0;                    // no current density
    unsigned int exp_deriv     = !gamma;               // use gamma or explicit derivatives
    if (!gga) exp_deriv = 0;                           // fall back to gamma-type if LDA

    xcfun_user_eval_setup(xcfun_p.get(), order, func_type, dens_type,
                          mode, laplacian, kinetic, current, exp_deriv);

    // ---------------------------------------------------------------------
    // MW derivative operator if GGA
    // ---------------------------------------------------------------------
    if (gga) {
        if      (diff_s == "bspline") diff_p = std::make_unique<mrcpp::BSOperator<3>>(mra, 1);
        else if (diff_s == "abgv_00") diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
        else if (diff_s == "abgv_55") diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.5, 0.5);
        // else leave null; Functional classes can validate as needed
    }

    // ---------------------------------------------------------------------
    // Build Functional (choose LibXC adapters when requested & available)
    // ---------------------------------------------------------------------
    std::unique_ptr<Functional> func_p{nullptr};

    if (spin) {
        if (gga) {
#ifdef MRCHEM_ENABLE_LIBXC
            if (use_libxc) {
                if (ids.empty()) {
                    MSG_ABORT("LibXC backend selected but no LibXC IDs specified. "
                              "Set MRCHEM_LIBXC_IDS (e.g. 402 for B3LYP, or '106,131' for BLYP), "
                              "or call Factory::setLibXCIDs().");
                }
                func_p = std::make_unique<LibXCSpinGGA>(order, xcfun_p, diff_p, ids);
            } else
#endif
            {
                func_p = std::make_unique<SpinGGA>(order, xcfun_p, diff_p);
            }
        } else if (lda) {
#ifdef MRCHEM_ENABLE_LIBXC
            if (use_libxc) {
                if (ids.empty()) {
                    MSG_ABORT("LibXC backend selected but no LibXC IDs specified. "
                              "Set MRCHEM_LIBXC_IDS (e.g. '1,9' for LDA X+PW92 C), "
                              "or call Factory::setLibXCIDs().");
                }
                func_p = std::make_unique<LibXCSpinLDA>(order, xcfun_p, ids);
            } else
#endif
            {
                func_p = std::make_unique<SpinLDA>(order, xcfun_p);
            }
        }
    } else {
        if (gga) {
#ifdef MRCHEM_ENABLE_LIBXC
            if (use_libxc) {
                if (ids.empty()) {
                    MSG_ABORT("LibXC backend selected but no LibXC IDs specified. "
                              "Set MRCHEM_LIBXC_IDS (e.g. 402 for B3LYP, or '106,131' for BLYP), "
                              "or call Factory::setLibXCIDs().");
                }
                func_p = std::make_unique<LibXCGGA>(order, xcfun_p, diff_p, ids);
            } else
#endif
            {
                func_p = std::make_unique<GGA>(order, xcfun_p, diff_p);
            }
        } else if (lda) {
#ifdef MRCHEM_ENABLE_LIBXC
            if (use_libxc) {
                if (ids.empty()) {
                    MSG_ABORT("LibXC backend selected but no LibXC IDs specified. "
                              "Set MRCHEM_LIBXC_IDS (e.g. '1,9' for LDA X+PW92 C), "
                              "or call Factory::setLibXCIDs().");
                }
                func_p = std::make_unique<LibXCLDA>(order, xcfun_p, ids);
            } else
#endif
            {
                func_p = std::make_unique<LDA>(order, xcfun_p);
            }
        }
    }

    if (!func_p) {
        MSG_ABORT("Invalid or unsupported functional type (spin/LDA/GGA combination).");
    }

    // ---------------------------------------------------------------------
    // Keep original behavior: set derivative operator & knobs on the functional
    // ---------------------------------------------------------------------
    // (Upstream MRChem sets ABGV(0,0) as a default fallback.)
    diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
    func_p->setDerivOp(diff_p);
    func_p->setLogGradient(log_grad);
    func_p->setDensityCutoff(cutoff);

    auto mrdft_p = std::make_unique<MRDFT>(grid_p, func_p);
    return mrdft_p;
}

} // namespace mrdft
