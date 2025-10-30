/*
* MRChem — XC Factory (LibXC/XCFun backend selection)
*/
#include "Factory.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <limits>
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

#ifdef MRCHEM_ENABLE_LIBXC
#  include "LibXCBackend.h"
#  include <xc.h>
#endif

namespace mrdft {

// ---------- tiny helpers ----------
static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    return s;
}
static const char* getenv_nn(const char* k) {
    const char* v = std::getenv(k);
    return (v && *v) ? v : nullptr;
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

#ifdef MRCHEM_ENABLE_LIBXC
// Helpers for LibXC ID family filtering
static int libxc_family_of_id(int id) {
    int fam = XC_FAMILY_UNKNOWN, num = -1;
    const int rc = xc_family_from_id(id, &fam, &num);
    return (rc == 0) ? fam : XC_FAMILY_UNKNOWN;
}
static bool fam_is_gga_like(int fam) {
    return (fam == XC_FAMILY_GGA || fam == XC_FAMILY_HYB_GGA);
}
// Fail early if user passed MGGA/HYB-MGGA IDs to the GGA backend
static void guard_gga_ids_family(const std::vector<int>& ids, int nspin) {
    for (int id : ids) {
        xc_func_type f;
        if (xc_func_init(&f, id, nspin) != 0) {
            throw std::runtime_error("LibXC: could not initialize functional id=" + std::to_string(id));
        }
        const int fam = f.info ? f.info->family : XC_FAMILY_UNKNOWN;
        xc_func_end(&f);
        if (!fam_is_gga_like(fam)) {
            throw std::runtime_error(
                "LibXC: selected GGA backend but LibXC id " + std::to_string(id) +
                " is not GGA/HYB-GGA. Use an LDA or MGGA backend class accordingly.");
        }
    }
}
// LDA family guard
static void guard_lda_ids_family(const std::vector<int>& ids, int nspin) {
    for (int id : ids) {
        xc_func_type f;
        if (xc_func_init(&f, id, nspin) != 0) {
            throw std::runtime_error("LibXC: could not initialize functional id=" + std::to_string(id));
        }
        const int fam = f.info ? f.info->family : XC_FAMILY_UNKNOWN;
        xc_func_end(&f);
        if (fam != XC_FAMILY_LDA) {
            throw std::runtime_error(
                "LibXC: selected LDA backend but LibXC id " + std::to_string(id) +
                " is not LDA.");
        }
    }
}
// precedence:
//   - HYB-GGA: MRCHEM_LIBXC_IDS_HYB_GGA
//   - GGA    : MRCHEM_LIBXC_IDS_GGA
//   - LDA    : MRCHEM_LIBXC_IDS_LDA
// fallback: MRCHEM_LIBXC_IDS  (filtered to the requested family)
static std::vector<int> pick_ids_for_family(int requested_family, bool /*spin*/) {
    if (requested_family == XC_FAMILY_LDA) {
        auto v = parse_ids_csv(getenv_nn("MRCHEM_LIBXC_IDS_LDA"));
        if (!v.empty()) return v;
    } else if (fam_is_gga_like(requested_family)) {
        auto vhyb = parse_ids_csv(getenv_nn("MRCHEM_LIBXC_IDS_HYB_GGA"));
        if (!vhyb.empty()) return vhyb;
        auto vgga = parse_ids_csv(getenv_nn("MRCHEM_LIBXC_IDS_GGA"));
        if (!vgga.empty())  return vgga;
    }
    auto all = parse_ids_csv(getenv_nn("MRCHEM_LIBXC_IDS"));
    if (all.empty()) return all;
    std::vector<int> filtered;
    for (int id : all) {
        int fam = libxc_family_of_id(id);
        if (requested_family == XC_FAMILY_LDA) {
            if (fam == XC_FAMILY_LDA) filtered.push_back(id);
        } else if (fam_is_gga_like(requested_family)) {
            if (fam_is_gga_like(fam)) filtered.push_back(id);
        }
    }
    return filtered;
}

// ---- NEW: Human-readable LibXC tokens -> numeric IDs ----
static int libxc_id_from_token(const std::string &tok) {
    // 1) try integer form
    try {
        size_t pos = 0;
        long id = std::stol(tok, &pos, 10);
        if (pos == tok.size() && id > 0 && id <= std::numeric_limits<int>::max()) {
            return static_cast<int>(id);
        }
    } catch (...) { /* not an integer */ }
    // 2) LibXC’s name->id lookup
    int id = xc_functional_get_number(tok.c_str());
    if (id <= 0) {
        throw std::runtime_error("LibXC: unknown functional token '" + tok + "'");
    }
    return id;
}
#endif // MRCHEM_ENABLE_LIBXC

Factory::Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA)
   : mra(MRA), xcfun_p(xcfun_new(), xcfun_delete) {}

std::unique_ptr<MRDFT> Factory::build() {
    // Grid
    auto grid_p = std::make_unique<Grid>(mra);

    // Backend selection (env overrides setter)
    std::string backend_eff = to_lower(backend);
    if (const char* be = std::getenv("MRCHEM_XC_BACKEND")) backend_eff = to_lower(std::string(be));
    const bool want_libxc  = (backend_eff == "libxc");

#ifdef MRCHEM_ENABLE_LIBXC
    const bool libxc_available = true;
#else
    const bool libxc_available = false;
#endif
    const bool use_libxc = (want_libxc && libxc_available);

    // XCFun user-mode setup (still used for masks/meta & EXX coefficient)
    const bool gga = xcfun_is_gga(xcfun_p.get());
    const bool lda = !gga;
    unsigned int mode = 1;                // partial derivatives
    unsigned int func_type = gga ? 1u : 0u;
    unsigned int dens_type = 1 + (spin ? 1u : 0u); // 1: unpol, 2: spin
    unsigned int laplacian = 0, kinetic = 0, current = 0;
    unsigned int exp_deriv = gga ? (gamma ? 0u : 1u) : 0u;
    xcfun_user_eval_setup(xcfun_p.get(), order, func_type, dens_type,
                          mode, laplacian, kinetic, current, exp_deriv);

    // Derivative operator for GGA
    if (gga) {
        if      (diff_s == "bspline") diff_p = std::make_unique<mrcpp::BSOperator<3>>(mra, 1);
        else if (diff_s == "abgv_00") diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
        else if (diff_s == "abgv_55") diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.5, 0.5);
    }

    // Build functional
    std::unique_ptr<Functional> func_p;

#ifdef MRCHEM_ENABLE_LIBXC
    // Determine the requested LibXC family for this call
    int req_family = XC_FAMILY_UNKNOWN;
    if (lda) req_family = XC_FAMILY_LDA;
    else     req_family = XC_FAMILY_GGA; // treat HYB-GGA & GGA both through the GGA backend

    // Gather the correct IDs for this call (by family), independent of other calls
    std::vector<int> ids_for_this_call;
    if (use_libxc) {
        // priority: user-setter -> env family vars -> env generic
        ids_for_this_call = libxc_ids;
        if (ids_for_this_call.empty()) {
            ids_for_this_call = pick_ids_for_family(req_family, spin);
        }
    }
#endif

    if (spin) {
        if (gga) {
#ifdef MRCHEM_ENABLE_LIBXC
            if (use_libxc && !ids_for_this_call.empty()) {
                guard_gga_ids_family(ids_for_this_call, XC_POLARIZED);
                func_p = std::make_unique<LibXCSpinGGA>(order, xcfun_p, diff_p, ids_for_this_call);
            } else
#endif
            {
                func_p = std::make_unique<SpinGGA>(order, xcfun_p, diff_p);
            }
        } else { // LDA spin
#ifdef MRCHEM_ENABLE_LIBXC
            if (use_libxc && !ids_for_this_call.empty()) {
                guard_lda_ids_family(ids_for_this_call, XC_POLARIZED);
                func_p = std::make_unique<LibXCSpinLDA>(order, xcfun_p, ids_for_this_call);
            } else
#endif
            {
                func_p = std::make_unique<SpinLDA>(order, xcfun_p);
            }
        }
    } else { // unpolarized
        if (gga) {
#ifdef MRCHEM_ENABLE_LIBXC
            if (use_libxc && !ids_for_this_call.empty()) {
                guard_gga_ids_family(ids_for_this_call, XC_UNPOLARIZED);
                func_p = std::make_unique<LibXCGGA>(order, xcfun_p, diff_p, ids_for_this_call);
            } else
#endif
            {
                func_p = std::make_unique<GGA>(order, xcfun_p, diff_p);
            }
        } else { // LDA unpolarized
#ifdef MRCHEM_ENABLE_LIBXC
            if (use_libxc && !ids_for_this_call.empty()) {
                guard_lda_ids_family(ids_for_this_call, XC_UNPOLARIZED);
                func_p = std::make_unique<LibXCLDA>(order, xcfun_p, ids_for_this_call);
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

    // Finalize (default derivative op upcast to base ptr)
    std::unique_ptr<mrcpp::DerivativeOperator<3>> deriv_base =
        std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
    func_p->setDerivOp(deriv_base);  // setDerivOp will std::move() from it
    func_p->setLogGradient(log_grad);
    func_p->setDensityCutoff(cutoff);

    return std::make_unique<MRDFT>(grid_p, func_p);
}

#ifdef MRCHEM_ENABLE_LIBXC
void Factory::setLibXCTokens(const std::vector<std::string> &tokens) {
    std::vector<int> ids;
    ids.reserve(tokens.size());
    for (const auto &t : tokens) {
        ids.push_back(libxc_id_from_token(t));
    }
    libxc_ids = std::move(ids);
}
#else
void Factory::setLibXCTokens(const std::vector<std::string> &) {
    MSG_ABORT("LibXC token support requested but MRCHEM_ENABLE_LIBXC is not enabled at build time.");
}
#endif

} // namespace mrdft
