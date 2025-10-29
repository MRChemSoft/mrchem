// MRChem — LibXC backend (batched + thread-local)
#include "LibXCBackend.h"

#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdlib>

namespace mrdft {

// =================== small helpers ===================

static inline double finite_or_zero(double x) { return std::isfinite(x) ? x : 0.0; }

static inline void assert_gga_like(const xc_func_type& f) {
    if (!f.info) {
        throw std::runtime_error("LibXC: functional has null info pointer");
    }
    const int fam = f.info->family;
    if (fam != XC_FAMILY_GGA && fam != XC_FAMILY_HYB_GGA) {
        throw std::runtime_error(
            std::string("LibXC: functional family mismatch in GGA backend. Got family=") +
            std::to_string(fam) + " (expected GGA or HYB-GGA). "
            "Use an LDA/MGGA backend class for other families.");
    }
}

static double read_env_thresh() {
    if (const char* v = std::getenv("MRCHEM_LIBXC_DENS_THRESH")) {
        try {
            double t = std::stod(v);
            return (t > 0.0 ? t : 0.0);
        } catch (...) { /* ignore */ }
    }
    return 1e-12; // sensible default
}

// Thread-local LibXC context (one per OpenMP thread)
struct XCTLS {
    std::vector<xc_func_type> funcs;
    std::vector<int> ids;
    int nspin = -1;
    double dens_thresh = -1.0;

    void clear() {
        for (auto &f : funcs) xc_func_end(&f);
        funcs.clear();
        ids.clear();
        nspin = -1;
        dens_thresh = -1.0;
    }

    void ensure(const std::vector<int>& want_ids, int want_nspin, double want_thresh) {
        bool same_ids = (ids.size() == want_ids.size()) &&
                        std::equal(ids.begin(), ids.end(), want_ids.begin());
        if (!same_ids || nspin != want_nspin || dens_thresh != want_thresh) {
            clear();
            ids = want_ids;
            nspin = want_nspin;
            dens_thresh = want_thresh;

            funcs.resize(ids.size());
            for (size_t i = 0; i < ids.size(); ++i) {
                if (xc_func_init(&funcs[i], ids[i], nspin) != 0) {
                    throw std::runtime_error("LibXC: could not initialize functional id=" + std::to_string(ids[i]));
                }
                if (dens_thresh > 0.0) {
                    xc_func_set_dens_threshold(&funcs[i], dens_thresh);
                }
            }
        }
    }

    ~XCTLS() { clear(); }
};

static thread_local XCTLS xc_tls;

// =====================================================
// =============== LDA (unpolarized) ===================
// =====================================================

Eigen::MatrixXd LibXCLDA::eval_lda_transposed(Eigen::MatrixXd &inp) const {
    const int n = static_cast<int>(inp.rows());
    const int m = static_cast<int>(inp.cols());
    if (m != 1) {
        throw std::runtime_error("LibXCLDA: unsupported input layout. Expected 1 column (rho), got " +
                                 std::to_string(m));
    }
    if (n == 0) return Eigen::MatrixXd(0, 2);

    // Prepare arrays
    static thread_local std::vector<double> rho, eps, vrho, tmp_eps, tmp_vrho;
    rho.resize(n);
    eps.assign(n, 0.0);
    vrho.assign(n, 0.0);
    tmp_eps.resize(n);
    tmp_vrho.resize(n);

    for (int i = 0; i < n; ++i) rho[i] = finite_or_zero(inp(i, 0));

    // Ensure TLS contexts
    xc_tls.ensure(ids_, nspin_, read_env_thresh());

    // Accumulate all functionals
    for (auto &f : xc_tls.funcs) {
        xc_lda_exc_vxc(&f, n, rho.data(), tmp_eps.data(), tmp_vrho.data());
        for (int i = 0; i < n; ++i) {
            eps[i]  += tmp_eps[i];
            vrho[i] += tmp_vrho[i];
        }
    }

    // Output: [F, dF/dn] with F = n * εxc
    Eigen::MatrixXd out(n, 2);
    for (int i = 0; i < n; ++i) {
        out(i, 0) = rho[i] * eps[i];
        out(i, 1) = vrho[i];
    }
    return out;
}

// =====================================================
// ================ GGA (unpolarized) ==================
// =====================================================

Eigen::MatrixXd LibXCGGA::eval_gga_transposed(Eigen::MatrixXd &inp) const {
    const int n = static_cast<int>(inp.rows());
    const int m = static_cast<int>(inp.cols());
    if (n == 0) {
        return (m == 2 ? Eigen::MatrixXd(0, 3) : Eigen::MatrixXd(0, 5));
    }
    if (m != 2 && m != 4) {
        throw std::runtime_error(
            "LibXCGGA: unsupported input layout. "
            "Expected (rho,sigma) or (rho,gradx,grady,gradz). "
            "Got " + std::to_string(m) + " columns.");
    }

    // Prepare arrays
    static thread_local std::vector<double> rho, sigma, gx, gy, gz;
    static thread_local std::vector<double> eps, vrho, vsigma;
    static thread_local std::vector<double> teps, tvrho, tvsigma;

    rho.resize(n);
    if (m == 2) {
        sigma.resize(n);
    } else {
        gx.resize(n); gy.resize(n); gz.resize(n);
        sigma.resize(n);
    }
    eps.assign(n, 0.0);
    vrho.assign(n, 0.0);
    vsigma.assign(n, 0.0);
    teps.resize(n);
    tvrho.resize(n);
    tvsigma.resize(n);

    // Fill inputs
    if (m == 2) {
        for (int i = 0; i < n; ++i) {
            rho[i]   = finite_or_zero(inp(i, 0));
            sigma[i] = std::max(0.0, finite_or_zero(inp(i, 1)));
        }
    } else { // m == 4
        for (int i = 0; i < n; ++i) {
            rho[i] = finite_or_zero(inp(i, 0));
            gx[i]  = finite_or_zero(inp(i, 1));
            gy[i]  = finite_or_zero(inp(i, 2));
            gz[i]  = finite_or_zero(inp(i, 3));
            sigma[i] = std::max(0.0, gx[i]*gx[i] + gy[i]*gy[i] + gz[i]*gz[i]);
        }
    }

    // TLS contexts
    xc_tls.ensure(ids_, nspin_, read_env_thresh());
    for (auto &f : xc_tls.funcs) {
        assert_gga_like(f);
        xc_gga_exc_vxc(&f, n, rho.data(), sigma.data(),
                       teps.data(), tvrho.data(), tvsigma.data());
        for (int i = 0; i < n; ++i) {
            eps[i]    += teps[i];
            vrho[i]   += tvrho[i];
            vsigma[i] += tvsigma[i];
        }
    }

    if (m == 2) {
        // sigma-layout output: [F, dF/dn, dF/dsigma]
        Eigen::MatrixXd out(n, 3);
        for (int i = 0; i < n; ++i) {
            out(i, 0) = rho[i] * eps[i];
            out(i, 1) = vrho[i];
            out(i, 2) = vsigma[i];
        }
        return out;
    } else {
        // gradient layout output: [F, dF/dn, dF/dgx, dF/dgy, dF/dgz]
        Eigen::MatrixXd out(n, 5);
        for (int i = 0; i < n; ++i) {
            const double s2 = 2.0 * vsigma[i];
            out(i, 0) = rho[i] * eps[i];
            out(i, 1) = vrho[i];
            out(i, 2) = s2 * gx[i];
            out(i, 3) = s2 * gy[i];
            out(i, 4) = s2 * gz[i];
        }
        return out;
    }
}

// =====================================================
// ============== Spin LDA (polarized) =================
// =====================================================

Eigen::MatrixXd LibXCSpinLDA::eval_lda_transposed(Eigen::MatrixXd &inp) const {
    const int n = static_cast<int>(inp.rows());
    const int m = static_cast<int>(inp.cols());
    if (m != 2) {
        throw std::runtime_error(
            "LibXCSpinLDA: unsupported input layout. "
            "Expected 2 columns (rho_up, rho_dn). Got " + std::to_string(m));
    }
    if (n == 0) return Eigen::MatrixXd(0, 3);

    // Prepare arrays
    static thread_local std::vector<double> rho2, eps, vrho2, teps, tvrho2;
    rho2.resize(2*n);
    eps.assign(n, 0.0);
    vrho2.assign(2*n, 0.0);
    teps.resize(n);
    tvrho2.resize(2*n);

    for (int i = 0; i < n; ++i) {
        const double up_raw = inp(i, 0);
        const double dn_raw = inp(i, 1);
        const double up = std::max(0.0, finite_or_zero(up_raw));
        const double dn = std::max(0.0, finite_or_zero(dn_raw));
        rho2[2*i + 0] = up;
        rho2[2*i + 1] = dn; 
    }

    xc_tls.ensure(ids_, nspin_, read_env_thresh());

    for (auto &f : xc_tls.funcs) {
        xc_lda_exc_vxc(&f, n, rho2.data(), teps.data(), tvrho2.data());
        for (int i = 0; i < n; ++i) {
            eps[i]         += teps[i];
            vrho2[2*i + 0] += tvrho2[2*i + 0];
            vrho2[2*i + 1] += tvrho2[2*i + 1];
        }
    }

    // Output: [F, dF/d rho_up, dF/d rho_dn] with F = (rho_up+rho_dn)*εxc
    Eigen::MatrixXd out(n, 3);
    for (int i = 0; i < n; ++i) {
        const double up = rho2[2*i + 0];
        const double dn = rho2[2*i + 1];
        out(i, 0) = (up + dn) * eps[i];
        out(i, 1) = vrho2[2*i + 0];
        out(i, 2) = vrho2[2*i + 1];
    }
    return out;
}

// =====================================================
// ============== Spin GGA (polarized) =================
// =====================================================

Eigen::MatrixXd LibXCSpinGGA::eval_gga_transposed(Eigen::MatrixXd &inp) const {
    const int n = static_cast<int>(inp.rows());
    const int m = static_cast<int>(inp.cols());
    if (n == 0) {
        return (m == 5 ? Eigen::MatrixXd(0, 6)
                       : (m == 8 ? Eigen::MatrixXd(0, 9)
                                 : Eigen::MatrixXd(0, 0)));
    }
    if (m != 5 && m != 8) {
        throw std::runtime_error(
            "LibXCSpinGGA: unsupported input layout. "
            "Expected 5 columns (rho_up,rho_dn,sigma_uu,sigma_ud,sigma_dd) "
            "or 8 columns (rho_up,rho_dn,gx_u,gy_u,gz_u,gx_d,gy_d,gz_d). "
            "Got " + std::to_string(m));
    }

    // Prepare arrays
    static thread_local std::vector<double> rho2;      // [up, dn] per point
    static thread_local std::vector<double> sigma3;    // [uu, ud, dd] per point
    static thread_local std::vector<double> gax, gay, gaz, gbx, gby, gbz;

    static thread_local std::vector<double> eps, vrho2, vs3;
    static thread_local std::vector<double> teps, tvrho2, tvs3;

    rho2.resize(2*n);
    sigma3.resize(3*n);

    if (m == 8) {
        gax.resize(n); gay.resize(n); gaz.resize(n);
        gbx.resize(n); gby.resize(n); gbz.resize(n);
    }

    eps.assign(n, 0.0);
    vrho2.assign(2*n, 0.0);
    vs3.assign(3*n, 0.0);
    teps.resize(n);
    tvrho2.resize(2*n);
    tvs3.resize(3*n);

    if (m == 5) {
        for (int i = 0; i < n; ++i) {
            const double up = std::max(0.0, finite_or_zero(inp(i, 0)));
            const double dn = std::max(0.0, finite_or_zero(inp(i, 1)));
            rho2[2*i + 0] = up;
            rho2[2*i + 1] = dn;

            sigma3[3*i + 0] = std::max(0.0, finite_or_zero(inp(i, 2))); // uu
            sigma3[3*i + 1] =              finite_or_zero(inp(i, 3));   // ud (can be negative)
            sigma3[3*i + 2] = std::max(0.0, finite_or_zero(inp(i, 4))); // dd
        }
    } else { // m == 8 (gradients)
        for (int i = 0; i < n; ++i) {
            const double up = std::max(0.0, finite_or_zero(inp(i, 0)));
            const double dn = std::max(0.0, finite_or_zero(inp(i, 1)));
            rho2[2*i + 0] = up;
            rho2[2*i + 1] = dn;

            gax[i] = finite_or_zero(inp(i, 2));
            gay[i] = finite_or_zero(inp(i, 3));
            gaz[i] = finite_or_zero(inp(i, 4));
            gbx[i] = finite_or_zero(inp(i, 5));
            gby[i] = finite_or_zero(inp(i, 6));
            gbz[i] = finite_or_zero(inp(i, 7));

            sigma3[3*i + 0] = std::max(0.0, gax[i]*gax[i] + gay[i]*gay[i] + gaz[i]*gaz[i]); // uu
            sigma3[3*i + 1] = gax[i]*gbx[i] + gay[i]*gby[i] + gaz[i]*gbz[i];                // ud
            sigma3[3*i + 2] = std::max(0.0, gbx[i]*gbx[i] + gby[i]*gby[i] + gbz[i]*gbz[i]); // dd
        }
    }

    xc_tls.ensure(ids_, nspin_, read_env_thresh());
    for (auto &f : xc_tls.funcs) {
        assert_gga_like(f);
        xc_gga_exc_vxc(&f, n, rho2.data(), sigma3.data(),
                       teps.data(), tvrho2.data(), tvs3.data());
        for (int i = 0; i < n; ++i) {
            eps[i]         += teps[i];
            vrho2[2*i + 0] += tvrho2[2*i + 0];
            vrho2[2*i + 1] += tvrho2[2*i + 1];
            vs3[3*i + 0]   += tvs3[3*i + 0];
            vs3[3*i + 1]   += tvs3[3*i + 1];
            vs3[3*i + 2]   += tvs3[3*i + 2];
        }
    }

    if (m == 8) {
        // gradient form output:
        // [F, dF/drho_up, dF/drho_dn, dF/dgx_up, dF/dgy_up, dF/dgz_up, dF/dgx_dn, dF/dgy_dn, dF/dgz_dn]
        Eigen::MatrixXd out(n, 9);
        for (int i = 0; i < n; ++i) {
            const double up = rho2[2*i + 0];
            const double dn = rho2[2*i + 1];

            const double vuu = vs3[3*i + 0];
            const double vud = vs3[3*i + 1];
            const double vdd = vs3[3*i + 2];

            // ∂F/∂∇ρa = 2*vuu*∇ρa + vud*∇ρb
            const double dFax = 2.0*vuu*gax[i] + vud*gbx[i];
            const double dFay = 2.0*vuu*gay[i] + vud*gby[i];
            const double dFaz = 2.0*vuu*gaz[i] + vud*gbz[i];
            // ∂F/∂∇ρb = 2.0*vdd*∇ρb + vud*∇ρa
            const double dFbx = 2.0*vdd*gbx[i] + vud*gax[i];
            const double dFby = 2.0*vdd*gby[i] + vud*gay[i];
            const double dFbz = 2.0*vdd*gbz[i] + vud*gaz[i];

            out(i, 0) = (up + dn) * eps[i];
            out(i, 1) = vrho2[2*i + 0];
            out(i, 2) = vrho2[2*i + 1];
            out(i, 3) = dFax;
            out(i, 4) = dFay;
            out(i, 5) = dFaz;
            out(i, 6) = dFbx;
            out(i, 7) = dFby;
            out(i, 8) = dFbz;
        }
        return out;
    } else { // m == 5 (sigma layout)
        // [F, dF/drho_up, dF/drho_dn, dF/dsigma_uu, dF/dsigma_ud, dF/dsigma_dd]
        Eigen::MatrixXd out(n, 6);
        for (int i = 0; i < n; ++i) {
            const double up = rho2[2*i + 0];
            const double dn = rho2[2*i + 1];
            out(i, 0) = (up + dn) * eps[i];
            out(i, 1) = vrho2[2*i + 0];
            out(i, 2) = vrho2[2*i + 1];
            out(i, 3) = vs3[3*i + 0];
            out(i, 4) = vs3[3*i + 1];
            out(i, 5) = vs3[3*i + 2];
        }
        return out;
    }
}

} // namespace mrdft
