// src/mrdft/LibXCBackend.cpp
#include "LibXCBackend.h"

#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>   // std::max

// If xc.h isn't included by the header, uncomment the next line:
// #include <xc.h>

namespace mrdft {

// ======================= helpers (unpolarized) =======================

static void split_gga_inputs_unpolarized(const Eigen::MatrixXd &inp,
                                         Eigen::VectorXd &rho,
                                         Eigen::VectorXd &sigma)
{
    const int n = static_cast<int>(inp.rows());
    const int m = static_cast<int>(inp.cols());

    if (n == 0) { rho.resize(0); sigma.resize(0); return; }

    if (m == 2) {
        // layout: [rho, sigma] where sigma = |∇rho|^2 (or gamma)
        rho   = inp.col(0);
        sigma = inp.col(1);
        return;
    }

    if (m == 4) {
        // layout: [rho, gradx, grady, gradz]  -> sigma = |∇rho|^2
        rho   = inp.col(0);
        sigma.resize(n);
        for (int i = 0; i < n; ++i) {
            const double gx = inp(i, 1);
            const double gy = inp(i, 2);
            const double gz = inp(i, 3);
            sigma(i) = gx*gx + gy*gy + gz*gz;
        }
        return;
    }

    throw std::runtime_error(
        "LibXCGGA: unsupported input layout. "
        "Expected (rho,sigma) or (rho,gradx,grady,gradz). "
        "Got " + std::to_string(m) + " columns.");
}

// Accumulate one LibXC GGA functional into (exc, vrho, vsigma)
// note: LibXC API takes a non-const xc_func_type*; we const_cast safely.
static void add_gga_unpolarized(const xc_func_type &f,
                                int n,
                                const double *rho,   // len n
                                const double *sigma, // len n
                                double *exc,         // len n (accum)
                                double *vrho,        // len n (accum)
                                double *vsigma)      // len n (accum)
{
    std::vector<double> exc_i(n, 0.0), vrho_i(n, 0.0), vsigma_i(n, 0.0);
    xc_func_type *fp = const_cast<xc_func_type*>(&f);
    xc_gga_exc_vxc(fp, n, rho, sigma, exc_i.data(), vrho_i.data(), vsigma_i.data());
    for (int i = 0; i < n; ++i) {
        exc[i]    += exc_i[i];
        vrho[i]   += vrho_i[i];
        vsigma[i] += vsigma_i[i];
    }
}

// ======================= helpers (polarized) =======================

static void split_lda_inputs_polarized(const Eigen::MatrixXd &inp,
                                       Eigen::VectorXd &rho_up,
                                       Eigen::VectorXd &rho_dn)
{
    const int n = static_cast<int>(inp.rows());
    const int m = static_cast<int>(inp.cols());
    if (n == 0) { rho_up.resize(0); rho_dn.resize(0); return; }

    if (m == 2) {           // [rho_up, rho_dn]
        rho_up = inp.col(0);
        rho_dn = inp.col(1);
        return;
    }

    throw std::runtime_error(
        "LibXCSpinLDA: unsupported input layout. "
        "Expected 2 columns (rho_up, rho_dn). Got " + std::to_string(m));
}

static void split_gga_inputs_polarized(const Eigen::MatrixXd &inp,
                                       Eigen::VectorXd &rho_up,
                                       Eigen::VectorXd &rho_dn,
                                       Eigen::VectorXd &sigma_uu,
                                       Eigen::VectorXd &sigma_ud,
                                       Eigen::VectorXd &sigma_dd)
{
    const int n = static_cast<int>(inp.rows());
    const int m = static_cast<int>(inp.cols());
    if (n == 0) {
        rho_up.resize(0); rho_dn.resize(0);
        sigma_uu.resize(0); sigma_ud.resize(0); sigma_dd.resize(0);
        return;
    }

    if (m == 5) {
        // [rho_up, rho_dn, sigma_uu, sigma_ud, sigma_dd]
        rho_up   = inp.col(0);
        rho_dn   = inp.col(1);
        sigma_uu = inp.col(2);
        sigma_ud = inp.col(3);
        sigma_dd = inp.col(4);
        return;
    }

    if (m == 8) {
        // [rho_up, rho_dn, gx_u, gy_u, gz_u, gx_d, gy_d, gz_d]
        rho_up = inp.col(0);
        rho_dn = inp.col(1);
        const Eigen::VectorXd gx_u = inp.col(2), gy_u = inp.col(3), gz_u = inp.col(4);
        const Eigen::VectorXd gx_d = inp.col(5), gy_d = inp.col(6), gz_d = inp.col(7);

        sigma_uu.resize(n);
        sigma_ud.resize(n);
        sigma_dd.resize(n);
        for (int i = 0; i < n; ++i) {
            const double gux = gx_u(i), guy = gy_u(i), guz = gz_u(i);
            const double gdx = gx_d(i), gdy = gy_d(i), gdz = gz_d(i);
            sigma_uu(i) = gux*gux + guy*guy + guz*guz;
            sigma_dd(i) = gdx*gdx + gdy*gdy + gdz*gdz;
            sigma_ud(i) = gux*gdx + guy*gdy + guz*gdz; // can be negative
        }
        return;
    }

    throw std::runtime_error(
        "LibXCSpinGGA: unsupported input layout. "
        "Expected 5 columns (rho_up,rho_dn,sigma_uu,sigma_ud,sigma_dd) "
        "or 8 columns (rho_up,rho_dn,gx_u,gy_u,gz_u,gx_d,gy_d,gz_d). "
        "Got " + std::to_string(m));
}

// ======================= LDA (unpolarized) =======================

Eigen::MatrixXd LibXCLDA::evaluate_transposed(Eigen::MatrixXd &inp) const {
    // Expect a single column: rho
    const int n = static_cast<int>(inp.rows());
    const int m = static_cast<int>(inp.cols());
    if (m != 1) {
        throw std::runtime_error(
            "LibXCLDA: unsupported input layout. Expected 1 column (rho), got " +
            std::to_string(m));
    }

    std::vector<double> rho(n), exc(n, 0.0), vrho(n, 0.0);
    for (int i = 0; i < n; ++i) rho[i] = inp(i, 0);

    for (auto &func : handle.funcs) {
        std::vector<double> exc_i(n, 0.0), vrho_i(n, 0.0);
        xc_lda_exc_vxc(&func, n, rho.data(), exc_i.data(), vrho_i.data());
        for (int i = 0; i < n; ++i) {
            exc[i]  += exc_i[i];
            vrho[i] += vrho_i[i];
        }
    }

    // MRChem expects energy density per volume: F = n * exc
    Eigen::MatrixXd out(n, 2);
    for (int i = 0; i < n; ++i) {
        out(i, 0) = rho[i] * exc[i];  // F
        out(i, 1) = vrho[i];          // dF/dn
    }
    return out;
}

// ======================= GGA (unpolarized) =======================

Eigen::MatrixXd LibXCGGA::evaluate_transposed(Eigen::MatrixXd &inp) const {
    const int n = static_cast<int>(inp.rows());
    const int m = static_cast<int>(inp.cols());

    Eigen::VectorXd rho_v, sigma_v;
    Eigen::VectorXd gx_v, gy_v, gz_v;

    if (m == 2) {
        rho_v   = inp.col(0);
        sigma_v = inp.col(1);                 // σ = |∇ρ|² (or gamma)
    } else if (m == 4) {
        rho_v = inp.col(0);
        gx_v  = inp.col(1);
        gy_v  = inp.col(2);
        gz_v  = inp.col(3);
        sigma_v.resize(n);
        for (int i = 0; i < n; ++i) {
            const double gx = gx_v(i), gy = gy_v(i), gz = gz_v(i);
            sigma_v(i) = gx*gx + gy*gy + gz*gz;
        }
    } else {
        throw std::runtime_error(
            "LibXCGGA: unsupported input layout. "
            "Expected (rho,sigma) or (rho,gradx,grady,gradz). "
            "Got " + std::to_string(m) + " columns.");
    }

    std::vector<double> rho(n), sigma(n), exc(n, 0.0), vrho(n, 0.0), vsigma(n, 0.0);
    for (int i = 0; i < n; ++i) { rho[i] = rho_v(i); sigma[i] = sigma_v(i); }

    for (auto &func : handle.funcs) {
        std::vector<double> exc_i(n, 0.0), vrho_i(n, 0.0), vsigma_i(n, 0.0);
        xc_gga_exc_vxc(&func, n, rho.data(), sigma.data(),
                       exc_i.data(), vrho_i.data(), vsigma_i.data());
        for (int i = 0; i < n; ++i) {
            exc[i]    += exc_i[i];
            vrho[i]   += vrho_i[i];
            vsigma[i] += vsigma_i[i];
        }
    }

    if (m == 2) {
        // sigma layout: [F, dF/dn, dF/dsigma]
        Eigen::MatrixXd out(n, 3);
        for (int i = 0; i < n; ++i) {
            out(i, 0) = rho[i] * exc[i];  // F = n * exc
            out(i, 1) = vrho[i];          // dF/dn
            out(i, 2) = vsigma[i];        // dF/dsigma
        }
        return out;
    } else {
        // gradient layout: [F, dF/dn, dF/dgx, dF/dgy, dF/dgz]
        Eigen::MatrixXd out(n, 5);
        for (int i = 0; i < n; ++i) {
            const double gx = gx_v(i), gy = gy_v(i), gz = gz_v(i);
            out(i, 0) = rho[i] * exc[i];       // F
            out(i, 1) = vrho[i];               // dF/dn
            // chain rule: sigma = gx^2 + gy^2 + gz^2  => dF/dgk = 2*vsigma*gk
            out(i, 2) = 2.0 * vsigma[i] * gx;  // dF/dgx
            out(i, 3) = 2.0 * vsigma[i] * gy;  // dF/dgy
            out(i, 4) = 2.0 * vsigma[i] * gz;  // dF/dgz
        }
        return out;
    }
}

// ======================= Spin LDA (polarized) =======================

Eigen::MatrixXd LibXCSpinLDA::evaluate_transposed(Eigen::MatrixXd &inp) const {
    const int n = static_cast<int>(inp.rows());

    Eigen::VectorXd rho_up_v, rho_dn_v;
    split_lda_inputs_polarized(inp, rho_up_v, rho_dn_v);

    std::vector<double> rho(2*n);
    for (int i = 0; i < n; ++i) {
        rho[2*i + 0] = std::max(0.0, rho_up_v(i));
        rho[2*i + 1] = std::max(0.0, rho_dn_v(i));
    }

    std::vector<double> exc(n, 0.0);
    std::vector<double> vrho(2*n, 0.0);

    for (auto &func : handle.funcs) {
        std::vector<double> exc_i(n, 0.0), vrho_i(2*n, 0.0);
        xc_lda_exc_vxc(&func, n, rho.data(), exc_i.data(), vrho_i.data());
        for (int i = 0; i < n; ++i) {
            exc[i]      += exc_i[i];
            vrho[2*i+0] += vrho_i[2*i+0];
            vrho[2*i+1] += vrho_i[2*i+1];
        }
    }

    // Output: [F, dF/d(rho_up), dF/d(rho_dn)] with F = (rho_up+rho_dn)*exc
    Eigen::MatrixXd out(n, 3);
    for (int i = 0; i < n; ++i) {
        const double n_tot = rho[2*i+0] + rho[2*i+1];
        out(i, 0) = n_tot * exc[i];
        out(i, 1) = vrho[2*i + 0];
        out(i, 2) = vrho[2*i + 1];
    }
    return out;
}

// ======================= Spin GGA (polarized) =======================

Eigen::MatrixXd LibXCSpinGGA::evaluate_transposed(Eigen::MatrixXd &inp) const {
    const int n = static_cast<int>(inp.rows());
    const int m = static_cast<int>(inp.cols());

    Eigen::VectorXd rho_up_v, rho_dn_v, sigma_uu_v, sigma_ud_v, sigma_dd_v;
    split_gga_inputs_polarized(inp, rho_up_v, rho_dn_v, sigma_uu_v, sigma_ud_v, sigma_dd_v);

    // Flattened LibXC layout:
    //   rho[2*i+0]=rho_up, rho[2*i+1]=rho_dn
    //   sigma[3*i+0]=sigma_uu, [3*i+1]=sigma_ud, [3*i+2]=sigma_dd
    std::vector<double> rho(2*n);
    std::vector<double> sigma(3*n);
    for (int i = 0; i < n; ++i) {
        rho[2*i + 0] = std::max(0.0, rho_up_v(i));
        rho[2*i + 1] = std::max(0.0, rho_dn_v(i));
        sigma[3*i + 0] = std::max(0.0, sigma_uu_v(i));
        sigma[3*i + 1] =           sigma_ud_v(i);  // can be negative
        sigma[3*i + 2] = std::max(0.0, sigma_dd_v(i));
    }

    std::vector<double> exc(n, 0.0);
    std::vector<double> vrho(2*n, 0.0);
    std::vector<double> vsigma(3*n, 0.0);

    for (auto &func : handle.funcs) {
        std::vector<double> exc_i(n, 0.0), vrho_i(2*n, 0.0), vsigma_i(3*n, 0.0);
        xc_gga_exc_vxc(&func, n, rho.data(), sigma.data(),
                       exc_i.data(), vrho_i.data(), vsigma_i.data());
        for (int i = 0; i < n; ++i) {
            exc[i]        += exc_i[i];
            vrho[2*i+0]   += vrho_i[2*i+0];
            vrho[2*i+1]   += vrho_i[2*i+1];
            vsigma[3*i+0] += vsigma_i[3*i+0];
            vsigma[3*i+1] += vsigma_i[3*i+1];
            vsigma[3*i+2] += vsigma_i[3*i+2];
        }
    }

    // If caller provided gradients (8-col input), expose gradient derivatives via chain rule:
    //   dF/d(∇ρ↑) = 2*vsigma_uu*∇ρ↑ + vsigma_ud*∇ρ↓
    //   dF/d(∇ρ↓) = 2*vsigma_dd*∇ρ↓ + vsigma_ud*∇ρ↑
    if (m == 8) {
        const Eigen::VectorXd gx_u = inp.col(2), gy_u = inp.col(3), gz_u = inp.col(4);
        const Eigen::VectorXd gx_d = inp.col(5), gy_d = inp.col(6), gz_d = inp.col(7);

        Eigen::MatrixXd out(n, 9);
        for (int i = 0; i < n; ++i) {
            const double n_tot = rho[2*i+0] + rho[2*i+1];
            const double vuu = vsigma[3*i+0];
            const double vud = vsigma[3*i+1];
            const double vdd = vsigma[3*i+2];

            // dF/d(grad up)
            const double dFx_u = 2.0*vuu*gx_u(i) + vud*gx_d(i);
            const double dFy_u = 2.0*vuu*gy_u(i) + vud*gy_d(i);
            const double dFz_u = 2.0*vuu*gz_u(i) + vud*gz_d(i);
            // dF/d(grad down)
            const double dFx_d = 2.0*vdd*gx_d(i) + vud*gx_u(i);
            const double dFy_d = 2.0*vdd*gy_d(i) + vud*gy_u(i);
            const double dFz_d = 2.0*vdd*gz_d(i) + vud*gz_u(i);

            out(i, 0) = n_tot * exc[i];     // F
            out(i, 1) = vrho[2*i + 0];      // dF/d rho_up
            out(i, 2) = vrho[2*i + 1];      // dF/d rho_dn
            out(i, 3) = dFx_u;              // dF/d gx_up
            out(i, 4) = dFy_u;              // dF/d gy_up
            out(i, 5) = dFz_u;              // dF/d gz_up
            out(i, 6) = dFx_d;              // dF/d gx_dn
            out(i, 7) = dFy_d;              // dF/d gy_dn
            out(i, 8) = dFz_d;              // dF/d gz_dn
        }
        return out;
    }

    // Otherwise (5-col sigma layout): [F, dF/drho_up, dF/drho_dn, dF/dsigma_uu, dF/dsigma_ud, dF/dsigma_dd]
    Eigen::MatrixXd out(n, 6);
    for (int i = 0; i < n; ++i) {
        const double n_tot = rho[2*i+0] + rho[2*i+1];
        out(i, 0) = n_tot * exc[i];
        out(i, 1) = vrho[2*i + 0];
        out(i, 2) = vrho[2*i + 1];
        out(i, 3) = vsigma[3*i + 0];
        out(i, 4) = vsigma[3*i + 1];
        out(i, 5) = vsigma[3*i + 2];
    }
    return out;
}

} // namespace mrdft
