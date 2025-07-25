/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

 #include "LibXC.h"
 #include <MRCPP/Printer>
 #include <algorithm>
 #include <map>
 
//  namespace mrchem {
 namespace mrdft {
 
 LibXC::LibXC()
     : order(1)
     , spin_polarized(false)
     , use_gamma(false)
     , has_lda(false)
     , has_gga(false)
     , has_mgga(false)
     , exx_coef(0.0) {
 }
 
 
 bool LibXC::setFunctional(const std::string &name, double weight) {
     // Map common XCFun functional names to LibXC identifiers
     static const std::map<std::string, int> func_map = {
          // LDA functionals
          {"slater", XC_LDA_X},
          {"vwn", XC_LDA_C_VWN},
          {"vwn5", XC_LDA_C_VWN},
          {"vwn3", XC_LDA_C_VWN_3},
          {"pw92c", XC_LDA_C_PW},          
          // GGA functionals
          {"becke88", XC_GGA_X_B88},
          {"lyp", XC_GGA_C_LYP},
          {"pbex", XC_GGA_X_PBE},
          {"pbec", XC_GGA_C_PBE},
          {"revpbex", XC_GGA_X_PBE_R},         
          // Hybrid functionals
          {"b3lyp", XC_HYB_GGA_XC_B3LYP},
          {"pbe0", XC_HYB_GGA_XC_APBE0},           
          // Meta-GGA functionals
          {"tpss", XC_MGGA_X_TPSS},
          {"m06l", XC_MGGA_X_M06_L}
     };
     
     auto it = func_map.find(name);
     if (it != func_map.end()) {
         return setFunctional(it->second, weight);
     } else {
         // Try to parse as a numeric ID
         try {
             int func_id = std::stoi(name);
             return setFunctional(func_id, weight);
         } catch (...) {
             MSG_ERROR("Unknown functional: " + name);
             return false;
         }
     }
 }
 
 bool LibXC::setFunctional(int func_id, double weight) {
     FunctionalUnit unit;
     unit.weight = weight;
     
     // Initialize with unpolarized by default, will be reinitialized during evaluation
     int nspin = spin_polarized ? XC_POLARIZED : XC_UNPOLARIZED;
     if (xc_func_init(&unit.func, func_id, nspin) != 0) {
         MSG_ERROR("Error initializing LibXC functional with ID: " + std::to_string(func_id));
         return false;
     }
     
     unit.initialized = true;
     unit.family = unit.func.info->family;
     
     // Track functional types
     if (unit.family == XC_FAMILY_LDA) has_lda = true;
     else if (unit.family == XC_FAMILY_GGA) has_gga = true;
     else if (unit.family == XC_FAMILY_MGGA) has_mgga = true;
     else if (unit.family == XC_FAMILY_HYB_GGA) {
         has_gga = true;
         double coef = xc_hyb_exx_coef(&unit.func);
         exx_coef += weight * coef;
     }
     else if (unit.family == XC_FAMILY_HYB_MGGA) {
         has_mgga = true;
         double coef = xc_hyb_exx_coef(&unit.func);
         exx_coef += weight * coef;
     }
     
     functionals.push_back(std::move(unit));
     return true;
 }
 
 void LibXC::setOrder(int k) {
     order = k;
 }
 
 void LibXC::setSpin(bool spin) {
     spin_polarized = spin;
 }
 
 void LibXC::setInputMode(bool gamma) {
     use_gamma = gamma;
 }
 
 int LibXC::getInputLength() const {
     if (!spin_polarized) {
         // Non-spin-polarized case
         if (isLDA()) return 1;       // rho
         if (isGGA() && use_gamma) return 2;  // rho, gamma
         if (isGGA()) return 4;       // rho, grad_x, grad_y, grad_z
         if (isMetaGGA()) return 5;   // rho, grad_x, grad_y, grad_z, tau
     } else {
         // Spin-polarized case
         if (isLDA()) return 2;       // rho_a, rho_b
         if (isGGA() && use_gamma) return 5;  // rho_a, rho_b, gamma_aa, gamma_ab, gamma_bb
         if (isGGA()) return 8;       // rho_a, rho_b, grad_x_a, grad_y_a, grad_z_a, grad_x_b, grad_y_b, grad_z_b
         if (isMetaGGA()) return 10;  // rho_a, rho_b, grad_x_a, ..., grad_z_b, tau_a, tau_b
     }
     return 0;
 }
 
 int LibXC::getOutputLength() const {
     if (order == 0) return 1;  // Just energy
     
     // First derivatives (order = 1)
     if (!spin_polarized) {
         if (isLDA()) return 2;      // eps, vrho
         if (isGGA() && use_gamma) return 3;  // eps, vrho, vgamma
         if (isGGA()) return 5;      // eps, vrho, vgrad_x, vgrad_y, vgrad_z
         if (isMetaGGA()) return 6;  // eps, vrho, vgrad_x, vgrad_y, vgrad_z, vtau
     } else {
         if (isLDA()) return 3;      // eps, vrho_a, vrho_b
         if (isGGA() && use_gamma) return 6;  // eps, vrho_a, vrho_b, vgamma_aa, vgamma_ab, vgamma_bb
         if (isGGA()) return 9;      // eps, vrho_a, vrho_b, vgrad_x_a, ..., vgrad_z_b
         if (isMetaGGA()) return 11; // eps, vrho_a, vrho_b, vgrad_x_a, ..., vgrad_z_b, vtau_a, vtau_b
     }
     
     // Higher order derivatives (order > 1) - implementation would depend on MRChem's needs
     // This is a placeholder and would need to be implemented based on specific requirements
     
     return 0;
 }
 
 void LibXC::eval(const double *rho_in, double *result_out) {
     // Clear the output array
     int outlen = getOutputLength();
     std::fill(result_out, result_out + outlen, 0.0);
     
     // Determine functional type and call appropriate evaluation method
     if (has_mgga) {
         // MetaGGA case - not implemented in this example
         // Would need to extract rho, sigma, tau from rho_in
         // evaluateMetaGGA(rho, sigma, lapl, tau, result_out);
     } else if (has_gga) {
         if (use_gamma) {
             // GGA with gamma-type derivatives
             const double *rho, *gamma;
             if (!spin_polarized) {
                 rho = rho_in;
                 gamma = rho_in + 1;
             } else {
                 rho = rho_in;
                 gamma = rho_in + 2;
             }
             
             // Convert gamma to sigma format required by LibXC
             double sigma[3];
             if (!spin_polarized) {
                 sigma[0] = gamma[0];  // |∇ρ|²
             } else {
                 sigma[0] = gamma[0];  // |∇ρₐ|²
                 sigma[1] = gamma[1];  // ∇ρₐ·∇ρᵦ
                 sigma[2] = gamma[2];  // |∇ρᵦ|²
             }
             
             evaluateGGA(rho, sigma, result_out);
         } else {
             // GGA with explicit derivatives
             const double *rho, *grad;
             double sigma[3] = {0.0, 0.0, 0.0};
             
             if (!spin_polarized) {
                 rho = rho_in;
                 grad = rho_in + 1;
                 
                 // Calculate sigma from gradients
                 sigma[0] = grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2];
             } else {
                 rho = rho_in;
                 double *grad_a = const_cast<double*>(rho_in + 2);
                 double *grad_b = const_cast<double*>(rho_in + 5);
                 
                 // Calculate sigma from gradients
                 sigma[0] = grad_a[0]*grad_a[0] + grad_a[1]*grad_a[1] + grad_a[2]*grad_a[2];
                 sigma[1] = grad_a[0]*grad_b[0] + grad_a[1]*grad_b[1] + grad_a[2]*grad_b[2];
                 sigma[2] = grad_b[0]*grad_b[0] + grad_b[1]*grad_b[1] + grad_b[2]*grad_b[2];
             }
             
             evaluateGGA(rho, sigma, result_out);
         }
     } else if (has_lda) {
         // LDA case
         evaluateLDA(rho_in, result_out);
     }
 }
 
 void LibXC::evaluateLDA(const double *rho, double *out) {
     // Temporary arrays for results from individual functionals
     std::vector<double> exc(1, 0.0);
     std::vector<double> vxc;
     
     if (order >= 1) {
         vxc.resize(spin_polarized ? 2 : 1, 0.0);
     }
     
     for (auto &func : functionals) {
         if (func.family != XC_FAMILY_LDA) continue;
         
         // Reinitialize with correct spin
         xc_func_end(&func.func);
         int nspin = spin_polarized ? XC_POLARIZED : XC_UNPOLARIZED;
         xc_func_init(&func.func, func.func.info->number, nspin);
         
         if (order == 0) {
             // Energy only
             xc_lda_exc(&func.func, 1, rho, exc.data());//What does this line do?
             out[0] += func.weight * exc[0];
         } else {
             // Energy and potential
             xc_lda_exc_vxc(&func.func, 1, rho, exc.data(), vxc.data());
             out[0] += func.weight * exc[0]; // What is the order here?
             
             if (spin_polarized) {
                 out[1] += func.weight * vxc[0];  // vrho_a
                 out[2] += func.weight * vxc[1];  // vrho_b
             } else {
                 out[1] += func.weight * vxc[0];  // vrho
             }
         }
     }
 }
 
 void LibXC::evaluateGGA(const double *rho, const double *sigma, double *out) {
     // Temporary arrays for results from individual functionals
     std::vector<double> exc(1, 0.0);
     std::vector<double> vrho, vsigma;
     
     if (order >= 1) {
         vrho.resize(spin_polarized ? 2 : 1, 0.0);
         vsigma.resize(spin_polarized ? 3 : 1, 0.0);
     }
     
     for (auto &func : functionals) {
         if (func.family != XC_FAMILY_GGA && func.family != XC_FAMILY_HYB_GGA) continue;
         
         // Reinitialize with correct spin
         xc_func_end(&func.func);
         int nspin = spin_polarized ? XC_POLARIZED : XC_UNPOLARIZED;
         xc_func_init(&func.func, func.func.info->number, nspin);
         
         if (order == 0) {
             // Energy only
             xc_gga_exc(&func.func, 1, rho, sigma, exc.data());
             out[0] += func.weight * exc[0];
         } else {
             // Energy and potential
             xc_gga_exc_vxc(&func.func, 1, rho, sigma, exc.data(), vrho.data(), vsigma.data());
             out[0] += func.weight * exc[0];
             
             if (spin_polarized) {
                 out[1] += func.weight * vrho[0];    // vrho_a
                 out[2] += func.weight * vrho[1];    // vrho_b
                 
                 if (use_gamma) {
                     // Gamma-type derivatives
                     out[3] += func.weight * vsigma[0];  // vgamma_aa
                     out[4] += func.weight * vsigma[1];  // vgamma_ab
                     out[5] += func.weight * vsigma[2];  // vgamma_bb
                 } else {
                     // Need to transform vsigma to gradient derivatives
                     // (Not implemented here - would require additional tensor algebra)
                 }
             } else {
                 out[1] += func.weight * vrho[0];    // vrho
                 
                 if (use_gamma) {
                     // Gamma-type derivatives
                     out[2] += func.weight * vsigma[0];  // vgamma
                 } else {
                     // Need to transform vsigma to gradient derivatives
                     // (Not implemented here - would require additional tensor algebra)
                 }
             }
         }
     }
 }
 
 void LibXC::evaluateMetaGGA(const double *rho, const double *sigma, 
                           const double *lapl, const double *tau, double *out) {
     // Not implemented in this example
     // Would follow similar pattern to LDA and GGA cases
     // using xc_mgga_exc() and xc_mgga_exc_vxc()
 }
 
 bool LibXC::isLDA() const {
     return has_lda && !has_gga && !has_mgga;
 }
 
 bool LibXC::isGGA() const {
     return has_gga && !has_mgga;
 }
 
 bool LibXC::isMetaGGA() const {
     return has_mgga;
 }
 
 int LibXC::getFunctionalLevel() const {
     if (has_mgga) return 2;  // Meta-GGA
     if (has_gga) return 1;   // GGA
     if (has_lda) return 0;   // LDA
     return -1;               // Unknown/error
 }
 
 } // namespace mrdft
//  } // namespace mrchem
