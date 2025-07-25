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

 #pragma once

 #include <MRCPP/trees/FunctionTreeVector.h>

 #include <memory>
 #include <string>
 #include <vector>
 #include <xc.h>
 #include <xc_funcs.h>

 
 
//  namespace mrchem {
 namespace mrdft {
 
 // Wrapper class for LibXC functionality
 class LibXC {
 public:
     LibXC();
    //  ~LibXC();
     virtual ~LibXC() = default;
 
     // Set up functional by name or ID
     bool setFunctional(const std::string &name, double weight = 1.0);
     bool setFunctional(int func_id, double weight = 1.0);
     
     // Evaluate functional and derivatives
     void eval(const double *rho_in, double *result_out);
     
     // Configuration methods
     void setOrder(int k);
     void setSpin(bool spin_polarized);
     void setInputMode(bool use_gamma);
     
     // Query methods
    //  int getInputLength() const;
    //  int getOutputLength() const;
     virtual int getCtrInputLength() const = 0;
     virtual int getCtrOutputLength() const = 0;
     double getEXX() const { return exx_coef; }

    virtual int getInputLength() const = 0;
    virtual int getOutputLength() const = 0;
     
     // Type checking methods
    //  bool isSpin() const { return spin_polarized; }
     virtual bool isSpin() const = 0;
     bool isLDA() const;
     bool isGGA() const;
     bool isMetaGGA() const;
     bool isHybrid() const { return exx_coef > 1.0e-10; }

    virtual void clear() = 0;
    virtual mrcpp::FunctionTreeVector<3> setupXCInput() = 0;
    virtual mrcpp::FunctionTreeVector<3> setupCtrInput() = 0;
    virtual void preprocess(mrcpp::FunctionTreeVector<3> &inp) = 0;
    virtual mrcpp::FunctionTreeVector<3> postprocess(mrcpp::FunctionTreeVector<3> &inp) = 0;

     
 private:
     struct FunctionalUnit {
         xc_func_type func;
         double weight;
         int family;
         bool initialized;
         
         FunctionalUnit() : weight(0.0), family(XC_FAMILY_UNKNOWN), initialized(false) {}
         ~FunctionalUnit() { if (initialized) xc_func_end(&func); }
     };
     
     std::vector<FunctionalUnit> functionals;
     
     // Configuration state
     int order;              // Derivative order 
     bool spin_polarized;    // Spin-polarized or unpolarized
     bool use_gamma;         // For GGA: use gamma or explicit gradient
     
     // Internal state
     bool has_lda;
     bool has_gga;
     bool has_mgga;
     double exx_coef;        // Exact exchange coefficient for hybrid functionals
     
     // Internal evaluation helpers
     void evaluateLDA(const double *rho, double *out);
     void evaluateGGA(const double *rho, const double *sigma, double *out);
     void evaluateMetaGGA(const double *rho, const double *sigma, 
                         const double *lapl, const double *tau, double *out);
                       
     // Helper method to identify highest functional level
     int getFunctionalLevel() const;
 };
 
 } // namespace mrdft
//  } // namespace mrchem
