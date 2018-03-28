#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#pragma GCC system_header
#include <Eigen/Eigenvalues>

#include <vector>
#include <string>

#include "QMTensorOperator.h"

class Getkw;

class Nuclei;
class Molecule;
class OrbitalVector;
class OrbitalOptimizer;
class EnergyOptimizer;
class GroundStateSolver;
class LinearResponseSolver;
class HelmholtzOperatorSet;
class KAIN;

class PoissonOperator;
template<int D> class PHOperator;
template<int D> class ABGVOperator;

class H_E_dip;
class H_B_dip;
class H_M_pso;

class QMOperator;
class FockOperator;
class CoulombPotential;
class KineticOperator;
class NuclearPotential;
class ExchangePotential;
class XCPotential;
class XCFunctional;


class ResponseCalculation {
public:
    ResponseCalculation(RankOneTensorOperator<3> *h, double w, bool im, int d)
        : pert(h), freq(w), imag(im), dir(d) { }

    bool isDynamic() const { if (fabs(this->freq) > MachineZero) return true; return false; }
    bool isImaginary() const { return this->imag; }

    RankOneTensorOperator<3> *pert;
    double freq;
    bool imag;
    int dir;
};

class ResponseCalculations {
public:
    void push_back(RankOneTensorOperator<3> *h, double w, bool im, int d) {
        ResponseCalculation rsp_calc(h, w, im, d);
        bool unique = true;
        for (int i = 0; i < this->calculations.size(); i++) {
            const ResponseCalculation &i_calc = calculations[i];
            if ((i_calc.pert == rsp_calc.pert) and
                (fabs(i_calc.freq - rsp_calc.freq) < MachineZero) and
                (i_calc.dir == rsp_calc.dir)) unique = false;
        }
        if (unique) {
            this->calculations.push_back(rsp_calc);
        }
    }

    void clear() { this->calculations.clear(); }

    int size() const { return this->calculations.size(); }
    const ResponseCalculation &operator[](int i) const { return this->calculations[i]; }

protected:
    std::vector<ResponseCalculation> calculations;
};

class SCFDriver {
public:
    SCFDriver(Getkw &input);
    virtual ~SCFDriver() { }

    void setup();
    void run();
    void clear();

protected:
    // Top level input
    int max_scale;
    double rel_prec;
    double nuc_prec;

    // World input
    bool center_of_mass;
    std::vector<double> gauge;

    // Derivative operators
    string diff_kin;
    string diff_orb;
    string diff_pso;

    // Run parameters
    bool calc_scf_energy;
    bool calc_dipole_moment;
    bool calc_quadrupole_moment;
    bool calc_polarizability;
    bool calc_hyperpolarizability;
    bool calc_optical_rotation;
    bool calc_magnetizability;
    bool calc_nmr_shielding;
    bool calc_spin_spin_coupling;
    bool calc_hyperfine_coupling;
    ResponseCalculations rsp_calculations;

    bool pol_velocity;
    bool optrot_velocity;
    string nmr_perturbation;
    string optrot_perturbation;
    std::vector<double> pol_frequency;
    std::vector<double> optrot_frequency;
    std::vector<int> nmr_nucleus_k;
    std::vector<int> sscc_nucleus_k;
    std::vector<int> sscc_nucleus_l;
    std::vector<int> hfcc_nucleus_k;

    // Molecule input
    int mol_charge;
    int mol_multiplicity;
    std::vector<std::string> mol_coords;

    // Wavefunction input
    bool wf_restricted;
    std::string wf_method;

    // DFT input
    bool dft_spin;
    double dft_x_fac;
    double dft_cutoff;
    std::vector<double> dft_func_coefs;
    std::vector<std::string> dft_func_names;

    // Ground state SCF input
    std::string scf_start;
    int scf_kain;
    int scf_max_iter;
    int scf_rotation;
    bool scf_run;
    bool scf_canonical;
    bool scf_write_orbitals;
    double scf_orbital_thrs;
    double scf_property_thrs;
    double scf_lambda_thrs;
    std::vector<double> scf_orbital_prec;

    // Kinetic free SCF input
    int kin_free_max_iter;
    bool kin_free_run;
    bool kin_free_canonical;
    double kin_free_orb_thrs;
    double kin_free_prop_thrs;

    // Linear response input
    std::string rsp_start;
    int rsp_kain;
    int rsp_max_iter;
    bool rsp_run;
    bool rsp_canonical;
    bool rsp_write_orbitals;
    double rsp_orbital_thrs;
    double rsp_property_thrs;
    std::vector<int> rsp_directions;
    std::vector<double> rsp_orbital_prec;

    // File input
    std::string file_start_orbitals;
    std::string file_final_orbitals;
    std::string file_basis_set;
    std::string file_dens_mat;
    std::string file_fock_mat;
    std::string file_energy_vec;
    std::string file_mo_mat_a;
    std::string file_mo_mat_b;

    // Gauge origin
    double r_O[3];

    // SCF machinery
    HelmholtzOperatorSet *helmholtz;
    KAIN *kain;
    KAIN *kain_x;
    KAIN *kain_y;

    // MRA operators
    PoissonOperator *P;
    PHOperator<3> *PH_1;
    PHOperator<3> *PH_2;
    ABGVOperator<3> *ABGV_00;
    ABGVOperator<3> *ABGV_55;

    // Unperturbed quantities
    Molecule *molecule;
    Nuclei *nuclei;
    OrbitalVector *phi;
    KineticOperator *T;
    NuclearPotential *V;
    CoulombPotential *J;
    ExchangePotential *K;
    XCPotential *XC;
    FockOperator *fock;
    Eigen::MatrixXd F;

    OrbitalVector *phi_np1;
    CoulombPotential *J_np1;
    ExchangePotential *K_np1;
    XCPotential *XC_np1;
    FockOperator *fock_np1;
    Eigen::MatrixXd F_np1;

    // Perturbed quantities
    OrbitalVector *phi_x;
    OrbitalVector *phi_y;
    CoulombPotential *dJ;
    ExchangePotential *dK;
    XCPotential *dXC;
    FockOperator *d_fock;

    // XCFun
    XCFunctional *xcfun;

    // Perturbation operators
    H_E_dip  *h_E; // dH/dE
    H_B_dip  *h_B; // dH/dB
    H_M_pso **h_M; // dH/dM[K]

    bool sanityCheck() const;

    bool runGroundState();
    void runLinearResponse(const ResponseCalculation &rsp_calc);

    void calcGroundStateProperties();
    void calcLinearResponseProperties(const ResponseCalculation &rsp_calc);

    void setupInitialGroundState();
    void setupPerturbedOperators(const ResponseCalculation &rsp_calc);
    void setupPerturbedOrbitals(bool dynamic);

    void clearPerturbedOperators();
    void clearPerturbedOrbitals(bool dynamic);

    GroundStateSolver *setupInitialGuessSolver();
    OrbitalOptimizer *setupOrbitalOptimizer();
    EnergyOptimizer *setupEnergyOptimizer();
    LinearResponseSolver* setupLinearResponseSolver(bool dynamic);

    void setup_np1();
    void clear_np1();

    void printEigenvalues(OrbitalVector &orbs, Eigen::MatrixXd &f_mat);

    void extendRotationMatrix(const OrbitalVector &orbs, Eigen::MatrixXd &O);
};
