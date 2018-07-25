#pragma once

#include "qmfunction_utils.h"
#include "qmoperators.h"

#include "mrdft/XCFunctional.h"

class Getkw;

namespace mrchem {

class Nuclei;
class Molecule;
class OrbitalOptimizer;
class EnergyOptimizer;
class GroundStateSolver;
class LinearResponseSolver;
class HelmholtzVector;
class KAIN;

class H_E_dip;
class H_B_dip;
class H_M_pso;

class QMOperator;
class FockOperator;
class CoulombOperator;
class KineticOperator;
class NuclearOperator;
class ExchangeOperator;
class XCOperator;
class XCFunctional;
class ElectricFieldOperator;
class MagneticFieldOperator;

class ResponseCalculation  final {
public:
    ResponseCalculation(RankOneTensorOperator<3> *h, double w, bool im, int d)
        : pert(h), freq(w), imag(im), dir(d) { }

    bool isDynamic() const { if (std::abs(this->freq) > mrcpp::MachineZero) return true; return false; }
    bool isImaginary() const { return this->imag; }

    RankOneTensorOperator<3> *pert;
    double freq;
    bool imag;
    int dir;
};

class ResponseCalculations final {
public:
    void push_back(RankOneTensorOperator<3> *h, double w, bool im, int d) {
        ResponseCalculation rsp_calc(h, w, im, d);
        bool unique = true;
        for (int i = 0; i < this->calculations.size(); i++) {
            const ResponseCalculation &i_calc = calculations[i];
            if ((i_calc.pert == rsp_calc.pert) and
                (std::abs(i_calc.freq - rsp_calc.freq) < mrcpp::MachineZero) and
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

class SCFDriver final {
public:
    SCFDriver(Getkw &input);
    ~SCFDriver() { }

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
    bool center_of_charge;
    std::vector<double> gauge;

    // Derivative operators
    string diff_kin;
    string diff_orb;
    string diff_pso;
    string diff_dft;

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
    bool dft_use_gamma;
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

    // External field input
    bool ext_electric;
    bool ext_magnetic;
    Eigen::Vector3d ext_electric_field;
    Eigen::Vector3d ext_magnetic_field;
    
    // File input
    std::string file_start_orbitals;
    std::string file_final_orbitals;
    std::string file_basis_set;
    std::string file_dens_mat_a;
    std::string file_dens_mat_b;
    std::string file_fock_mat;
    std::string file_energy_vec;
    std::string file_mo_mat_a;
    std::string file_mo_mat_b;

    // Gauge origin
    double r_O[3];

    // SCF machinery
    HelmholtzVector *helmholtz;
    KAIN *kain;
    KAIN *kain_x;
    KAIN *kain_y;

    // MRA operators
    mrcpp::PoissonOperator *P;
    mrcpp::PHOperator<3> *PH_1;
    mrcpp::PHOperator<3> *PH_2;
    mrcpp::ABGVOperator<3> *ABGV_00;
    mrcpp::ABGVOperator<3> *ABGV_55;

    // Unperturbed quantities
    Molecule *molecule;
    Nuclei *nuclei;
    OrbitalVector *phi;
    KineticOperator *T;
    NuclearOperator *V;
    CoulombOperator *J;
    ExchangeOperator *K;
    XCOperator *XC;
    ElectricFieldOperator *Vext;
    FockOperator *fock;
    ComplexMatrix F;

    OrbitalVector *phi_np1;
    CoulombOperator *J_np1;
    ExchangeOperator *K_np1;
    XCOperator *XC_np1;
    FockOperator *fock_np1;
    ComplexMatrix F_np1;

    // Perturbed quantities
    OrbitalVector *phi_x;
    OrbitalVector *phi_y;
    CoulombOperator *dJ;
    ExchangeOperator *dK;
    XCOperator *dXC;
    FockOperator *d_fock;

    // XCFun
    mrdft::XCFunctional *xcfun;

    // Perturbation operators
    H_E_dip  *h_E; // dH/dE
    H_B_dip  *h_B; // dH/dB
    H_M_pso **h_M; // dH/dM[K]

    bool sanityCheck() const;

    bool runGroundState();
    void runLinearResponse(const ResponseCalculation &rsp_calc);

    void calcGroundStateProperties();
    void calcLinearResponseProperties(const ResponseCalculation &rsp_calc);

    mrdft::XCFunctional* setupFunctional(int order);
    void setupInitialGrid(mrdft::XCFunctional &func, const Molecule &mol);
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

    void printEigenvalues(OrbitalVector &orbs, ComplexMatrix &f_mat);

    void extendRotationMatrix(const OrbitalVector &orbs, ComplexMatrix &O);

    mrcpp::DerivativeOperator<3>* useDerivative(string derivative_name);
};

} //namespace mrchem
