#pragma once

/**
 * @file projectorOperator.h
 * @brief Non-local pseudopotential projector operator classes.
 *
 * This file contains the classes needed to apply the non-local part of
 * pseudopotentials using projector functions. The projectors are organized
 * hierarchically by atom, angular momentum (l), and magnetic quantum number (m).
 */

#include "tensor/RankZeroOperator.h"
#include "mrchem.h"
#include "pseudopotential/pseudopotential.h"
#include "pseudopotential/projector.h"
#include "qmfunctions/Orbital.h"
#include "chemistry/Molecule.h"
#include "qmoperators/QMOperator.h"
#include <string>

#include "MRCPP/Printer"

/**
 * @class magneticQuantumNumberProjector
 * @brief Container for projector functions with the same magnetic quantum number m.
 *
 * For a given angular momentum l and magnetic quantum number m, this class
 * holds all projector functions (indexed by i) that share these quantum numbers.
 */
class magneticQuantumNumberProjector {
public:
    std::vector<ProjectorFunction> iProj; ///< Projector functions for this (l,m) pair
    int nProj;                             ///< Number of projectors
};

/**
 * @class angularMomentumProjector
 * @brief Container for projectors with the same angular momentum l.
 *
 * For a given angular momentum l, this class holds projectors for all
 * magnetic quantum numbers m ranging from -l to +l.
 */
class angularMomentumProjector {
public:
    std::vector<magneticQuantumNumberProjector> mProj; ///< Projectors indexed by m
    int nM;                                             ///< Number of m values (2l+1)
};

/**
 * @class AtomProjector
 * @brief Container for all projectors associated with a single atom.
 *
 * This class organizes the non-local pseudopotential projectors for one atom,
 * with projectors grouped by angular momentum l.
 */
class AtomProjector {
public:
    std::vector<angularMomentumProjector> lProj; ///< Projectors indexed by angular momentum l
    int numberOfAngMom;                           ///< Number of angular momentum channels
};

/**
 * @class ProjectorOperatorQM
 * @brief Quantum mechanical operator for non-local pseudopotential projectors.
 *
 * This class implements the non-local part of the pseudopotential using
 * the Goedecker-Teter-Hutter (GTH) form. It applies projector operators to orbitals
 * by computing inner products with projector functions and forming linear
 * combinations weighted by the pseudopotential matrix elements h_ij^l.
 */
class ProjectorOperatorQM final : public mrchem::QMOperator {

    std::vector<PseudopotentialData> pp; ///< Pseudopotential data for each atom
    std::vector<AtomProjector> proj;      ///< Projector functions for each atom
    double prec;                          ///< Numerical precision

public:
    /**
     * @brief Constructs the projector operator for a set of nuclei.
     * @param nucs The nuclei with associated pseudopotential data.
     * @param prec The numerical precision for projector functions.
     */
    ProjectorOperatorQM(mrchem::Nuclei const &nucs, double prec){

        // mrchem::Nuclei nucs = molecule.getNuclei();
        // for (int i = 0; i < nucs.size(); i++){
        //     std::string elem = nucs[i].getElement().getSymbol();
        //     std::string fname = "psppar." + elem;
        //     this->pp.push_back(PseudopotentialData(fname));
        //     std::cout << "Pseudopotential data for atom " << i << " loaded" << std::endl;
        //     pp[i].print();
        // }

        // std::cout << "Pseudopotential data loaded" << std::endl;
        for (int i = 0; i < nucs.size(); i++){
            if (!nucs[i].hasPseudopotential()){
                MSG_ABORT("No pseudopotential for atom " + std::to_string(i));
            }
            this->pp.push_back(*nucs[i].getPseudopotentialData());
            // nucs[i].getPseudopotentialData()->print();
        }

        this->pp = pp;
        this->prec = prec;
        int npp = 0;

        // loop over all atoms and create projectors
        for (int i = 0; i < nucs.size(); i++) {
            // std::cout << "Creating projectors for atom " << i << std::endl << std::endl;
            mrcpp::Coord<3> pos = nucs[i].getCoord();
            proj.push_back(AtomProjector());
            for (int l = 0; l < pp[i].nsep; l++) {
                // std::cout << "Creating angular momentum projectors for momentum " << l << std::endl;
                proj[i].lProj.push_back(angularMomentumProjector());
                for (int m = -l; m <= l; m++) {
                    // std::cout << "Creating magnetic quantum number projectors for magnetic quantum number " << m << std::endl;
                    int mIndex = m + l;
                    proj[i].lProj[l].mProj.push_back(magneticQuantumNumberProjector());
                    for (int idim = 0; idim < pp[i].dim_h[l]; idim++){
                        // proj.push_back(ProjectorFunction(pos, pp[i].rl[l], isep, l, m, prec));
                        // std::cout << "Creating ProjectorFunction " << l << " " << m << " " << idim << std::endl;
                        ProjectorFunction pppp(pos, pp[i].rl[l], idim, l, m, prec);
                        proj[i].lProj[l].mProj[mIndex].iProj.push_back(pppp);
                        // std::cout << "ProjectorFunction constructed " << i << std::endl;
                        // std::cout << "i = " << i << std::endl;
                        // std::cout << "nsep = " << pp[i].nsep << std::endl;
                        proj[i].numberOfAngMom = pp[i].nsep;
                        // std::cout << "ProjectorFunction added to projector" << std::endl;
                        proj[i].lProj[l].nM = 2*pp[i].nsep + 1;
                        // std::cout << "ProjectorFunction added to projector" << std::endl;
                        proj[i].lProj[l].mProj[mIndex].nProj = pp[i].dim_h[l];
                        // std::cout << "End of loopsss" << std::endl << std::endl;
                        npp++;
                    }
                }
            }
        }
        // std::cout << "ProjectorOperator constructed                      aasdfa asdf" << std::endl;
    }

    /**
     * @brief Sets up the operator with the given precision.
     * @param prec The numerical precision.
     */
    void setup(double prec) {
        this->prec = prec;
    }

    /**
     * @brief Clears the operator state.
     */
    void clear() {
    }

protected:

    /**
     * @brief Applies the non-local projector operator to an orbital.
     * @param phi The input orbital.
     * @return The resulting orbital after applying the projector operator.
     *
     * The application follows the GTH separable form: sum over atoms, angular
     * momenta, and projector indices of h_ij^l * <p_i|phi> * |p_j>.
     */
    mrchem::Orbital apply(mrchem::Orbital phi) {
    // std::cout << "Applying projector operator" << std::endl;
    // loop over all atoms
    ComplexDouble dotComplex;

    std::vector<ComplexDouble> complexCoefficients;
    // Note that mrcpp::CompFunctionVector is more than just a vector of CompFunction
    std::vector<mrcpp::CompFunction<3>> complexFunctionVector;

    for (int iat = 0; iat < proj.size(); iat++) {
        // loop over all angular momenta
        for (int l = 0; l < pp[iat].nsep; l++){
            // loop over all magnetic quantum numbers
            // std::cout << "h: " << pp[iat].h[l] << std::endl;
            for (int m = -l; m <= l; m++){
                int mm = m + l;
                // loop over all projectors
                Eigen::VectorXd dot_products(pp[iat].dim_h[l]);
                // std::cout << "Projector " << iat << " " << l << " " << m << std::endl;
                // std::cout << "Number of projectors " << pp[iat].dim_h[l] << std::endl;
                for (int ip = 0; ip < pp[iat].dim_h[l]; ip++){
                    // dotComplex = mrchem::qmfunction::dot(phi, proj[iat].lProj[l].mProj[m].iProj[ip]);
                    // std::cout << "computing dot product " << ip << std::endl;
                    mrcpp::Coord<3> r = {0.0, 0.0, 0.3};
                    // std::cout << "projector value at origin: " << proj[iat].lProj[l].mProj[mm].iProj[ip].real().evalf(r) << std::endl;
                    dotComplex = mrcpp::dot(phi, *proj[iat].lProj[l].mProj[mm].iProj[ip].projector_ptr);
                    dot_products(ip) = dotComplex.real();
                    // std::cout << "Dot product " << ip << " " << dotComplex << std::endl;
                }
                dot_products = pp[iat].h[l] * dot_products;
                // loop over all projectors
                for (int ip = 0; ip < pp[iat].dim_h[l]; ip++){
                    complexCoefficients.push_back(dot_products(ip));
                    complexFunctionVector.push_back(*proj[iat].lProj[l].mProj[mm].iProj[ip].projector_ptr);
                }
            }
        }

    }
    // convert complexCoefficients to Eigen Vector:
    mrchem::ComplexVector complexCoefficientsEigen = Eigen::Map<Eigen::VectorXcd>(complexCoefficients.data(), complexCoefficients.size());

    mrchem::Orbital result;
    // result.add()
    // mrchem::qmfunction::linear_combination(result, complexCoefficientsEigen, complexFunctionVector, prec);

    // std::cout << "size of complexCoefficients " << complexCoefficients.size() << std::endl;

    for (int i = 0; i < complexCoefficients.size(); i++){
        // std::cout << "Adding to result " << i << " " << complexCoefficients[i] << std::endl;
        result.add(complexCoefficients[i], complexFunctionVector[i]);
    }

    return result;
}

    /**
     * @brief Applies the adjoint (dagger) operator to an orbital.
     * @param phi The input orbital.
     * @return The resulting orbital (same as apply since the operator is Hermitian).
     */
    mrchem::Orbital dagger(mrchem::Orbital phi) {
        return apply(phi);
    }

    /**
     * @brief Evaluates the operator at a given coordinate (not implemented).
     * @param r The coordinate.
     * @return Zero (placeholder).
     */
    mrchem::ComplexDouble evalf(const mrcpp::Coord<3> &r) const {
        return ComplexDouble(0.0, 0.0);
    }

    /**
     * @brief Applies the operator to another operator (not implemented).
     * @param O The input operator.
     */
    mrchem::QMOperatorVector apply(std::shared_ptr<mrchem::QMOperator> &O) {
        NOT_IMPLEMENTED_ABORT;
    }
};

/**
 * @class ProjectorOperator
 * @brief High-level interface for the non-local pseudopotential projector operator.
 *
 * This class wraps ProjectorOperatorQM as a RankZeroOperator, providing an
 * interface for use in the MRChem framework.
 */
class ProjectorOperator : public mrchem::RankZeroOperator {

public:
    /**
     * @brief Constructs a projector operator for a set of nuclei.
     * @param nucs The nuclei with associated pseudopotential data.
     * @param prec The numerical precision.
     */
    ProjectorOperator(mrchem::Nuclei const &nucs, double prec) {
        auto qmOperator = std::make_shared<ProjectorOperatorQM>(nucs, prec);
        mrchem::RankZeroOperator &pp = (*this);
        pp = qmOperator;
    }
};
