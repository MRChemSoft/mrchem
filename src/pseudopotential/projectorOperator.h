#pragma once

#include "tensor/RankZeroOperator.h"
#include "mrchem.h"
#include "pseudopotential/pseudopotential.h"
#include "pseudopotential/projector.h"
#include "qmfunctions/Orbital.h"
#include "chemistry/Molecule.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMOperator.h"
#include <string>

#include "MRCPP/Printer"

class magneticQuantumNumberProjector {
    public:
    std::vector<ProjectorFunction> iProj;
    int nProj;
};

class angularMomentumProjector {
    public:
    std::vector<magneticQuantumNumberProjector> mProj;
    int nM;
};

class AtomProjector {
    public:
    std::vector<angularMomentumProjector> lProj;
    int numberOfAngMom;
};

class ProjectorOperatorQM final : public mrchem::QMOperator {

    std::vector<PseudopotentialData> pp;
    std::vector<AtomProjector> proj;
    double prec;

public:
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

    void setup(double prec) {
        this->prec = prec;
    }

    void clear() {
    }

protected:

mrchem::Orbital apply(mrchem::Orbital phi) {
    // std::cout << "Applying projector operator" << std::endl;
    // loop over all atoms
    ComplexDouble dotComplex;

    std::vector<ComplexDouble> complexCoefficients;
    mrchem::ComplexFunctionVector complexFunctionVector;

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
                    dotComplex = mrcpp::cplxfunc::dot(phi, *proj[iat].lProj[l].mProj[mm].iProj[ip].projector_ptr);
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

mrchem::Orbital dagger(mrchem::Orbital phi) {
    return apply(phi);
}

mrchem::ComplexDouble evalf(const mrcpp::Coord<3> &r) const {
    return ComplexDouble(0.0, 0.0);
}

mrchem::QMOperatorVector apply(std::shared_ptr<mrchem::QMOperator> &O) {
    NOT_IMPLEMENTED_ABORT;
}

};

class ProjectorOperator : public mrchem::RankZeroOperator {

public:
    ProjectorOperator(mrchem::Nuclei const &nucs, double prec) {
        auto qmOperator = std::make_shared<ProjectorOperatorQM>(nucs, prec);
        mrchem::RankZeroOperator &pp = (*this);
        pp = qmOperator;
    }


};