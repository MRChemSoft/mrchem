#include "DensityProjector.h"
#include "OrbitalVector.h"
#include "Density.h"
#include "SerialFunctionTree.h"

extern MultiResolutionAnalysis<3> *MRA;

using namespace std;

DensityProjector::DensityProjector(double prec, int max_scale)
    : add(prec, max_scale),
      mult(prec, max_scale),
      grid(max_scale) {
      }

void DensityProjector::setPrecision(double prec) {
    this->add.setPrecision(prec);
    this->mult.setPrecision(prec);
}

void DensityProjector::operator()(Density &rho, Orbital &phi) {
    if (rho.hasTotal()) MSG_ERROR("Density not empty");
    if (rho.hasSpin()) MSG_ERROR("Density not empty");
    if (rho.hasAlpha()) MSG_ERROR("Density not empty");
    if (rho.hasBeta()) MSG_ERROR("Density not empty");

    double occ = 1.0;
    if (not rho.isSpinDensity()) occ = (double) phi.getOccupancy();

    double prec = mult.getPrecision();
    if (prec < 0.0) MSG_ERROR("Adaptive multiplication with negative prec");

    FunctionTreeVector<3> sum_vec;
    if (phi.hasReal()) {
        FunctionTree<3> *real_2 = new FunctionTree<3>(*MRA);
        this->grid(*real_2, phi.real());
        this->mult(*real_2, occ, phi.real(), phi.real(), 1);
        sum_vec.push_back(real_2);
    }
    if (phi.hasImag()) {
        FunctionTree<3> *imag_2 = new FunctionTree<3>(*MRA);
        this->grid(*imag_2, phi.imag());
        this->mult(*imag_2, occ, phi.imag(), phi.imag(), 1);
        sum_vec.push_back(imag_2);
    }

    if (rho.isSpinDensity()) {
        if (phi.getSpin() == Orbital::Paired) {
            rho.allocAlpha();
            rho.allocBeta();
            this->grid(rho.alpha(), sum_vec);
            this->grid(rho.beta(), sum_vec);
            this->add(rho.alpha(), sum_vec, 0);
            this->add(rho.beta(), sum_vec, 0);
        }
        if (phi.getSpin() == Orbital::Alpha) {
            rho.allocAlpha();
            this->grid(rho.alpha(), sum_vec);
            this->add(rho.alpha(), sum_vec, 0);

            rho.allocBeta();
            this->grid(rho.beta(), sum_vec);
            rho.beta().setZero();
        }
        if (phi.getSpin() == Orbital::Beta) {
            rho.allocBeta();
            this->grid(rho.beta(), sum_vec);
            this->add(rho.beta(), sum_vec, 0);

            rho.allocAlpha();
            this->grid(rho.alpha(), sum_vec);
            rho.alpha().setZero();
        }
        FunctionTreeVector<3> tot_vec;
        tot_vec.push_back(1.0, &rho.alpha());
        tot_vec.push_back(1.0, &rho.beta());
        rho.allocTotal();
        this->grid(rho.total(), tot_vec);
        this->add(rho.total(), tot_vec, 0);

        FunctionTreeVector<3> spin_vec;
        spin_vec.push_back( 1.0, &rho.alpha());
        spin_vec.push_back(-1.0, &rho.beta());
        rho.allocSpin();
        this->grid(rho.spin(), spin_vec);
        this->add(rho.spin(), spin_vec, 0);
    } else {
        rho.allocTotal();
        this->grid(rho.total(), sum_vec);
        this->add(rho.total(), sum_vec, 0);
    }
    sum_vec.clear(true);
}

void DensityProjector::operator()(Density &rho, OrbitalVector &phi) {
    if (rho.hasTotal()) MSG_ERROR("Density not empty");
    if (rho.hasAlpha()) MSG_ERROR("Density not empty");
    if (rho.hasBeta()) MSG_ERROR("Density not empty");

    int nOrbs = phi.size();
    double prec = mult.getPrecision();
    if (prec < 0.0) MSG_ERROR("Adaptive addition with negative prec");
    add.setPrecision(prec/nOrbs);

    FunctionTreeVector<3> total_vec, spin_vec, alpha_vec, beta_vec;
    vector<Density *> dens_vec;
    vector<int> rho_i_Ix;
    Density* rho_tmp1 = 0;
    Density* rho_tmp2 = 0;
    Density* rho_tmp = 0;
    Density* rho_i = 0;
    if (mpiOrbSize>1) {
	//only master does the summation
	for (int i_Orb = 0; i_Orb < phi.size(); i_Orb++) {
	    rho_i = new Density(rho);	
	    if (i_Orb%mpiOrbSize == mpiOrbRank) {
		Orbital &phi_i = phi.getOrbital(i_Orb);
		//cout<<mpiOrbRank<<" size orbital "<< i_Orb<< " in MB "<<int((phi_i.real().getNNodes()* phi_i.real().getKp1_d()*8.0*8.0)/1024/1024)<<endl;
		(*this)(*rho_i, phi_i);
		if (mpiOrbRank != 0) rho_i->send_Density(0, i_Orb);
	    }
	    //add on the fly
	    if (mpiOrbRank == 0 and i_Orb == 0){
		//first iteration does not sum only receive in tmp1
		if (i_Orb%mpiOrbSize != mpiOrbRank) {
		    rho_tmp1->Rcv_Density(i_Orb%mpiOrbSize, i_Orb);
		}else{
		    rho_tmp1 = rho_i;
		    rho_i = rho_tmp2;
		}
	    }else if (mpiOrbRank == 0) {
		if (i_Orb%mpiOrbSize != mpiOrbRank) rho_i->Rcv_Density(i_Orb%mpiOrbSize, i_Orb);
		//exchange pointers. old result is in tmp1: tmp1 = rho_i+tmp2
		rho_tmp2 = rho_tmp1;
		rho_tmp1 = new Density(rho);
		if (i_Orb == phi.size()-1) {
		    //last iteration, put result into rho
		    rho_tmp1=&rho;
		}		
		if (rho_i->hasTotal()) {
		    if(not rho_tmp1->hasTotal())rho_tmp1->allocTotal();
		    total_vec.push_back(&rho_tmp2->total());
		    total_vec.push_back(&rho_i->total());
		    this->grid(rho_tmp1->total(), total_vec);// kopierer grid fra funksjonene i total_vec
		    this->add(rho_tmp1->total(), total_vec,0);
		    total_vec.clear(true);
		}
		if (rho_i->hasSpin()) {
		    if(not rho_tmp1->hasSpin())rho_tmp1->allocSpin();
		    spin_vec.push_back(&rho_tmp2->spin());
		    spin_vec.push_back(&rho_i->spin());
		    this->grid(rho_tmp1->spin(), spin_vec);
		    this->add(rho_tmp1->spin(), spin_vec,0);
		    spin_vec.clear(true);
		}
		if (rho_i->hasAlpha()) {
		    if(not rho_tmp1->hasAlpha())rho_tmp1->allocAlpha();
		    alpha_vec.push_back(&rho_tmp2->alpha());
		    alpha_vec.push_back(&rho_i->alpha());
		    this->grid(rho_tmp1->alpha(), alpha_vec);
		    this->add(rho_tmp1->alpha(), alpha_vec,0);
		    alpha_vec.clear(true);
		}
		if (rho_i->hasBeta()) {
		    if(not rho_tmp1->hasBeta())rho_tmp1->allocBeta();
		    beta_vec.push_back(&rho_tmp2->beta());
		    beta_vec.push_back(&rho_i->beta());
		    this->grid(rho_tmp1->beta(), beta_vec);
		    this->add(rho_tmp1->beta(), beta_vec,0);
		    beta_vec.clear(true);
		}
	    }
	}
    }else{      
	//Serial processing
	for (int i = 0; i < phi.size(); i++) {
	    Orbital &phi_i = phi.getOrbital(i);
	    Density *rho_i = new Density(rho);
	    (*this)(*rho_i, phi_i);
	    dens_vec.push_back(rho_i);
	    if (rho_i->hasTotal()) total_vec.push_back(&rho_i->total());
	    if (rho_i->hasSpin()) spin_vec.push_back(&rho_i->spin());
	    if (rho_i->hasAlpha()) alpha_vec.push_back(&rho_i->alpha());
	    if (rho_i->hasBeta()) beta_vec.push_back(&rho_i->beta());
	}
    } 

    if (mpiOrbRank == 0 and not rho.isShared()) {
	if (total_vec.size() > 5) {
	    rho.allocTotal();
	    this->add(rho.total(), total_vec);
	} else if (total_vec.size() > 0) {
	    rho.allocTotal();
	    this->grid(rho.total(), total_vec);
	    this->add(rho.total(), total_vec, 0);
	}
	if (spin_vec.size() > 5) {
	    rho.allocSpin();
	    this->add(rho.spin(), spin_vec);
	} else if (spin_vec.size() > 0) {
	    rho.allocSpin();
	    this->grid(rho.spin(), spin_vec);
	    this->add(rho.spin(), spin_vec, 0);
	}
	if (alpha_vec.size() > 5) {
	    rho.allocAlpha();
	    this->add(rho.alpha(), alpha_vec);
	} else if (alpha_vec.size() > 0) {
	    rho.allocAlpha();
	    this->grid(rho.alpha(), alpha_vec);
	    this->add(rho.alpha(), alpha_vec, 0);
	}
	if (beta_vec.size() > 5) {
	    rho.allocBeta();
	    this->add(rho.beta(), beta_vec);
	} else if (beta_vec.size() > 0) {
	    rho.allocBeta();
	    this->grid(rho.beta(), beta_vec);
	    this->add(rho.beta(), beta_vec, 0);
	}

	for (int i = 0; i < dens_vec.size(); i++) {
	    dens_vec[i]->clear();
	    delete dens_vec[i];
	    dens_vec[i] = 0;
	}
    }


    //if(mpiOrbRank==0)cout<<"size density MB "<<int((rho.getNNodes()* rho.total().getKp1_d()*8.0*8.0)/1024/1024)<<endl;
    if (mpiOrbSize > 1) {
	//we always broadcast density
	//If the density is shared, only metadata will be sent/received
	if (mpiOrbRank == 0) {
	    for (int i_MPI = 1; i_MPI < mpiOrbSize; i_MPI++) {
		rho.send_Density(i_MPI, 54);
	    }
	}else{
	    //do nothing, only receive Density
	    rho.Rcv_Density(0, 54);
	}
    }
}
