/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include <Eigen/Eigenvalues>

#include "MRCPP/Printer"

#include "utils/NonlinearMaximizer.h"
#include "parallel.h"

using namespace std;

namespace mrchem {

/** General method to find the maximum of a nonlinear functional,
 *  given the Gradient and Hessian functions.
 *  The definitions of functional, make_gradient, make_hessian and do_step
 *  can be defined in a subclass.
 */
// We consider only the diagonal of the Hessian, until we are close to the minimum
int NonlinearMaximizer::maximize() {
    int i,j,k,l,iter,dcount;
    double mu, h2;
    double old_norm,new_norm,gradient_norm,value_functional,expected_change,relative_change;
    double value_functional_old,step_norm2,first_order,second_order;

    int print=0;//0: print nothing, 1: print only one line, 2: print one line per iteration; >50 print entire matrices
    int maxIter=150; //max number of iterations
    int CG_maxiter=5;// max number of iterations for the Conjugated Gradient solver for Newton step
    bool converged=false;
    double threshold=1.0e-12;//convergence when norm of gradient is smaller than threshold
    double CG_threshold=1.0e-12;//convergence when norm of residue is smaller than threshold
    double h=0.1; //initial value of trust radius, should be set small enough.
    bool wrongstep=false;
    int newton_step=0;
    double mu_min=1.0E-12;

    DoubleVector eiVal(N2h),old_step(N2h),diag(N2h),step(N2h) ;
    double mu_Newton_init=2.0,acc_fac=1.0;
    double mu_Newton=mu_Newton_init;
    double lamb1,sq,a,c,x,direction=0.0;
    int N_newton_step=0,newton_step_exact=0;
    double maxEiVal;

    //value_functional is what should be maximized (i.e. the sum of <i R i>^2 for orbitals)
    value_functional = this->functional();
    value_functional_old = value_functional;
    if(print>10 and mpi::orb_rank==0)cout <<"size "<<N2h<<endl;
    gradient_norm=this->make_gradient()/value_functional/N2h;// make the first gradient matrix

    if(print>100 and mpi::orb_rank==0)cout <<"gradient "<<gradient<<endl;

    //Start of iterations
    for(iter=1; iter<maxIter+1 && !converged ; iter++) {
        if(print>5 and mpi::orb_rank==0)cout <<" iteration  "<< iter << endl;
        dcount=0;
        for (i=0; i<N2h; i++) {
            eiVal(i)=this->get_hessian(i,i);
        }

        maxEiVal=eiVal(0);
        for (i=0; i<N2h; i++) {
            if(eiVal(i)>maxEiVal)maxEiVal=eiVal(i);
        }
        //    cout << " maxEiVal: "  << maxEiVal << endl;

        //We shift the eigenvalues, such that all are <0 (since we want a maximum)

        //To find mu, the trust radius is approximated using only the largest contribution in the series
        mu=mu_min;
        for (i=0; i<N2h; i++) {
            if(eiVal(i)+std::abs(this->gradient(i))/h>mu)mu=eiVal(i)+std::abs(this->gradient(i))/h+1.E-16;
        }

        diag= DoubleVector::Constant(N2h,-mu);
        diag += eiVal; //shifted eigenvalues of Hessian
        if(print>100 and mpi::orb_rank==0){
            cout << "mu and The shifted eigenvalues of H are: "  <<  mu<<  " h= "  << h <<  "  "  << diag << endl;
        }
        if(print>10 and mpi::orb_rank==0)cout << "largest eigenvalues of H: "  << maxEiVal << endl;

        for (i=0; i<N2h; i++) {
            step(i) = -this->gradient(i)/diag(i);
        }

        if(print>100 and mpi::orb_rank==0)cout << "step: "  << step << endl;
        direction=0.0;
        old_norm=0.0;
        new_norm=0.0;
        if(iter>1){
            for (i=0; i<N2h; i++) {
                direction+=old_step(i)*step(i);
                old_norm+=old_step(i)*old_step(i);
                new_norm+=step(i)*step(i);
            }
            direction=direction/std::sqrt(new_norm*old_norm);
        }
        if(direction>0.95){
            //direction is not changing much: accelerate!
            acc_fac=1.5;
            if(direction>0.98)acc_fac=2.0;
            if(direction>0.99)acc_fac=5.0;
            if(direction>0.995)acc_fac=10.0;
            if(print>10 and mpi::orb_rank==0)cout << " acceleration: "  << acc_fac<<endl;
            step*=acc_fac;
        }
        step_norm2 = step.transpose()*step ;
        first_order = this->gradient.transpose()*step ;
        second_order = 0.0;//step.transpose()*this->hessian*step;
        for (i=0; i<N2h; i++) {
            second_order += step(i)*step(i)*eiVal(i);
        }
        if(print>10 and mpi::orb_rank==0)cout << " gradient magnitude: "  << first_order*first_order/step_norm2<<endl;

        //Newton step when all eigenvalues are <0, and gradient/h sufficiently small
        newton_step=0;
        if(mu<mu_min*1.1){
            newton_step=1;
            N_newton_step++;
            if((N_newton_step>2&&(step_norm2<1.E-3||newton_step_exact==1))){
                if(print>10 and mpi::orb_rank==0)cout << "Taking Newton step  "<<  endl;
                newton_step_exact=1;
                mu_Newton*=0.2;//level shift for positive and zero eigenvalues
                DoubleVector  r(N2h), Ap(N2h), p(N2h), z(N2h),err(N2h), precond(N2h);
                double pAp, Hij;
                multiply_hessian(step, Ap);
                for (int i = 0 ; i<N2h ; i++){
                    //we take the shifted diagonal of the Hessian diag(i) as preconditioner
                    precond(i)=diag(i);
                    r(i) = -this->gradient(i)- Ap(i);
                    z(i) = r(i)/precond(i);
                    p(i) = z(i);
                }
                for (int iCG = 0; iCG <CG_maxiter; iCG++){
                    double rr =  r.transpose()*r;
                    double rz =  r.transpose()*z;
                    pAp = 0.0;
                    multiply_hessian(p, Ap);
                    for (int i = 0 ; i<N2h ; i++){
                        pAp += p(i)*Ap(i);
                    }
                    double a = rz/pAp;
                    for (int i = 0 ; i<N2h ; i++){
                        step(i) = step(i) + a*p(i);
                    }
                    multiply_hessian(step, Ap);
                    for (int i = 0 ; i<N2h ; i++){
                        //Note: Should we compute from step to avoid accumulation of errors or rather reuse Ap?
                        r(i) = -this->gradient(i) - Ap(i);
                        z(i) = r(i)/precond(i);
                    }
                    double b = r.transpose()*r;
                    if(std::sqrt(b)<CG_threshold) break;
                    b = r.transpose()*z;
                    b = b/rz;
                    for (int i = 0 ; i<N2h ; i++){
                        p(i) = z(i) +  b*p(i);
                    }
                    if(print>10) {
                        //test
                        double ee=0.0;
                        for (int i = 0 ; i<N2h ; i++){
                            double ss=0.0;
                            for (int j = 0 ; j<N2h ; j++){
                                ss += this->get_hessian(i,j) *step(j);
                            }
                            ee += (ss+ this->gradient(i))*(ss+ this->gradient(i));
                        }
                        if(mpi::orb_rank==0)cout<<iCG<<" error this iteration "<<std::sqrt(ee)<<endl;
                    }
                }
                step_norm2 = step.transpose()*step ;
                first_order = this->gradient.transpose()*step ;
                second_order = 0.0;
                for (int i = 0 ; i<N2h ; i++){
                    second_order += Ap(i)*step(i);
                }
            }else{
                newton_step_exact=0;
                mu_Newton=mu_Newton_init;
            }

            if(print==20 and mpi::orb_rank==0)cout<<endl<<" 2nd  "<<second_order<<" fi*fi/d  "<<x<<" s g "<<step.transpose()*gradient<<endl;
        }else{
            N_newton_step=0;
            newton_step_exact=0;
            mu_Newton=mu_Newton_init;
        }

        //    if(print==2){cout<<setw(22)<< step_norm2;}
        if(print>100 and mpi::orb_rank==0){
            cout << "step : "  ;
            for (i=0; i<N2h; i++) {
                cout << " "<< step(i) ;
            }
            cout << endl ;
        }

        expected_change=first_order +0.5*second_order;

        if(print>10 and mpi::orb_rank==0)  cout << "step size  "<< std::sqrt(step_norm2)   <<endl;
        if(print>10 and mpi::orb_rank==0)cout << "expected first, second order and total change in r*r  ";
        if(print>10 and mpi::orb_rank==0)cout << first_order<< " " << 0.5*second_order <<" " <<expected_change  <<endl;

        this->do_step(step);
        value_functional = this->functional();

        if(print>10 and mpi::orb_rank==0) cout <<" r*r  "<< value_functional << " change in r*r  "<< value_functional-value_functional_old << endl;

        //relative_change is the size of  second order change compared to actual change
        // = 0 if no higher order contributions
        relative_change = std::abs(expected_change-(value_functional-value_functional_old))
                /(1.0E-25+std::abs(value_functional-value_functional_old ));
        gradient_norm=this->make_gradient()/value_functional/N2h; //update gradient
        direction=0.0;
        old_norm=0.0;
        new_norm=0.0;
        if(iter>1){
            for (i=0; i<N2h; i++) {
                direction+=old_step(i)*step(i);
                old_norm+=old_step(i)*old_step(i);
                new_norm+=step(i)*step(i);
            }
            direction=direction/std::sqrt(new_norm*old_norm);
        }
        //    if(print==2)cout<<setw(22) <<gradient_norm<< setw(22) <<direction;
        if(print==2)cout<< setprecision(3)<<setw(10)<<gradient_norm;
        if(print==2)cout<<setprecision(4)  <<setw(9) <<direction;

        if(value_functional-value_functional_old <-threshold){
            //jumped to far, value_functional is decreasing
            //jump back!
            h=0.05*h;
            if(print>=2)cout << endl;
            if(print>=2)cout <<"Wrong step! Returning ";
            wrongstep=true;
            N_newton_step=0;
            //go back
            this->do_step(-step);
            gradient_norm=this->make_gradient()/value_functional_old/N2h; //update gradient
        }else{
            //take this point as new reference
            value_functional_old = value_functional;
            wrongstep=false;
            old_step=step;
        }

        if(newton_step!=1){
            //change the trust radius h, based on how close the resulting r*r is to the
            //1st+2nd order predicted value
            if(relative_change<0.0001){
                h=10.0*h;
            }else if(relative_change<0.001) {
                h=5.0*h;
            }else if(relative_change<0.01) {
                h=2.0*h;
            }else if(relative_change<0.1) {
                h=1.2*h;
            }

            if(relative_change>2.0) {
                h=0.2*h;
            }else if(relative_change>0.5) {
                h=0.5*h;
            }else if(relative_change>0.2) {
                h=0.8*h;
            }
            if(h<1.E-8)h=1.E-8;
        }
        if(print>10 and mpi::orb_rank==0)cout <<"trust radius set  to "<< h << " test: "<<relative_change << " mu: "<<mu<<  " maxeival: "<<maxEiVal <<endl;
        if(print>10 and mpi::orb_rank==0)cout << "gradient norm " << gradient_norm<< endl;

        if(gradient_norm<threshold&&maxEiVal<10*std::sqrt(std::abs(threshold))){
            //finished
            converged=true;
        }
        if(print==2&&!wrongstep){cout<<setw(10)<< dcount;}
        if(print==2)cout<<endl;

    }//iterations

    if(print>15 and mpi::orb_rank==0){
        cout<<"Exact Hessian (upper left corner)"<<endl;
        for (int i = 0 ; i<20 and i<N2h ; i++){
            for (int j = 0 ; j<20 and j<N2h ; j++){
                printf (" %8.5f",this->get_hessian(i,j));
            }
            printout(0, std::endl);
        }
    }

    if (iter < maxIter) {
        iter = iter - 1;
    } else {
        iter = -1;
    }
    return iter;
}


} //namespace mrchem
