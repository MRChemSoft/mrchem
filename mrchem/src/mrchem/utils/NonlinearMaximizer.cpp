#include <Eigen/Eigenvalues>

#include "MRCPP/Printer"

#include "utils/NonlinearMaximizer.h"

using namespace std;

namespace mrchem {

/** General method to find the maximum of a nonlinear functional,
 *  given the Gradient and Hessian functions.
 *  The definitions of functional, make_gradient, make_hessian and do_step
 *  can be defined in a subclass.
 */
int NonlinearMaximizer::maximize() {
    int i,j,k,l,iter,dcount;
    double mu, h2;
    double old_norm,new_norm,gradient_norm,value_functional,expected_change,relative_change;
    double value_functional_old,step_norm2,first_order,second_order;

    int print=0;//0: print nothing, 1: print only one line, 2: print one line per iteration; >50 print entire matrices
    int maxIter=150; //max number of iterations
    bool converged=false;
    double threshold=1.0e-12;//convergence when norm of gradient is small than threshold
    double h=0.1; //initial value of trust radius, should be set small enough.
    bool wrongstep=false;
    int newton_step=0;
    double mu_min=1.0E-12;

    DoubleMatrix eiVec(N2h,N2h),Hess_tmp(N2h,N2h);
    DoubleVector fi(N2h),eiVal(N2h),old_step(N2h),gradient_tmp(N2h),diag(N2h),sigma(N2h),step(N2h) ;
    double mu_Newton_init=2.0,acc_fac=1.0;
    double mu_Newton=mu_Newton_init;
    double lamb1,sq,a,c,x,U00,U01,fac=10.0,direction=0.0;
    int N_newton_step=0,newton_step_exact=0;
    double maxEiVal;
    //value_functional is what should be maximized (i.e. the sum of <i R i>^2 for orbitals)
    value_functional = this->functional();
    value_functional_old = value_functional;
    if(print>10)cout <<"size "<<N2h<<endl;
    gradient_norm=this->make_gradient()/value_functional/N2h;// make the first gradient matrix

    if(print>100)cout <<"gradient "<<gradient<<endl;

    if(print==2){cout <<"iteration "<<" step_type "<<"     step  "<<"     r*r     "<<
                        "  2nd_order"<<"  actual_diff "<<"gradient "<< "dir_change" << " nondiag_used"<< endl;}

    //Start of iteration
    for(iter=1; iter<maxIter+1 && !converged ; iter++) {
        if(print>5)cout <<" iteration  "<< iter << endl;
        if(print==2){
            cout<<setw(6)<< iter;
        }
        if(!wrongstep){
            this->make_hessian();//update Hessian
            Hess_tmp=this->hessian;
        }
        dcount=0;
        if(N2h>100){
            //large system, use approximate diagonalizatio

            for (i=0; i<N2h; i++) {
                for (j=0; j<N2h; j++) {
                    eiVec(i,j)=0.0;
                }
                eiVec(i,i)=1.0;
            }
            //Roughly diagonalize Hessian. Diagnoalize all off-diagonal elements independently, but succesively.
            // can be iterated (happens only after Newton steps): iter_diag is the number of these iterations
            //fac is the threshold for elements to treat (most of them are very small or zero).
            int iter_diag=1;
            double iter_fac=2.0;
            int Newton_iter=10000;//for testing only
            iter_fac=fac;
            iter_diag=1;
            if(N_newton_step>0){
                iter_fac=fac/0.4;
                iter_diag=2;
            }
            if(newton_step_exact==1)iter_fac=10.0;//we don't need it anyway
            if(iter<Newton_iter){
                for (l=0; l<iter_diag  ; l++) {
                    for (i=0; i<N2h; i++) {
                        for (int j=i+1; j<N2h; j++){
                            a=Hess_tmp(i,i)-Hess_tmp(j,j);
                            c=Hess_tmp(j,i);
                            if(std::abs(c)>iter_fac*std::abs(a)+1.0E-10){
                                dcount++;
                                sq=std::sqrt(a*a+4*c*c);
                                if(a>0){
                                    lamb1=0.5*(sq-a);
                                    //	      lamb2=lamb1-sq;
                                }else{
                                    //	      lamb2=0.5*(sq-a);
                                    lamb1=-0.5*(sq+a);
                                }
                                x=lamb1/c;
                                U00=1.0/std::sqrt(x*x+1.0);
                                U01=-U00*x;

                                //update Eigenvectors eiVec=eiVec*U
                                for (k=0; k<N2h; k++) {
                                    x=eiVec(k,i)*U00-eiVec(k,j)*U01;
                                    eiVec(k,j)=eiVec(k,i)*U01+eiVec(k,j)*U00;
                                    eiVec(k,i)=x;
                                }

                                //update Hessian Hess_tmp=eiVec.transpose()*Hess*eiVec;
                                for (k=0; k<N2h; k++) {
                                    x=Hess_tmp(k,i)*U00-Hess_tmp(k,j)*U01;
                                    Hess_tmp(k,j)=Hess_tmp(k,i)*U01+Hess_tmp(k,j)*U00;
                                    Hess_tmp(k,i)=x;
                                }
                                //could be neglected! :
                                for (k=0; k<N2h; k++) {
                                    x=Hess_tmp(i,k)*U00-Hess_tmp(j,k)*U01;//U is transposed
                                    Hess_tmp(j,k)=Hess_tmp(i,k)*U01+Hess_tmp(j,k)*U00;//U is transposed
                                    Hess_tmp(i,k)=x;
                                }
                            }
                        }
                    }
                    for (i=0; i<N2h; i++) {
                        for (int j=i+1; j<N2h; j++){
                            Hess_tmp(i,j)=Hess_tmp(j,i);//ensure symmetri
                        }
                    }
                    iter_fac*=0.4;
                }
            }

            for (i=0; i<N2h; i++) {
                eiVal(i)=Hess_tmp(i,i);
            }

        }else{
            //use exact diagonalization
            Eigen::SelfAdjointEigenSolver<DoubleMatrix> eigensolver(this->hessian);
            eiVec=eigensolver.eigenvectors();
            eiVal=eigensolver.eigenvalues();
            //Hessian=eiVec* eiVal *eiVec.transpose()
        }
        if(print>100)cout <<"Hess "<<hessian<<endl;
        if(print>100)cout <<"Approximate EigenVec "<<eiVec<<endl;

        fi=eiVec.transpose()*this->gradient; //gradient in eigenvector basis
        if(print>100){
            cout << "fi: "  << fi << endl;
        }

        maxEiVal=eiVal(0);
        for (i=0; i<N2h; i++) {
            if(eiVal(i)>maxEiVal)maxEiVal=eiVal(i);
        }
        //    cout << " maxEiVal: "  << maxEiVal << endl;

        //We shift the eigenvalues, such that all are <0 (since we want a maximum)
        //eiVal(N2h-1) is highest eigenvalue

        //To find mu, the trust radius is approximated using only the largest contribution in the series
        mu=mu_min;
        for (i=0; i<N2h; i++) {
            if(eiVal(i)+std::abs(fi(i))/h>mu)mu=eiVal(i)+std::abs(fi(i))/h+1.E-16;
        }

        diag= DoubleVector::Constant(N2h,-mu);
        diag += eiVal; //shifted eigenvalues of Hessian
        if(print>100){
            cout << "mu and The shifted eigenvalues of H are: "  <<  mu<<  " h= "  << h <<  "  "  << diag << endl;
        }
        if(print>10)cout << "largest eigenvalues of H: "  << maxEiVal << endl;

        h2=0.0;
        for (i=0; i<N2h; i++) {
            sigma(i) = fi(i)/diag(i);
            h2+=sigma(i)*sigma(i);
        }
        if(print>100)cout << "sigma: "  << sigma << endl;

        //transform sigma back into original basis
        step = -eiVec*sigma;

        if(print>100)cout << "step: "  << step << endl;
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
            if(print>10)cout << " acceleration: "  << acc_fac<<endl;
            step*=acc_fac;
        }
        step_norm2 = step.transpose()*step ;
        first_order = this->gradient.transpose()*step ;
        second_order = step.transpose()*this->hessian*step;
        if(print>10)cout << " gradient magnitude: "  << first_order*first_order/step_norm2<<endl;

        fac=0.7*std::abs(first_order/second_order)*10*std::sqrt(step_norm2);
        if(fac>100.0)fac=100.0;
        if(fac<0.00001)fac=0.00001;

        //Newton step when all eigenvalues are <0, and gradient/h sufficiently small
        newton_step=0;
        if(mu<mu_min*1.1){
            if(print>10)cout << "Taking Newton step  "<<  endl;
            newton_step=1;
            N_newton_step++;
            if(N_newton_step>2&&(step_norm2<1.E-3||newton_step_exact==1)){
                if(print==2)cout<<"  Newton exact ";
                newton_step_exact=1;
                mu_Newton*=0.2;//level shift for positive and zero eigenvalues
                gradient_tmp=gradient;
                Hess_tmp=hessian;
                for (i=0; i<N2h; i++) {
                    Hess_tmp(i,i)-=mu_Newton;
                }
                //EIGEN SOLVER
                //      step=-Hess_tmp.colPivHouseholderQr().solve(gradient_tmp);
                step=-Hess_tmp.ldlt().solve(gradient_tmp);

                step_norm2 = step.transpose()*step ;
                first_order = this->gradient.transpose()*step ;
                second_order = step.transpose()*this->hessian*step;
            }else{
                if(print==2)cout<<" Newton approx ";
                newton_step_exact=0;
                mu_Newton=mu_Newton_init;
            }

            if(print==20)cout<<endl<<" 2nd  "<<second_order<<" fi*fi/d  "<<x<<" s g "<<step.transpose()*gradient<<endl;
        }else{
            if(print==2)cout<<"    restricted ";
            N_newton_step=0;
            newton_step_exact=0;
            mu_Newton=mu_Newton_init;
        }

        //    if(print==2){cout<<setw(22)<< step_norm2;}
        if(print==2)cout<< setiosflags(ios::fixed) <<setprecision(3)<< setw(10)<<step_norm2;
        if(print>100){
            cout << "step : "  ;
            for (i=0; i<N2h; i++) {
                cout << " "<< step(i) ;
            }
            cout << endl ;
        }

        expected_change=first_order +0.5*second_order;

        if(print>10)  cout << "step size  "<< std::sqrt(step_norm2)   <<endl;
        if(print>10)cout << "expected first, second order and total change in r*r  ";
        if(print>10)cout << first_order<< " " << 0.5*second_order <<" " <<expected_change  <<endl;

        this->do_step(step);
        value_functional = this->functional();

        //    if(print==2){cout<< setw(22)<< value_functional<< setw(22)<< expected_change;}
        //    if(print==2){cout<< setw(22)<< value_functional-value_functional_old;}
        if(print==2)cout<< setprecision(12)<< setw(15)<< value_functional;
        if(print==2)cout<< setprecision(4)<< setw(11)<< expected_change;
        if(print==2)cout<< setprecision(4)<< setw(11)<< value_functional-value_functional_old;

        if(print>10) cout <<" r*r  "<< value_functional << " change in r*r  "<< value_functional-value_functional_old << endl;

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
            //    if(h>20.0)h=20.0;
        }
        if(print>10)cout <<"trust radius set  to "<< h << " test: "<<relative_change << " mu: "<<mu<<  " maxeival: "<<maxEiVal <<endl;
        if(print>10)cout << "gradient norm " << gradient_norm<< endl;

        if(gradient_norm<threshold&&maxEiVal<10*std::sqrt(std::abs(threshold))){
            //finished
            converged=true;
        }
        if(print==2&&!wrongstep){cout<<setw(10)<< dcount;}
        if(print==2)cout<<endl;

    }//iterations

    if (iter < maxIter) {
        iter = iter - 1;
    } else {
        iter = -1;
    }
    return iter;
}

} //namespace mrchem
