#include <iostream>
#include <fstream>
#include <list>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <ctime>
#include <Rcpp.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;
using namespace Rcpp;


#define ARMA_DONT_PRINT_ERRORS

/*
// define a class 
class model_params{
public:
  const int max_iter = 0;
  const int min_iter = 30;
  const int num_rept = 10;
  int init_select = 1;
  const double diff_err = 1e-4;
  const double mean_fitratio = 0.1;
};
*/

void printProgBar( int percent ) {
  string bar;
  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }

  cout<< "\r" "[" << bar << "] ";
  cout.width( 3 );
  cout<< percent << "%     " << std::flush;
}

// function to compute objective function
double cal_obj_func(mat X, mat U, mat V, mat L1, mat L2){
	mat residual = X - U*V.t();
	// *NMF term
	double obj_value = sum(sum(residual%residual));
	
	// *add penalty terms
	obj_value += sum(sum((V.t()*L1)%V.t()));
	obj_value += sum(sum((U.t()*L2)%U.t()));
	
	return(obj_value);
}

//************************************************************//
//  Single-Omic Nonnegative Matrix Factorization using ADMM   //
//************************************************************//
// [[Rcpp::export]]
SEXP SOMF_cpp(SEXP Xin, SEXP Win, SEXP Hin, SEXP betain, SEXP rhoin, SEXP num_iterin){// *
	try {
		mat X = as<mat>(Xin);  // *dim = num_samp x num_fea
		mat W = as<mat>(Win);  // *dim = num_samp x num_pc
		mat H = as<mat>(Hin);  // *dim = num_pc x num_fea
 ///
		// *coefficient of regularization
		double rho = Rcpp::as<double>(rhoin);
		int beta = Rcpp::as<int>(betain);
		int num_iter = Rcpp::as<int>(num_iterin);
		
		int num_pc = W.n_cols;
		mat E(num_pc, num_pc, fill::eye);
		//int num_samp = X.n_cols;
		//int num_fea = X.n_rows;
		
		// *initial values
		mat Xprime = W*H;
		mat Wplus = W;
		mat Hplus = H;
		
		mat Xdual = zeros<mat>(size(Xprime));
		mat Wdual = zeros<mat>(size(W));
		mat Hdual = zeros<mat>(size(H));
		
		int iter = 0;
		mat ZH = zeros<mat>(size(H));
		mat ZW = zeros<mat>(size(W));
		
		while(iter < num_iter){
			//cout<<"SONMF:: iter  = "<< iter <<endl;
			printProgBar(( (iter+1) / (float) num_iter) * (float) 100);
			// *updating H
			H = inv_sympd(W.t()*W + E) * (W.t()*Xprime + Hplus + 1/rho *( W.t()*Xdual - Hdual) );
			
			
			// *updating W
			mat P = H*H.t() + E;
			W = (Xprime*H.t() + Wplus + 1/rho *(Xdual*H.t() - Wdual)) * inv_sympd(P);
			
			// *updating X
			mat Xtmp = W*H;
			if(beta == 1){
				mat B = rho*Xtmp - Xdual -1;
				Xprime = ( B + sqrt(B%B + 4.0*rho*X) )/(2.0*rho);
			}else if(beta == 0){
				mat A = Xdual/rho - Xtmp;
				mat B = 1.0/(3*rho) - A%A/9.0;
				mat C = -A%A%A/27.0 + A/(6.0*rho) + X/(2.0*rho);
				mat D = B%B%B + C%C;
				/*
				uvec idxlaz = find(D>=0);
				uvec idxlez = find(D<0);
				Xprime.elem(idxlaz) = pow( C.elem(idxlaz) + sqrt(D.elem(idxlaz)), 1.0/3.0) + 
									  pow( C.elem(idxlaz) - sqrt(D.elem(idxlaz)), 1.0/3.0) - 
									  A.elem(idxlaz)/3.0;
				
				Xprime.elem(idxlez) = 2.0*sqrt(-B.elem(idxlez))%cos(acos(C.elem(idxlez)/(pow(-B.elem(idxlez),3.0/2.0)) )/3.0 ) - 
				                      A.elem(idxlez)/3.0;
				*/
				for(size_t i=0;i<X.n_rows;i++){
					for(size_t j=0;j<X.n_cols;j++){
						if(D(i,j)>=0){
							Xprime(i,j) = cbrt( C(i,j) + sqrt(D(i,j))) +
										  cbrt( C(i,j) - sqrt(D(i,j))) -
										  A(i,j)/3.0;
						}else{
							Xprime(i,j) = 2.0*sqrt( -B(i,j))*
							              cos( acos(C(i,j)/sqrt(-B(i,j)*B(i,j)*B(i,j))) /3.0 ) - 
				                           A(i,j)/3.0;
						}
						
					}
				}
			}else if(beta == 2){
				cout<<"add it later ..." <<endl;
			}
			
			// *updating Hplus, Wplus
			Hplus = arma::max( H + 1.0/rho*Hdual, ZH);
			Wplus = arma::max( W + 1.0/rho*Wdual, ZW);
			
			// *updating dual variables
			Xdual = Xdual + rho*(Xprime - Xtmp);
			Hdual = Hdual + rho*(H - Hplus);
			Wdual = Wdual + rho*(W - Wplus);
			
			// *update iterator
			iter++;
		}
		
		cout<<endl;
		return List::create(Named("W") = Wplus, Named("H") = Hplus, Named("iteration") = iter );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}


// ***********************************************************//
//  Multi-Omic Nonnegative Matrix Factorization using ADMM   //
//***********************************************************//
// *define structure for each data set
class scXWClass{
	public:
	mat X;
	mat Xprime;
	mat Xtmp;
	
	rowvec zero_ratio; // *num_sample x 1
	
	mat W;
	
	double Weta; // *l1 penalty parameter
	mat Wzero;
	
	mat Wplus;
	
	mat Xdual;
	mat Wdual;
	
};

// *unified gene expression weight U
class scUClass{
	public:
	mat U;

	mat Uzero;
	
	mat Uplus;
	mat Udual;
};


// function to normalization
void normal_wu(struct scXWClass *scXW, struct scUClass &scU, int num_data){
	
	// *U dim k x num_sample 
	colvec norms = sqrt(sum(scU.Uplus%scU.Uplus,1)); // 1 x num_pc
	
	norms.elem( find(norms <1e-15) ).fill(1e-15);

	mat norms_factorU = ones<mat>( size(scU.Uplus) );
	// *H.n_cols = num_pc
	for(size_t i=0;i<scU.Uplus.n_rows;i++){
		norms_factorU.row(i) = norms_factorU.row(i)*(1/norms[i]);
	}
	scU.Uplus = scU.Uplus%norms_factorU;
	
	for(size_t id=0; id<num_data; id++){
		mat norms_factorW = ones<mat>( size(scXW[id].Wplus) );
		for(size_t i=0;i<scU.Uplus.n_rows;i++){
			norms_factorW.col(i) = norms_factorW.col(i)*norms[i];
		}
		scXW[id].Wplus = scXW[id].Wplus%norms_factorW;
	}// *end for id

	//return 0;
}


double cal_obj_func(struct scXWClass *scXW, struct scUClass scU, int num_data, int beta){
	
	double obj_value = 0.0;
	
	for(size_t id=0; id<num_data; id++){
		mat Xappr = scXW[id].Wplus * scU.Uplus;
		mat ratio_XX = scXW[id].X/Xappr;
	
		if(beta == 2){ // *euclidean distance
			obj_value += accu( (scXW[id].X - Xappr) % (scXW[id].X - Xappr) );
		}else if(beta == 1){ // *Kullback-Leibler divergence
		
			for(size_t i=0;i<Xappr.n_rows;i++){
				for(size_t j=0;j<Xappr.n_cols;j++){
					if(std::isfinite(ratio_XX(i,j)) && ratio_XX(i,j) != 0){
						obj_value += Xappr(i,j) * log(ratio_XX(i,j)) - scXW[id].X(i,j) + 1;
					}
				}
			}

		}else if(beta == 0){ // *Itakura-Saito divergence
			for(size_t i=0;i<Xappr.n_rows;i++){
				for(size_t j=0;j<Xappr.n_cols;j++){
					if(std::isfinite(ratio_XX(i,j)) && ratio_XX(i,j) != 0){
						obj_value += ratio_XX(i,j) - log(ratio_XX(i,j)) - 1;
					}
				}
			}
		}// end fi
	}// *end for id
	return obj_value;
}// end func


void update_w(scXWClass &scXW, scUClass scU, mat E, double rho, int beta){
	
	// *updating W with unified U
	mat P = scU.U*scU.U.t() + E;
	scXW.W = (scXW.Xprime*scU.U.t() + scXW.Wplus + 1/rho *(scXW.Xdual*scU.U.t() - scXW.Wdual)) * inv_sympd(P);
	//return 0;
}

void update_x(scXWClass &scXW, scUClass scU, mat E, double rho, int beta){
	
	// *updating X with unified U
	scXW.Xtmp = scXW.W*scU.U;
	if(beta == 1){ // *Kullback-Leibler divergence
		mat B = rho*scXW.Xtmp - scXW.Xdual -1;
		scXW.Xprime = ( B + sqrt(B%B + 4.0*rho*scXW.X) )/(2.0*rho);
	}else if(beta == 0){// *Itakura-Saito divergence
		mat A = scXW.Xdual/rho - scXW.Xtmp;
		mat B = 1.0/(3*rho) - A%A/9.0;
		mat C = -A%A%A/27.0 + A/(6.0*rho) + scXW.X/(2.0*rho);
		mat D = B%B%B + C%C;
			
		for(size_t i=0;i<scXW.X.n_rows;i++){
			for(size_t j=0;j<scXW.X.n_cols;j++){
				if(D(i,j)>=0){
					scXW.Xprime(i,j) = cbrt( C(i,j) + sqrt(D(i,j))) +
								  cbrt( C(i,j) - sqrt(D(i,j))) -
								  A(i,j)/3.0;
				}else{
					scXW.Xprime(i,j) = 2.0*sqrt( -B(i,j))*
							      cos( acos(C(i,j)/sqrt(-B(i,j)*B(i,j)*B(i,j))) /3.0 ) - 
				                   A(i,j)/3.0;
				}// *end if
						
			}// *end for j
		}// *end for i
	}else if(beta == 2){
		cout<<"add it later ..." <<endl;
	}

	//return 0;
}
void update_u(scXWClass *scXW, scUClass &scU, mat E, int num_data, double rho){
	
	// *updating U
	//cout<<"update::U before = "<<scU.U[1,2]<<endl;
	mat WtW = E;
	scU.U = scU.Uplus - 1/rho * scU.Udual;
	for(size_t i=0; i<num_data; i++){
		WtW += scXW[i].W.t()*scXW[i].W;
		//scU.U += scXH[i].H.t()*scXH[i].Xprime - scXH[i].H.t()*scXH[i].H*scXH[i].Wplus + 1/rho *scXH[i].H.t()*scXH[i].Xdual; 
		scU.U += scXW[i].W.t()*scXW[i].Xprime + 1/rho *scXW[i].W.t()*scXW[i].Xdual; 
	}
	scU.U = inv_sympd(WtW) * scU.U;
	//cout<<"inv_sympd(HtH) = "<<inv_sympd(HtH)<<endl;
	//cout<<"update::U after = "<<scU.U[1,2]<<endl;
	// return 0;
}

// [[Rcpp::export]]
SEXP MOMF_cpp(SEXP Xin, SEXP Win, SEXP Uin, SEXP betain, SEXP rhoin, SEXP num_datain, SEXP num_iterin){// *
try{
		// *convert data format
		const Rcpp::List multiX(Xin);  // *dim = num_samp x num_fea x num_omic, data set
		Rcpp::List multiW(Win);        // *dim = num_samp x num_pc x num_omic, initial values
		//mat U = as<mat>(Uin);        // *dim = num_pc x num_fea
		
		// *coefficient of regularization
		double rho = Rcpp::as<double>(rhoin);
		int beta = Rcpp::as<int>(betain);
		int num_iter = Rcpp::as<int>(num_iterin);
		int num_data = Rcpp::as<int>(num_datain);
		//int num_pc = Rcpp::as<int>(num_pcin); // *all datasets should share the same number of components
		
		// *convert to struct
		scXWClass scXW[num_data];
		
		for(size_t i=0; i<num_data; i++){
			stringstream datX, datW;
			datX <<"X"<<(i+1);
			scXW[i].X = as<mat>(multiX[datX.str()]);

			datW <<"W"<<(i+1);
			scXW[i].W = as<mat>(multiW[datW.str()]);

			scXW[i].Wzero = zeros<mat>(size(scXW[i].W));
		}// *end for i		
		
		
		scUClass scU;
		scU.U = as<mat>(Uin);  
		scU.Udual = zeros<mat>(size(scU.U));
		//scU.Uzero = zeros<mat>(size(scU.U));
		scU.Uplus = scU.U;
		// *initial values
		int num_pc = scU.U.n_rows;
		mat E(num_pc, num_pc, fill::eye);
		
		for(size_t i=0; i<num_data; i++){
			scXW[i].Xprime = scXW[i].W*scU.U;
			scXW[i].Wplus = scXW[i].W;
		
			scXW[i].Xdual = zeros<mat>(size(scXW[i].X));
			scXW[i].Wdual = zeros<mat>(size(scXW[i].W));
		}
		
		
		// *main loop iteration
		int iter = 0;
		vec hist_obj = zeros<vec>(num_iter);
		while(iter < num_iter){
			//cout<<"monmf_admm:: iter  = "<< iter <<endl;	
			printProgBar(( (iter+1) / (float) num_iter) * (float) 100);
			// *updating H for each data set
			for(size_t i=0; i<num_data; i++){
				//cout<<"H before = "<<scXH[i].H[1,2]<<endl;
				update_w(scXW[i], scU, E, rho, beta);
				//cout<<"H after = "<<scXH[i].H[1,2]<<endl;
			}
			
			// *updating U
			
			update_u(scXW, scU, E, num_data, rho);
				
			
			// *updating X for each data set
			for(size_t i=0; i<num_data; i++){
				update_x(scXW[i], scU, E, rho, beta);
			}
			
			
			// *updating dual variables, Xdual, Hdual, Udual
			for(size_t i=0; i<num_data; i++){
				scXW[i].Xdual = scXW[i].Xdual + rho*(scXW[i].Xprime - scXW[i].Xtmp);
				scXW[i].Wdual = scXW[i].Wdual + rho*(scXW[i].W - scXW[i].Wplus);
			}
			scU.Udual = scU.Udual + rho*(scU.U - scU.Uplus);
			
			// *updating Hplus and  Uplus
			for(size_t i=0; i<num_data; i++){
				scXW[i].Wplus = arma::max( scXW[i].W + 1.0/rho*scXW[i].Wdual, scXW[i].Wzero);
			}
			scU.Uplus = arma::max(scU.U + 1.0/rho*scU.Udual, scU.Uzero);
			
			// *normalize the Hplus, and Uplus
			//normal_hu(scXH, scU, num_data);
			// *stop rule
			hist_obj(iter) = cal_obj_func(scXW, scU, num_data, beta);
			// *update iterator
			iter++;
		}
		cout<<endl;
		// *restore back
		for(size_t i=0; i<num_data; i++){
			stringstream datW;
			datW <<"W"<<(i+1);
			multiW[datW.str()] = scXW[i].Wplus;
		}// *end for i	
		return List::create(Named("U") = scU.Uplus, Named("multiW") = Rcpp::wrap(multiW), Named("iteration") = iter, Named("hist_obj") = hist_obj );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}// funcs


// [[Rcpp::export]]
SEXP MOMF_fixU_cpp(SEXP Xin, SEXP Win, SEXP Uin, SEXP Uplusin, SEXP betain, SEXP rhoin, SEXP num_datain, SEXP num_iterin){// *
try{
		// *convert data format
		const Rcpp::List multiX(Xin);  // *dim = num_samp x num_fea x num_omic, data set
		Rcpp::List multiW(Win);        // *dim = num_samp x num_pc x num_omic, initial values
		//mat U = as<mat>(Uin);        // *dim = num_pc x num_fea
		
		// *coefficient of regularization
		double rho = Rcpp::as<double>(rhoin);
		int beta = Rcpp::as<int>(betain);
		int num_iter = Rcpp::as<int>(num_iterin);
		int num_data = Rcpp::as<int>(num_datain);
		//int num_pc = Rcpp::as<int>(num_pcin); // *all datasets should share the same number of components
		
		// *convert to struct
		scXWClass scXW[num_data];
		
		for(size_t i=0; i<num_data; i++){
			stringstream datX, datW;
			datX <<"X"<<(i+1);
			scXW[i].X = as<mat>(multiX[datX.str()]);

			datW <<"W"<<(i+1);
			scXW[i].W = as<mat>(multiW[datW.str()]);

			scXW[i].Wzero = zeros<mat>(size(scXW[i].W));
		}// *end for i		
		
		
		scUClass scU;
		scU.U = as<mat>(Uin);  
		scU.Udual = zeros<mat>(size(scU.U));
		scU.Uzero = zeros<mat>(size(scU.U));
		//scU.Uplus = scU.U;
		scU.Uplus = as<mat>(Uplusin);
		// *initial values
		int num_pc = scU.U.n_rows;
		mat E(num_pc, num_pc, fill::eye);
		
		for(size_t i=0; i<num_data; i++){
			scXW[i].Xprime = scXW[i].W*scU.U;
			//scXW[i].Wplus = scXW[i].W;
			scXW[i].Wplus = arma::max(scXW[i].W, scXW[i].Wzero);
		
			scXW[i].Xdual = zeros<mat>(size(scXW[i].X));
			scXW[i].Wdual = zeros<mat>(size(scXW[i].W));
		}
		
		
		// *main loop iteration
		int iter = 0;
		// vec hist_obj = zeros<vec>(num_iter);
		while(iter < num_iter){
			//cout<<"monmf_admm:: iter  = "<< iter <<endl;	
			printProgBar(( (iter+1) / (float) num_iter) * (float) 100);
			// *updating H for each data set
			for(size_t i=0; i<num_data; i++){
				//cout<<"H before = "<<scXH[i].H[1,2]<<endl;
				update_w(scXW[i], scU, E, rho, beta);
				//cout<<"H after = "<<scXH[i].H[1,2]<<endl;
			}
			
			// *updating U
			
			update_u(scXW, scU, E, num_data, rho);
				
			
			// *updating X for each data set
			for(size_t i=0; i<num_data; i++){
				update_x(scXW[i], scU, E, rho, beta);
			}
			
			
			// *updating dual variables, Xdual, Hdual, Udual
			for(size_t i=0; i<num_data; i++){
				scXW[i].Xdual = scXW[i].Xdual + rho*(scXW[i].Xprime - scXW[i].Xtmp);
				scXW[i].Wdual = scXW[i].Wdual + rho*(scXW[i].W - scXW[i].Wplus);
			}// end for
			scU.Udual = scU.Udual + rho*(scU.U - scU.Uplus);
			
			// *updating Hplus and  Uplus
			for(size_t i=0; i<num_data; i++){
				scXW[i].Wplus = arma::max( scXW[i].W + 1.0/rho*scXW[i].Wdual, scXW[i].Wzero);
			}// end for
			
			//==============================
			// *because we have provided Uplus from prior knowledge, no need to update further
			//scU.Uplus = arma::max(scU.U + 1.0/rho*scU.Udual, scU.Uzero);
			//=============================================

			// *normalize the Hplus, and Uplus
			//normal_hu(scXH, scU, num_data);
			// *stop rule
			// hist_obj(iter) = cal_obj_func(scXW, scU, num_data, beta);
			// cout<<"hist_obj = " << hist_obj(iter)<<endl;
			// *update iterator
			iter++;
		}
		cout<<endl;
		// *restore back
		for(size_t i=0; i<num_data; i++){
			stringstream datW;
			datW <<"W"<<(i+1);
			multiW[datW.str()] = scXW[i].Wplus;
		}// *end for i	
		return List::create(Named("U") = scU.U+1.0/rho*scU.Udual, Named("multiW") = Rcpp::wrap(multiW), Named("iteration") = iter );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}// end func


/////////////////////////////////////////////////////////////////////////////////////////
//                             CODE END HERE                                           //
/////////////////////////////////////////////////////////////////////////////////////////











