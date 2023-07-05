// Base dependencies
#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>
#include <RcppArmadillo.h>
// #include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
double CPP_ESA_C(const arma::mat& d2mask, 
				const double& sigma,
				const double& g0,
				const double& area,
				const int& detfn)
{
	
	int n_traps = d2mask.n_cols;
	int n_mask = d2mask.n_rows;
	arma::mat p_trap(n_mask, n_traps);
	arma::vec p_capt = arma::ones(n_mask);
	
	
	if(detfn == 1){
		p_trap =  g0*exp(-d2mask/(2*sigma*sigma));	// Make it easy to update the detection function here.
	}else{
		p_trap = 1 - exp(-(g0*exp(-d2mask/(2*sigma*sigma))));
	}
	
	arma::vec miss = arma::ones(n_mask);		
	for(int j = 0; j < n_traps; j++)
	{
		miss %= (1 - p_trap.col(j));
	}
		
	p_capt -= miss;
	double ESA = sum(p_capt)*area;
	return ESA;
}

// [[Rcpp::export]]
double CPP_ESA_A(const arma::mat& d2mask, 
				const double& lambda,
				const double& sigma,
				const double& g0,
				const double& area,
				const int& detfn)
{
	
	int n_traps = d2mask.n_cols;
	int n_mask = d2mask.n_rows;
	arma::mat p_trap(n_mask, n_traps);
	arma::vec p_capt = arma::ones(n_mask);
	
	
	if(detfn == 1){
		p_trap =  g0*exp(-d2mask/(2*sigma*sigma));	// Make it easy to update the detection function here.
	}else{
		p_trap = 1 - exp(-(g0*exp(-d2mask/(2*sigma*sigma))));
	}
	
	arma::vec miss = arma::ones(n_mask);		
	for(int j = 0; j < n_traps; j++)
	{
		miss %= (1 - p_trap.col(j));
	}
		
	p_capt -= miss;
	double ESA = sum((1-exp(-p_capt*lambda)))*area;
	return ESA;
}

// [[Rcpp::export]]
double CPP_ESA_Cam(const arma::mat& d2mask, 
				const double& lambda,
				const double& sigma,
				const double& area)
{
	
	int n_traps = d2mask.n_cols;
	int n_mask = d2mask.n_rows;
	arma::mat p_trap(n_mask, n_traps);
	
	
	p_trap =  lambda*exp(-d2mask/(2*sigma*sigma));	// Make it easy to update the detection function here.
	
	arma::vec Hk = arma::zeros(n_mask);		
	for(int j = 0; j < n_traps; j++)
	{
		Hk += p_trap.col(j);
	}
		
	double ESA = sum(1-exp(-Hk*lambda))*area;
	return ESA;
}