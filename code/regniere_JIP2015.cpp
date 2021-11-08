#include <TMB.hpp> 
#include <iostream>
#include <string>

// Function to calculate TA from other model parameters
template<class Type>
Type calc_TA(Type HL, Type HH, Type TL, Type TH) {
  Type c3 = Type(0.0001987);
  Type den, tphi;
  den = c3*log(-HL/HH) + (HL/TL) - (HH/TH);
  tphi = (HL - HH)/den;
  return tphi;
}

// Development time function
template<class Type>
Type calc_pred(Type temp, Type rho25, Type HA, Type TL, Type HL, Type TH, Type HH, Type TA) {
	Type c1 = Type(273);
  Type c2 = TA;
	Type c3 = Type(0.0001987);
	Type tK, num, den1, den2, tau;
	
	tK = temp + c1;
	num = rho25 * tK/c2 * exp(HA/c3* (1/c2-1/tK));
	den1 = exp(HL/c3 *(1/TL-1/tK));
	den2 = exp(HH/c3 *(1/TH-1/tK));
	tau = 1/(num/(1 + den1 + den2));
	return tau;
}

// Define Normal pdf
template<class Type>
Type dnorm1(Type x){
  return Type(1.0/sqrt(2.0*M_PI)) * exp(-Type(.5)*x*x);
}

// Define Normal cdf (log)
TMB_ATOMIC_VECTOR_FUNCTION(
  pnorm_log,
  1,
  ty[0] = atomic::Rmath::Rf_pnorm5(tx[0],0,1,1,1);
,
Type W  = ty[0];                    
Type DW = dnorm1(tx[0])/exp(W);
px[0] = DW * py[0];
)
  
template<class Type> 
Type pnorm_log(Type x){
  CppAD::vector<Type> tx(1);
  tx[0] = x;
  return pnorm_log(tx)[0];
}  

// Define cauchy quantile function  
template<class Type>
Type qcauchy(Type y, Type location, Type scale){
  Type q = location + scale*tan(M_PI*(y - Type(0.5)));
  return q;
}

// Data and parameters  
template<class Type>
Type objective_function<Type>::operator() ()
{
	DATA_IVECTOR(nobs);
	DATA_IVECTOR(block);
	DATA_VECTOR(time1);
	DATA_VECTOR(time2);
	DATA_VECTOR(temp1);
	DATA_VECTOR(temp2);
	DATA_VECTOR(time2d);
	DATA_INTEGER(use_prior);

	PARAMETER(rho25);
	PARAMETER(HA);
	PARAMETER(TL);
	PARAMETER(HL);
	PARAMETER(TH);
	PARAMETER(HH);
	PARAMETER(s_eps);
	PARAMETER(s_upsilon);
	PARAMETER_VECTOR(upsilon);

	
	Type jnll = 0;
	Type tpred1;
	Type tpred2;
	Type epsm1;
	Type epsij;
	Type pnormdiff;
	Type nlogp;
	Type epsm1_std;
	Type epsij_std;
	Type TA;
	Type u_upsilon;
	Type c_upsilon;

	// Loop over observations
	for (int i=0; i<time1.size(); i++) {
		
		TA = calc_TA(-HL, HH, TL, TH);
	  
	  // Calculate dev times for both sustainable and lethal temps
		tpred1 = calc_pred(temp1(i), rho25, HA, TL, -HL, TH, HH, TA);
		tpred2 = calc_pred(temp2(i), rho25, HA, TL, -HL, TH, HH, TA);
		
		// Quantile match Cauchy
		u_upsilon = pnorm(upsilon(block(i)), Type(0), Type(1));
		c_upsilon = qcauchy(u_upsilon, Type(1), s_upsilon);
		
		// Transform for bias reduction
		tpred1 *= c_upsilon;
		tpred2 *= c_upsilon;

		// Calculate and standardize observed values of epsilon
		epsm1 = log(time1(i)/tpred1 + time2d(i)/tpred2);
		epsij = log(time1(i)/tpred1 + time2(i)/tpred2);
		
		epsm1_std = epsm1/s_eps;
		epsij_std = epsij/s_eps;

		// Calculate log probability
		Type pnorm_ij = pnorm_log(epsij_std);
		Type pnorm_m1 = pnorm_log(epsm1_std);
		
		if (time2d(i) == 0) {
		  pnormdiff = pnorm_ij;
		}
		else {
		  pnormdiff =  logspace_sub(pnorm_ij, pnorm_m1);
		}
		
		// Subtract log prob times number of observations
		nlogp = nobs(i)*pnormdiff;

		
		jnll -= nlogp;

	} 
	// Subtract prior probabilities
  if (use_prior == 1) {
    jnll -= sum(dnorm(upsilon, Type(0), Type(1), 1));
    
    jnll -= dgamma(rho25, Type(2), Type(0.25), 1); 
    
    jnll -= dgamma(HL, Type(6), Type(2), 1);
    jnll -= dgamma(HA, Type(5), Type(0.2), 1);
    jnll -= dgamma(HH, Type(10), Type(3), 1);
    
    jnll -= dnorm(TL, Type(284), Type(2), 1);
    jnll -= dnorm(TH, Type(304), Type(2), 1);
    
    jnll -= dnorm(log(s_eps), Type(-1.5), Type(0.1), 1);
    jnll -= dnorm(log(s_upsilon), Type(-2.5), Type(0.05), 1);
  }
	
	return jnll;

}
