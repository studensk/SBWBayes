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

//Define Normal cdf (log)
TMB_ATOMIC_VECTOR_FUNCTION(
  pnorm_log1,
  1,
  ty[0] = atomic::Rmath::Rf_pnorm5(tx[0],0,1,1,1);
,
Type W  = ty[0];
Type DW = dnorm1(tx[0])/exp(W);
px[0] = DW * py[0];
)
  
  template<class Type> 
  Type pnorm_log1(Type x){
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return pnorm_log1(tx)[0];
  }  
  
  // Data and parameters
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    DATA_IVECTOR(nobs);
    DATA_IVECTOR(stage);
    DATA_IVECTOR(t_block1);
    DATA_IVECTOR(t_block2);
    DATA_VECTOR(k_vec);

    DATA_VECTOR(time1);
    DATA_VECTOR(time2);
    DATA_VECTOR(temp1);
    DATA_VECTOR(temp2);
    DATA_VECTOR(time2d);
    DATA_VECTOR(time1d);
    DATA_INTEGER(use_prior);
    
    PARAMETER(HA);
    PARAMETER(TL);
    PARAMETER(HL);
    PARAMETER(TH);
    PARAMETER(HH);
    PARAMETER_VECTOR(rho);
    PARAMETER_VECTOR(s_eps);
    PARAMETER_VECTOR(s_upsilon);
    PARAMETER_VECTOR(upsilon);
    
    
    Type jnll = 0;
    Type tpred1;
    Type tpred2;
    Type epsm1;
    Type epsij;
    Type pnormdiff;
    Type nlogp;
    Type TA;
    Type epsm1_std;
    Type epsij_std;
    Type c_upsilon1;
    Type c_upsilon2;
    Type tdiff;
    
    // Loop over observations
    for (int i=0; i<time1.size(); i++) {

      // Calculate rho from quadratic function of dev stage and
      //   transform for bias reduction

      TA = calc_TA(-HL, HH, TL, TH);

      // Calculate dev times for both sustainable and lethal temps
      tpred1 = calc_pred(temp1(i), rho(stage(i)), HA, TL, -HL, TH, HH, TA);
      tpred2 = calc_pred(temp2(i), rho(stage(i)), HA, TL, -HL, TH, HH, TA);

      c_upsilon1 = exp(upsilon(t_block1(i))*s_upsilon(stage(i)));
      c_upsilon2 = exp(upsilon(t_block2(i))*s_upsilon(stage(i)));

      // Transform for bias reduction
      tpred1 *= c_upsilon1;
      tpred2 *= c_upsilon2;

      // Calculate and standardize observed values of epsilon
      epsm1 = log(time1d(i)/tpred1 + time2d(i)/tpred2);
      epsij = log(time1(i)/tpred1 + time2(i)/tpred2);

      epsm1_std = epsm1/s_eps(stage(i));
      epsij_std = epsij/s_eps(stage(i));

      Type pnorm_ij = pnorm_log1(epsij_std);
      Type pnorm_m1 = pnorm_log1(epsm1_std);

      Type timesum = time1d(i) + time2d(i);
      if (timesum == 0) {
        pnormdiff = pnorm_ij;
      }
      else {
        pnormdiff =  logspace_sub(pnorm_ij, pnorm_m1);
      }

      // Subtract log prob times number of observations
      nlogp = nobs(i)*pnormdiff;
      jnll -= nlogp;

      // if (stage(i) == 2) {
      //   r0 = 1/calc_pred(Type(0), rho(stage(i)), HA, TL, -HL, TH, HH, TA);
      //   r40 = 1/calc_pred(Type(40), rho(stage(i)), HA, TL, -HL, TH, HH, TA);
      //   jnll -= dexp(r0, Type(20), 1);
      //   jnll -= dexp(r40, Type(20), 1);
      // }

    }
    
    // Subtract prior probabilities
    if (use_prior == 1) {
      tdiff = TH - TL;
      jnll -= sum(dnorm(upsilon, Type(0), Type(1), 1));
      
      for (int j=0; j<k_vec.size(); j++) {
        jnll -= dgamma(rho(j), k_vec(j), Type(0.045), 1);
      }
      
      jnll -= dgamma(HL, Type(3.6), Type(2.253), 1);
      jnll -= dgamma(HA, Type(5.4), Type(0.134), 1);
      jnll -= dgamma(HH, Type(7.6), Type(3.12), 1);
      
      jnll -= dnorm(log(TL), Type(5.64), Type(0.0067), 1);
      jnll -= dgamma(tdiff, Type(112), Type(0.226), 1);
      
      jnll -= sum(dnorm(log(s_eps), Type(-1.5), Type(0.1), 1));
      jnll -= sum(dnorm(log(s_upsilon), Type(-2.5), Type(0.05), 1));
    }
    return jnll;
    
  }
  