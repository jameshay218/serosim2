#include <string>
#include <iostream>
#include <fstream>
#include <Rmath.h>
#include <math.h>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

//' Calculate likelihood
//'
//' Calculates the likelihood of a given set of observed titres given predicted titres. Based on truncated discritised normal.
//' @param expected NumericVector, as returned by \code{\link{infection_model_indiv}}
//' @param data NumericVector, the vector of observed titres
//' @param theta NumericVector, the vector of named model parameters, requiring MAX_TITRE and error
//' @return a single log likelihood
//' @export
//[[Rcpp::export]]
double likelihood_titre(NumericVector expected, NumericVector data, double error, int MAX_TITRE){
  int n = expected.size();
  double lnlike = 0;
  for(int i=0; i < n; ++i){
    if(!NumericVector::is_na(data[i])){
      if(data[i] > MAX_TITRE){
	lnlike += R::pnorm(MAX_TITRE, expected[i], error,0,1);
      } else if(data[i] <= 0){
	lnlike += R::pnorm(1, expected[i], error,1,1);
      } else {
	lnlike += log(R::pnorm(data[i]+1, expected[i],error, 1,0) - 
		      R::pnorm(data[i], expected[i], error,1,0));
      }
    }
  }
  return(lnlike);
}

//[[Rcpp::export]]
NumericVector vector_sort(NumericVector x){
  NumericVector y = clone(x);
  return(y.sort());
}
//[[Rcpp::export]]
IntegerVector vector_order(NumericVector x){
  NumericVector y = clone(x).sort();
  return(match(y,x));
}

//[[Rcpp::export]]
NumericVector individual_sim(
			     NumericVector mu_pars, 
			     NumericVector tp_pars, 
			     NumericVector m_pars, 
			     NumericVector ti_pars, 
			     double  y0b,
			     double lower_titre_bound, 
			     NumericVector times
			     ){
  double y0 =  y0b;
  double final_t, t_i, mu, tp, m, tmp, max_t, tmp2;
  NumericVector::iterator t = times.begin();
  NumericVector Y(times.size(), y0);

  int j = 0;
  int no_infections = ti_pars.size();
  max_t = times[times.size()-1];
  
  for(int i = 1; i <= no_infections; ++i){
    /*while(ti_pars[i-1] < times[0]){
      ++i;
      if(i > no_infections) break;
    }
    */
    if(i == no_infections){
      final_t = max_t;
    }
    else {
      final_t = ti_pars[i];
    }
    t_i = ti_pars[i-1];
    mu = mu_pars[i-1];
    tp = tp_pars[i-1];
    m = m_pars[i-1];
    
    while(t != times.end() && *t <= final_t){
      tmp = 0;
      if(*t <= t_i) tmp = y0;
      else if(*t > t_i && *t <= (t_i + tp)) tmp = (mu/tp)*(*t) - (mu/tp)*(t_i) + y0;
      else tmp = -m*((*t)) + m*(t_i + tp) + mu + y0;
      Y[j] = tmp;
      
      if(Y[j] < lower_titre_bound) Y[j] = lower_titre_bound;
      ++t;
      ++j;
    }
    
    if(i < no_infections){
      tmp2 = ti_pars[i];
      if(tmp2 <= t_i) y0 = y0;
      else if (tmp2 <= t_i + tp) y0 = (mu/tp)*tmp2 - (mu/tp)*t_i + y0;
      else y0 = m*(t_i + tp - tmp2) + mu + y0;
    }
   
    if(y0 < lower_titre_bound){
      y0 = lower_titre_bound;
    }
  }
  //out(_,0) = times;
  //out(_,1)=Y;
  return Y;
}

//[[Rcpp::export]]
NumericMatrix multiple_sim(
			   NumericVector ti_pars, 
			   NumericVector y0s, 
			   NumericMatrix mu_pars, 
			   NumericMatrix tp_pars, 
			   NumericMatrix m_pars, 
			   NumericVector times
			   ){
  double lower_bound = 0;
  double y0;
  NumericMatrix y(times.size(), ti_pars.size()+1);
  y(_,0) = times;
  NumericVector new_ti = clone(ti_pars);

  IntegerVector indices(ti_pars.size());
  NumericVector mus(ti_pars.size());
  NumericVector tps(ti_pars.size());
  NumericVector ms(ti_pars.size());

  new_ti = vector_sort(ti_pars);
  indices = vector_order(ti_pars) - 1;
  y0s = y0s[indices];

  for(int i = 0; i < ti_pars.size(); ++i){
    mus = mu_pars(i,_);
    mus = mus[indices];
    tps = tp_pars(i,_);
    tps = tps[indices];
    ms = m_pars(i,_);
    ms = ms[indices];
    y0 = y0s[i];
    y(_,i+1) = individual_sim(mus,tps,ms, new_ti,y0, lower_bound,times);
  }
  return(y);
}


//[[Rcpp::export]]
double obs_error(int actual, int obs, double S, double EA){
  int MAX_TITRE = 11;
  //  double S = 0.95;
  //  double EA = 0.02;
  if(actual == (MAX_TITRE) && obs == (MAX_TITRE)) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==0 && obs==0) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==obs) return(S);
  else if(actual == (obs + 1) || actual==(obs-1)) return(EA/2.0);
  return((1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
}


// Version 1 - only mu is mixed effects, rest is fixed
//[[Rcpp::export]]
double posterior_1(
			   NumericVector ti_pars, 
			   NumericVector y0s, 
			   double mu,
			   NumericMatrix cr_pars, 
			   NumericMatrix tp_pars, 
			   NumericMatrix m_pars, 
			   NumericMatrix data,
			   NumericVector error_pars,
			   double mu_pop,
			   double mu_pop_sigma
		 ){
  double pop_sigma, S, EA;
  pop_sigma = error_pars[0];
  S = error_pars[1];
  EA = error_pars[2];

  double lower_bound = 0;
  double y0;
  double lik = 0;
  NumericMatrix y(data.nrow(), ti_pars.size()+1);
  y(_,0) = data(_,0);
  NumericVector new_ti = clone(ti_pars);

  IntegerVector indices(ti_pars.size());
  NumericVector mus(ti_pars.size());
  NumericVector tps(ti_pars.size());
  NumericVector ms(ti_pars.size());

  new_ti = vector_sort(ti_pars);
  indices = vector_order(ti_pars) - 1;

  y0s = y0s[indices];

  for(int i = 0; i < ti_pars.size(); ++i){
    mus = cr_pars(i,_)*mu;
    mus = mus[indices];
    tps = tp_pars(i,_);
    tps = tps[indices];
    ms = m_pars(i,_);
    ms = ms[indices];
    y0 = y0s[i];
    y(_,i+1) = individual_sim(mus,tps,ms, new_ti,y0, lower_bound,data(_,0));
    
    /* ================================
       HERE TO CHANGE ERROR FUNCTION
       ================================ */
    //lik += likelihood_titre(y(_,i+1),data(_,i+1),S,EA);
    for(int j = 0; j < y.nrow();++j){
      lik += log(obs_error(data(j,i+1),floor(y(j,i+1)),S,EA));
      //lik += R::dnorm(data(j,i+1),y(j,i+1),pop_sigma,true);
    }
  }

  // ALL MIXED EFFECTS
  lik += R::dnorm(mu, mu_pop,mu_pop_sigma,true);
  //lik += R::dpois(floor(mu),floor(mu_pop),true);
  
  return(lik);
}



// Version 2 - mixed effects CR
//[[Rcpp::export]]
double posterior_2(
			   NumericVector ti_pars, 
			   NumericVector y0s, 
			   double mu,
			   NumericMatrix cr_pars, 
			   NumericMatrix tp_pars, 
			   NumericMatrix m_pars, 
			   NumericMatrix data,
			   double pop_sigma,
			   double mu_pop,
			   double mu_pop_sigma,
			   NumericVector cr_pop,
			   NumericVector cr_sigmas
		   ){
  double lower_bound = 0;
  double y0;
  double lik = 0;
  NumericMatrix y(data.nrow(), ti_pars.size()+1);
  y(_,0) = data(_,0);
  NumericVector new_ti = clone(ti_pars);

  IntegerVector indices(ti_pars.size());
  NumericVector mus(ti_pars.size());
  NumericVector tps(ti_pars.size());
  NumericVector ms(ti_pars.size());

  new_ti = vector_sort(ti_pars);
  indices = vector_order(ti_pars) - 1;

  y0s = y0s[indices];

  for(int i = 0; i < ti_pars.size(); ++i){
    mus = cr_pars(i,_)*mu;
    mus = mus[indices];
    tps = tp_pars(i,_);
    tps = tps[indices];
    ms = m_pars(i,_);
    ms = ms[indices];
    y0 = y0s[i];
    y(_,i+1) = individual_sim(mus,tps,ms, new_ti,y0, lower_bound,data(_,0));
    for(int j = 0; j < y.nrow();++j){
      lik += R::dnorm(data(j,i+1),y(j,i+1),pop_sigma,true);
    }
  }

  // MIXED EFFECTS
  lik += R::dnorm(mu, mu_pop,mu_pop_sigma,true);
 
  lik += R::dnorm(cr_pars(0,1),cr_pop[0],cr_sigmas[0],true);
  lik += R::dnorm(cr_pars(0,2),cr_pop[1],cr_sigmas[1],true);
  lik += R::dnorm(cr_pars(1,2),cr_pop[2],cr_sigmas[2],true);

  return(lik);
}


// Version 3 - mixed effects CR and m
//[[Rcpp::export]]
double posterior_3(
			   NumericVector ti_pars, 
			   NumericVector y0s, 
			   double mu,
			   NumericMatrix cr_pars, 
			   NumericMatrix tp_pars, 
			   NumericMatrix m_pars, 
			   NumericMatrix data,
			   double pop_sigma,
			   double mu_pop,
			   double mu_pop_sigma,
			   NumericVector cr_pop,
			   NumericVector m_pop,
			   NumericVector cr_sigmas,
			   NumericVector m_sigmas
		 ){
  double lower_bound = 0;
  double y0;
  double lik = 0;
  NumericMatrix y(data.nrow(), ti_pars.size()+1);
  y(_,0) = data(_,0);
  NumericVector new_ti = clone(ti_pars);

  IntegerVector indices(ti_pars.size());
  NumericVector mus(ti_pars.size());
  NumericVector tps(ti_pars.size());
  NumericVector ms(ti_pars.size());

  new_ti = vector_sort(ti_pars);
  indices = vector_order(ti_pars) - 1;

  y0s = y0s[indices];


  for(int i = 0; i < ti_pars.size(); ++i){
    mus = cr_pars(i,_)*mu;
    mus = mus[indices];
    tps = tp_pars(i,_);
    tps = tps[indices];
    ms = m_pars(i,_);
    ms = ms[indices];
    y0 = y0s[i];
    y(_,i+1) = individual_sim(mus,tps,ms, new_ti,y0, lower_bound,data(_,0));

    //lik += likelihood_titre(y(_,i+1),data(_,i+1),S,EA);
    for(int j = 0; j < y.nrow();++j){
      lik += R::dnorm(data(j,i+1),y(j,i+1),pop_sigma,true);
    }
    
  }
  lik += R::dnorm(mu, mu_pop,mu_pop_sigma,true);
  
  lik += R::dnorm(cr_pars(0,1),cr_pop[0],cr_sigmas[0],true);
  lik += R::dnorm(cr_pars(0,2),cr_pop[1],cr_sigmas[1],true);
  lik += R::dnorm(cr_pars(1,2),cr_pop[2],cr_sigmas[2],true); 
  
  lik += R::dnorm(m_pars(0,0),m_pop[0],m_sigmas[0],true);
  lik += R::dnorm(m_pars(0,1),m_pop[1],m_sigmas[1],true);
  lik += R::dnorm(m_pars(0,2),m_pop[2],m_sigmas[2],true);
  lik += R::dnorm(m_pars(1,2),m_pop[3],m_sigmas[3],true);

  return(lik);
}


// Version 1 - only mu is mixed effects, rest is fixed
//[[Rcpp::export]]
double posterior_4(
			   NumericVector ti_pars, 
			   NumericVector y0s, 
			   double mu,
			   NumericMatrix cr_pars, 
			   NumericMatrix tp_pars, 
			   NumericMatrix m_pars, 
			   NumericMatrix data,
			   NumericVector error_pars,
			   double mu_pop,
			   double mu_pop_sigma
		 ){
  double pop_sigma, S, EA;
  pop_sigma = error_pars[0];
  S = error_pars[1];
  EA = error_pars[2];

  double lower_bound = 0;
  double y0;
  double lik = 0;
  NumericMatrix y(data.nrow(), ti_pars.size()+1);
  y(_,0) = data(_,0);
  NumericVector new_ti = clone(ti_pars);

  IntegerVector indices(ti_pars.size());
  NumericVector mus(ti_pars.size());
  NumericVector tps(ti_pars.size());
  NumericVector ms(ti_pars.size());

  new_ti = vector_sort(ti_pars);
  indices = vector_order(ti_pars) - 1;

  y0s = y0s[indices];

  for(int i = 0; i < ti_pars.size(); ++i){
    mus = cr_pars(i,_)*mu;
    mus = mus[indices];
    tps = tp_pars(i,_);
    tps = tps[indices];
    ms = m_pars(i,_);
    ms = ms[indices];
    y0 = y0s[i];
    y(_,i+1) = individual_sim(mus,tps,ms, new_ti,y0, lower_bound,data(_,0));
    
    /* ================================
       HERE TO CHANGE ERROR FUNCTION
       ================================ */
    for(int j = 0; j < y.nrow();++j){
      lik += log(obs_error(data(j,i+1),floor(y(j,i+1)),S,EA));
      //lik += R::dnorm(data(j,i+1),y(j,i+1),pop_sigma,true);
    }
  }

  // ALL MIXED EFFECTS
  lik += R::dnorm(mu, mu_pop,mu_pop_sigma,true);
  
  for(int i = 0; i < (y.ncol()-1);++i){
    lik += log(obs_error(data(0,i+1),floor(y0s[i]),S,EA));
  }
  
  return(lik);
}
