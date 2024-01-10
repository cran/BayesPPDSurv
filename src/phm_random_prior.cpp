// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>





using namespace Rcpp;

//------------------------------------ Class Specific Functions -------------------------------------------//
class pch{

public:

  // data
  Rcpp::List hist_tables;
  arma::vec a0_vec;
  arma::vec n_intervals;
  int P; // number of columns of X
  
  

  // model definition

  // priors
  arma::vec prior_beta_mu;
  arma::vec prior_beta_sigma;
  arma::vec prior_lambdah_hp1; // baseline hazard parameters are not shared
  arma::vec prior_lambdah_hp2;
  
  
  arma::vec                lower_limits;
  arma::vec                upper_limits;
  arma::vec                slice_widths;
  int m;


  // public member functions;
  pch(Rcpp::List & hist_tables0, arma::vec a0_vec0,
      arma::vec & n_intervals0, int P0,
      arma::vec prior_beta_mu0, arma::vec prior_beta_sigma0,
      arma::vec prior_lambdah_hp10, arma::vec prior_lambdah_hp20,
      arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0);


  double logFC(arma::vec & parm0, const int & p);

};

pch::pch( Rcpp::List & hist_tables0, arma::vec a0_vec0,
          arma::vec & n_intervals0, int P0,
          arma::vec prior_beta_mu0, arma::vec prior_beta_sigma0, 
          arma::vec prior_lambdah_hp10, arma::vec prior_lambdah_hp20,
         arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0)
{

  hist_tables = hist_tables0;
  a0_vec = a0_vec0;
  
  n_intervals = n_intervals0;
  P = P0;

  
  prior_beta_mu = prior_beta_mu0;
  prior_beta_sigma = prior_beta_sigma0;
  prior_lambdah_hp1 = prior_lambdah_hp10;
  prior_lambdah_hp2 = prior_lambdah_hp20;


  lower_limits = lower_limits0;
  upper_limits = upper_limits0;
  slice_widths = slice_widths0;

  m=10;
}

// Define the log likelihood 
double pch::logFC(arma::vec & parm0, const int & p)
{

  // extract regression parameters;
  arma::vec beta = parm0;

  // compute likelihood contribution;
  double ll = 0;
  

    
  int Sn = n_intervals.size();// number of strata
  Rcpp::List nu_list(Sn); 
  Rcpp::List r_list(Sn); 
  for(int s = 0; s < Sn; s++){
    nu_list[s] = rep(0,n_intervals[s]);
    r_list[s] = rep(0,n_intervals[s]);
  }
  
  for(int j = 0; j < hist_tables.size(); j++){
    Rcpp::List histdata = hist_tables[j];
    double a0 = a0_vec[j];
    for(int s = 0; s < Sn; s++){
        arma::mat tb = histdata[s];
        
        int K = n_intervals[s];
        NumericVector v = as<NumericVector>(nu_list[s]);
        NumericVector w = as<NumericVector>(r_list[s]);
        
        for(int k = 0; k < K; k++){
          
          arma::mat tb_k = tb.rows(find(tb.col(P)==(k+1)));
          arma::mat X = tb_k.cols(0,P-1);
          arma::mat stats = tb_k.cols(arma::span(tb_k.n_cols-2, tb_k.n_cols-1));
          
          v[k] += a0*sum(stats.col(1));
          w[k] += a0*sum(exp(X * beta) % stats.col(0));
          
          ll += a0 * sum((X * beta) % stats.col(1));
        }
        nu_list[s] = v;
        r_list[s] = w;
    }
  }
  
  // create list for hyperparameters of lambda
  Rcpp::List c_list(n_intervals.size());
  Rcpp::List d_list(n_intervals.size());
  arma::vec cumu = cumsum(n_intervals);
  
  for(int s = 0; s < Sn; s++){
    int J = n_intervals[s];
    if (s==0){
      c_list[s] = prior_lambdah_hp1.subvec(0,J-1);
      d_list[s] = prior_lambdah_hp2.subvec(0,J-1);
    }else{
      c_list[s] = prior_lambdah_hp1.subvec(cumu[s-1],cumu[s-1]+J-1);
      d_list[s] = prior_lambdah_hp2.subvec(cumu[s-1],cumu[s-1]+J-1);
    }
  }
  
  
  for(int s = 0; s < Sn; s++){
    
    int K = n_intervals[s];
    NumericVector nu = nu_list[s];
    NumericVector r = r_list[s];
    NumericVector c = c_list[s];
    NumericVector d = d_list[s];
    
    for(int k = 0; k < K; k++){
      ll += - (nu[k] + c[k]) * log(r[k] + d[k]);
    }
    
   }
  
   // Add indepedent normal priors 
   ll += R::dnorm(parm0[p], prior_beta_mu[p], prior_beta_sigma[p], TRUE);


   return  ll;
}


// slice sampler
void slice( arma::vec & parms, pch & b)
{

  double b0, f0, f0_L, f0_R, f0_x1, h0, L, R, V, J, K,w,lower,upper;
  arma::vec parm0;

  
  for (int p = 0; p < (b.P); p++)
  {
    //Rcout << "p: " << p << "\n";
    //Rcout << "parm " << parm0 << "\n";
    // create vector of parameters to modify for slice sampling;
    parm0 = parms;

    // extract slice width and parameter bounds;
    w     = b.slice_widths[p];
    lower = b.lower_limits[p];
    upper = b.upper_limits[p];

    // skip over fixed parameter values;
    if (lower==upper){parms(p) = lower;}
    else
    {
      // current value of the parameter in question;
      b0 = parm0(p);

      // calculate current full conditional value;
      f0 = b.logFC(parm0,p);

      // calculate height of the horizontal slice;
      h0 = f0 - R::rexp(1.0);

      // Calculate initial horizontal interval;
      L = parm0(p) - R::runif(0.0,1.0)*w;
      R = L+w;

      // Truncate bounds to support of the parameter space;
      L = std::max(L,lower);
      R = std::min(R,upper);

      // Step out;
      V = R::runif(0.0,1.0);
      J = floor(b.m*V);
      K = (b.m-1)-J;

      // compute log of full conditional at current boundaries;
      parm0(p) = L; f0_L = b.logFC(parm0,p);
      parm0(p) = R; f0_R = b.logFC(parm0,p);

      while(J>0 and h0<f0_L and L>=lower)
      {
        L        = L-w; if (L<=lower) {L=lower;}
        J        = J-1;
        parm0(p) = L;
        f0_L     = b.logFC(parm0,p);
      }
      while(K>0 and h0<f0_R and R<=upper)
      {
        R        = R+w; if (R>=upper) {R=upper;}
        K        = K-1;
        parm0(p) = R;
        f0_R     = b.logFC(parm0,p);
      }

      // perform rejection sampling;
      int stop  = 0;
      while(stop == 0)
      {
        parm0(p)     = L + R::runif(0.0,1.0)*(R-L);
        f0_x1        = b.logFC(parm0,p);

        if      ( h0       <  f0_x1 ) { parms(p) = parm0(p); stop = 1;  }
        else if ( parm0(p) <  b0    ) { L = parm0(p);                     }
        else if ( parm0(p) >= b0    ) { R = parm0(p);                     }

        if (-0.0000000001 <= L-R and L-R <= 0.0000000001)
        {
          parms(p)= 0.5*(L+R);
          stop      = 1;
        }
      }
    }
  }
}






// [[Rcpp::export]]
Rcpp::List npp_beta(int L, Rcpp::List & hist_tables0, 
                       arma::vec & n_intervals0, int P0, 
                       arma::vec prior_a0_shape1, arma::vec prior_a0_shape2,
                       arma::vec prior_beta_mu0, arma::vec prior_beta_sigma0,
                       arma::vec prior_lambdah_hp10, arma::vec prior_lambdah_hp20,
                       arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0, 
                       int nBI){
  
  // declare object and set values;
  
  arma::mat beta_matrix(L, P0);
  
  for(int l = 0; l < L; l++){
    
    arma::vec a0_vec0(hist_tables0.size());
    
    for(int j=0; j < hist_tables0.size(); j++){
      a0_vec0[j] = R::rbeta(prior_a0_shape1[j], prior_a0_shape2[j]);
    }
    
    pch b(hist_tables0,a0_vec0,n_intervals0,P0,prior_beta_mu0,prior_beta_sigma0,
          prior_lambdah_hp10, prior_lambdah_hp20,
          lower_limits0,upper_limits0,slice_widths0);

    
    // Construct container for mcmc samples;
    //arma::mat samples(nMC,b.D);

    // create parameter vector container and initial values;
    arma::vec parms(P0);
    for (int p=0; p<P0; p++)
    {

      parms[p]= R::runif(0,1);

    }

    for (int s=-nBI;s<1;s++)
    {
      slice(parms,b);

      if (s>=0){	beta_matrix.row(l) = parms.t();	}
    }
    
  }
  return(Rcpp::List::create(
      Rcpp::Named("beta_matrix") = beta_matrix));

  
}


