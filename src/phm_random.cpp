// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <mvnorm.h>





using namespace Rcpp;

//------------------------------------ Class Specific Functions -------------------------------------------//
class phm_random{

public:

  // data
  Rcpp::List curr_tables; 
  Rcpp::List hist_tables;
  arma::vec n_intervals;
  int P; // number of columns of X
  int D; // total number of parameters
  
  

  // model definition

 
  // priors
  Rcpp::List prior_beta_mvn;
  arma::vec prior_lambda_hp1;
  arma::vec prior_lambda_hp2;
  // arma::vec prior_lambdah_hp1; 
  // arma::vec prior_lambdah_hp2;
  // arma::vec prior_a0_shape1;
  // arma::vec prior_a0_shape2;
  
  
  arma::vec                lower_limits;
  arma::vec                upper_limits;
  arma::vec                slice_widths;
  int m;


  // public member functions;
  phm_random(Rcpp::List & curr_tables0, Rcpp::List & hist_tables0,  
      arma::vec & n_intervals0, int P0,
      Rcpp::List prior_beta_mvn0,
      arma::vec prior_lambda_hp10, arma::vec prior_lambda_hp20,
      arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0);


  double logFC(arma::vec & parm0, const int & p);

};

phm_random::phm_random( Rcpp::List & curr_tables0, Rcpp::List & hist_tables0, 
          arma::vec & n_intervals0, int P0,
          Rcpp::List prior_beta_mvn0,
          arma::vec prior_lambda_hp10, arma::vec prior_lambda_hp20,
          arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0)
{

  curr_tables = curr_tables0;
  hist_tables = hist_tables0;

  
  n_intervals = n_intervals0;
  P = P0;
  
  D = P0 + sum(n_intervals0);
  

  prior_beta_mvn = prior_beta_mvn0;
  prior_lambda_hp1 = prior_lambda_hp10;
  prior_lambda_hp2 = prior_lambda_hp20;
  //prior_lambdah_hp1 = prior_lambdah_hp10;
  //prior_lambdah_hp2 = prior_lambdah_hp20;
  //prior_a0_shape1 = prior_a0_shape10;
  //prior_a0_shape2 = prior_a0_shape20;
  


  lower_limits = lower_limits0;
  upper_limits = upper_limits0;
  slice_widths = slice_widths0;

  m=10;
}

// Define the log likelihood 
double phm_random::logFC(arma::vec & parm0, const int & p)
{
  //Rcout << "parm0: " << parm0 << "\n";
  // extract regression parameters;
  int n_int = sum(n_intervals);
  arma::vec beta = parm0.subvec(0,P-1);
  arma::vec lambda = parm0.subvec(P,P+n_int-1);
  //arma::vec lambda_h = parm0.subvec(P+n_int,P+2*n_int-1);
  //arma::vec a0_vec = parm0.subvec(P+2*n_int, D-1);
  
  // compute likelihood contribution;
  double ll = 0;

  // current data likelihood;
  int Sn = curr_tables.size(); // number of strata
  for(int s = 0; s < Sn; s++){
    arma::mat tb = curr_tables[s];
    //Rcout << "tb: " << tb << "\n";
    arma::vec cumu = cumsum(n_intervals);
    int J = n_intervals[s];
    
    arma::vec lmd;
    if (s==0){
      lmd = lambda.subvec(0,J-1);
    }else{
      lmd = lambda.subvec(cumu[s-1],cumu[s-1]+J-1);
    }
    //Rcout << "lmd: " << lmd << "\n";
    for(int j = 0; j < J; j++){
      arma::mat tb_j = tb.rows(find(tb.col(P)==(j+1)));
      if(tb_j.n_rows != 0){
        arma::mat X = tb_j.cols(0,P-1);
        arma::mat stats = tb_j.cols(arma::span(tb_j.n_cols-2, tb_j.n_cols-1));
        double sum_ev = sum(stats.col(1));
        ll +=  sum_ev  * log(lmd[j])
              - lmd[j] * sum(exp(X * beta) % stats.col(0))
              +  sum((X * beta) % stats.col(1)); 
      }
    }
  }
    //Rcout << "ll: " << ll << "\n";


  // Rcpp::List nu_list(Sn); 
  // Rcpp::List r_list(Sn); 
  // for(int s = 0; s < Sn; s++){
  //   nu_list[s] = rep(0,n_intervals[s]);
  //   r_list[s] = rep(0,n_intervals[s]);
  // }
  // // compute p_sk and q_sk 
  // for(int j = 0; j < hist_tables.size(); j++){
  //   Rcpp::List histdata = hist_tables[j];
  // 
  //   for(int s = 0; s < Sn; s++){
  //     arma::mat tb = histdata[s];
  //     double a0 = a0_vec[j];
  //     
  //     int K = n_intervals[s];
  //     NumericVector v = as<NumericVector>(nu_list[s]);
  //     NumericVector w = as<NumericVector>(r_list[s]);
  //     
  //     for(int k = 0; k < K; k++){
  //       
  //       arma::mat tb_k = tb.rows(find(tb.col(P)==(k+1)));
  //       arma::mat X = tb_k.cols(0,P-1);
  //       arma::mat stats = tb_k.cols(arma::span(tb_k.n_cols-2, tb_k.n_cols-1));
  //       
  //       v[k] += a0*sum(stats.col(1));
  //       w[k] += a0*sum(exp(X * beta) % stats.col(0));
  //       
  //       ll += a0 * sum((X * beta) % stats.col(1));
  //     }
  //     nu_list[s] = v;
  //     r_list[s] = w;
  //   }
  // }
  // 
  // // flatten out nu_list and r_list
  // arma::vec cumu = cumsum(n_intervals);
  // arma::vec nu_list_flat(sum(n_intervals));
  // arma::vec r_list_flat(sum(n_intervals));
  // 
  // for(int s = 0; s < Sn; s++){
  //   int J = n_intervals[s];
  //   if (s==0){
  //     nu_list_flat.subvec(0,J-1) = as<arma::vec>(nu_list[s]);
  //     r_list_flat.subvec(0,J-1) = as<arma::vec>(r_list[s]);
  //   }else{
  //     nu_list_flat.subvec(cumu[s-1],cumu[s-1]+J-1) = as<arma::vec>(nu_list[s]);
  //     r_list_flat.subvec(cumu[s-1],cumu[s-1]+J-1) = as<arma::vec>(r_list[s]);
  //   }
  // }
  
  arma::mat parm0_mat(1, beta.size());
  parm0_mat.row(0) = beta.t();
  
    
  if(p < P){

    // MVN prior for beta
    for(int i = 0; i < prior_beta_mvn.size(); i++){
      Rcpp::List mvn = prior_beta_mvn[i];
      arma::vec mvn_mean = mvn[0];
      arma::mat mvn_cov = mvn[1];
      double w = mvn[2];
      //Rcout << "w: " << w << "\n";
      ll += w * dmvnorm(parm0_mat, mvn_mean, mvn_cov, TRUE)[0];
    }
  }else{
    
    if(p < (P+n_int)){
      
      // gamma prior for lambda
      ll +=  R::dgamma(parm0[p], prior_lambda_hp1[p-P], 1/(prior_lambda_hp2[p-P]), TRUE);
    }
    
  }

  //Rcout << "post: " << ll << "\n";
  return  ll;
}


// slice sampler
void slice( arma::vec & parms, phm_random & b)
{

  double b0, f0, f0_L, f0_R, f0_x1, h0, L, R, V, J, K,w,lower,upper;
  arma::vec parm0;

  
  for (int p = 0; p < (b.D); p++)
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
Rcpp::List phm_random_a0(Rcpp::List & curr_tables0, Rcpp::List & hist_tables0, 
                       arma::vec & n_intervals0, int P0, 
                       Rcpp::List & prior_beta_mvn0,
                       arma::vec prior_lambda_hp10, arma::vec prior_lambda_hp20, 
                       arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0, 
                       int nMC, int nBI){
  
   // declare object and set values;
   phm_random b(curr_tables0,hist_tables0,n_intervals0,P0,prior_beta_mvn0,
          prior_lambda_hp10,prior_lambda_hp20,
          lower_limits0,upper_limits0,slice_widths0);

    
    // Construct container for mcmc samples;
    arma::mat samples(nMC,b.D);

    // create parameter vector container and initial values;
    arma::vec parms(b.D);
    for (int p=0; p<b.D; p++)
    {

      parms[p]= R::runif(0,1);

    }

    for (int s=-nBI;s<nMC;s++)
    {
      slice(parms,b);

      if (s>=0){	samples.row(s) = parms.t();	}
    }
    
    arma::mat beta_samps = samples.cols(0,P0-1);
    arma::mat lambda_samps;
    //arma::mat lambda_h_samps;
    //arma::mat a0_samps;
    lambda_samps = samples.cols(P0,P0+sum(n_intervals0)-1);
    //lambda_h_samps = samples.cols(P0+sum(n_intervals0), P0+2*sum(n_intervals0)-1);
    //a0_samps = samples.cols(P0+2*sum(n_intervals0), b.D-1);
    
    Rcpp::List lams(n_intervals0.size());
    arma::vec cumu = cumsum(n_intervals0);
    
    for(int s = 0; s < n_intervals0.size(); s++){
      int J = n_intervals0[s];
      if (s==0){
        lams[s] = lambda_samps.cols(0,J-1);
      }else{
        lams[s] = lambda_samps.cols(cumu[s-1],cumu[s-1]+J-1);
      }
    }
    
    // Rcpp::List lams_h(n_intervals0.size());
    //   
    // for(int s = 0; s < n_intervals0.size(); s++){                            
    //   int J = n_intervals0[s];
    //   if (s==0){
    //     lams_h[s] = lambda_h_samps.cols(0,J-1);
    //   }else{
    //     lams_h[s] = lambda_h_samps.cols(cumu[s-1],cumu[s-1]+J-1);
    //   }
    // }
    
    return(Rcpp::List::create(
          Rcpp::Named("beta_samples") = beta_samps,
          Rcpp::Named("lambda_samples") = lams));


  
}


