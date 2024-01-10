
utils::globalVariables(c("r", "ev", "interval", "rtime", "evtime"))

#' Approximating the normalized power prior for \eqn{\beta} for the proportional hazards model with piecewise constant hazard and random a0
#'
#' @description Approximation of the normalized power prior for \eqn{\beta} for the proportional hazards model with piecewise constant hazard and random \eqn{a_0}.
#' The function returns discrete samples of \eqn{\beta} from the normalized power prior, and the user can use any mixture of multivariate normal distributions as an 
#' approximation for the normalized power prior for \eqn{\beta}.
#' This function is used to produce \code{prior.beta.mvn} in the function \code{\link{power.phm.random.a0}}.
#' 
#' @param prior.a0.shape1 Vector of the first shape parameters of the independent beta priors for \eqn{a_0}. The length of the vector should be equal to the number of historical datasets. The default is a vector of one's.
#' @param prior.a0.shape2 Vector of the second shape parameters of the independent beta priors for \eqn{a_0}. The length of the vector should be equal to the number of historical datasets. The default is a vector of one's.
#' @param prior.beta.mean Vector of means of the normal initial prior on \eqn{\beta}. The default value is zero for all the elements of \eqn{\beta}.
#' @param prior.beta.sd Vector of standard deviations of the normal initial prior on \eqn{\beta}. The default value is 10^3 for all the elements of \eqn{\beta}. 
#' @param prior.lambda0.hp1 Vector of first hyperparameters of the Gamma initial prior on \eqn{\lambda_0}. 
#' The length of the vector should be equal to the dimension of \eqn{\lambda_0}, i.e., the total number of intervals for all strata. The default value is 10^(-5) for all the elements of \eqn{\lambda_0}.
#' @param prior.lambda0.hp2 Vector of second hyperparameters of the Gamma initial prior on \eqn{\lambda_0}. 
#' The length of the vector should be equal to the dimension of \eqn{\lambda_0}, i.e., the total number of intervals for all strata. The default value is 10^(-5) for all the elements of \eqn{\lambda_0}.
#' @param lower.limits Vector of lower limits for \eqn{\beta} to be used by the slice sampler. The length of the vector should be equal to the length of \eqn{\beta}. The default is -100 for all the elements of \eqn{\beta} (may not be appropriate for all situations). 
#' @param upper.limits Vector of upper limits for \eqn{\beta} to be used by the slice sampler. The length of the vector should be equal to the length of \eqn{\beta}. The default is 100 for all the elements of \eqn{\beta} (may not be appropriate for all situations).
#' @param slice.widths Vector of initial slice widths for \eqn{\beta} to be used by the slice sampler. The length of the vector should be equal to the total number of parameters. The default is 0.1 for all the elements of \eqn{\beta} (may not be appropriate for all situations). 
#' 
#' @inheritParams power.phm.fixed.a0
#' @details 
#' 
#' This function is used to produce \code{prior.beta.mvn} in the function \code{\link{power.phm.random.a0}}. It approximates the normalized power prior for \eqn{\beta} when \eqn{a_0} is modeled as random.
#' The function returns discrete samples of \eqn{\beta} from the normalized power prior, and the user can use any mixture of multivariate normal distributions as an 
#' approximation for the normalized power prior for \eqn{\beta}.
#' 
#' Baseline hazard parameters for the 
#' current and historical data are NOT shared. 
#' The baseline hazards of the historical data are denoted by \eqn{\lambda_0}. We assume Gamma initial priors for \eqn{\lambda_0}
#' and independent normal initial priors for \eqn{\beta}.   
#' 
#' Posterior samples are obtained through slice sampling. 
#' The default lower limits are -100 for \eqn{\beta}. The default upper limits 
#' for the parameters are 100. The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#' 
#' @return Samples of \eqn{\beta} (approximating the normalized power prior) are returned.
#' @references Ibrahim, J. G., Chen, M.-H. and Sinha, D. (2001). Bayesian Survival Analysis. New York: Springer Science & Business Media.
#'
#' Psioda, M. A. and Ibrahim, J. G. (2019). Bayesian clinical trial design using historical data that inform the treatment effect. Biostatistics 20, 400–415.
#' 
#' Shen, Y., Psioda, M. A., and Joseph, J. G. (2023). BayesPPD: an R package for Bayesian sample size determination using the power and normalized power prior for generalized linear models. The R Journal, 14(4).
#' 
#' @seealso \code{\link{phm.random.a0}} and \code{\link{power.phm.random.a0}}
#' @examples 
#' 
#' # Simulate two historical datasets
#' n <- 100
#' P <- 4
#' time1 <- round(rexp(n, rate=0.5),1)
#' event1 <- rep(1,n)
#' X1 <- matrix(rbinom(n*P,prob=0.5,size=1), ncol=P)
#' S1 <- c(rep(1,n/2),rep(2,n/2))
#' time2 <- round(rexp(n, rate=0.7),1)
#' event2 <- rep(1,n)
#' X2 <- matrix(rbinom(n*P,prob=0.5,size=1), ncol=P)
#' S2 <- c(rep(1,n/2),rep(2,n/2))
#' historical <- list(list(time=time1, event=event1, X=X1, S=S1),
#'                    list(time=time2, event=event2, X=X2, S=S2))
#' 
#' # We choose three intervals for the first stratum and two intervals for the second stratum
#' n.intervals <- c(3,2) 
#' change.points <- list(c(1,2), 2)
#' 
#' # Get samples from the approximate normalized power prior for beta
#' nMC <- 100 # nMC should be larger in practice
#' nBI <- 50
#' prior.beta <- approximate.prior.beta(historical, n.intervals, change.points=change.points,
#'                                      prior.a0.shape1=c(1,1), prior.a0.shape2=c(1,1), 
#'                                      nMC=nMC, nBI=nBI)
#' prior_beta_mu=colMeans(prior.beta)
#' prior_beta_sigma=cov(prior.beta) 
#' 
#' # Aprroximate the discrete sames with a single multivariate normal with weight one.
#' # The user can use any mixture of multivariate normal distributions as an 
#' # approximation for the normalized power prior for beta.
#' prior.beta.mvn <- list(list(prior_beta_mu, prior_beta_sigma, 1))
#' # prior.beta.mvn is a parameter for phm.random.a0() and power.phm.random.a0()
#' 
#' 
#' 
#' @export
approximate.prior.beta <- function(historical, n.intervals, change.points=NULL,
                       prior.a0.shape1=rep(1,10), prior.a0.shape2=rep(1,10), 
                       prior.beta.mean=rep(0,50), prior.beta.sd=rep(1000,50),
                       prior.lambda0.hp1=rep(10^(-5),50), prior.lambda0.hp2=rep(10^(-5),50), 
                       lower.limits=rep(-100, 50), upper.limits=rep(100, 50), slice.widths=rep(0.1, 50),
                       nMC=10000, nBI=250){
  
  # add zero and infinity to change.points
  change.points.new <- list()
  if(is.null(change.points)){
      change.points.new <- create_intervals_historical(historical, n.intervals)
  }else{
    for(i in 1:length(change.points)){
      l <- change.points[[i]]
      l1 <- unique(c(0, l, Inf))
      change.points.new[[i]] <- l1
    }
  }
  tables <- collapse_data(historical=historical, n.intervals=n.intervals, change.points=change.points.new, dCurrent=FALSE)
  t0 <- tables[["hist_tables"]]
  P <- ncol(historical[[1]]$X)
  
  samples <- npp_beta(nMC, t0, n.intervals, P, 
                      prior.a0.shape1, prior.a0.shape2,
                      prior.beta.mean, prior.beta.sd,
                      prior.lambda0.hp1, prior.lambda0.hp2,
                      lower.limits, upper.limits, slice.widths, nBI)
  
  res <- samples$beta_matrix
  
  return(res)
  
}



#' Power/type I error calculation for the proportional hazards model with piecewise constant hazard and random a0
#'
#' @description Power/type I error calculation using the normalized power prior for the proportional hazards model with piecewise constant hazard and random \eqn{a_0}
#'
#' @param prior.beta.mvn List of vectors of multivariate normal approximations of the normalized power prior for \eqn{\beta}. Each vector has three elements, 
#' the mean vector, the covariance matrix and the weight of the multivariate normal distribution. The normalized power prior for \eqn{\beta} 
#' is approximated by the weighted mixture of the multivariate normal distributions provided. By default, a single multivariate normal distribution is assumed.
#' The user can use the \code{\link{approximate.prior.beta}} function to obtain samples of \eqn{\beta} from the normalized power prior, and use any mixture of multivariate normals to approximate 
#' the normalized power prior for \eqn{\beta}.
#' @param prior.lambda.hp1 Vector of first hyperparameters of the Gamma initial prior on \eqn{\lambda}. 
#' The length of the vector should be equal to the dimension of \eqn{\lambda}, i.e., the total number of intervals for all strata. The default value is 10^(-5) for all the elements of \eqn{\lambda}.
#' @param prior.lambda.hp2 Vector of second hyperparameters of the Gamma initial prior on \eqn{\lambda}. 
#' The length of the vector should be equal to the dimension of \eqn{\lambda}, i.e., the total number of intervals for all strata. The default value is 10^(-5) for all the elements of \eqn{\lambda}.
#' @param lower.limits Vector of lower limits for parameters (\eqn{\beta} and \eqn{\lambda}, in this order) to be used by the slice sampler. The length of the vector should be equal to the total number of parameters. The default is -100 for \eqn{\beta} and 0 for \eqn{\lambda} (may not be appropriate for all situations). 
#' @param upper.limits Vector of upper limits for parameters (\eqn{\beta} and \eqn{\lambda}, in this order) to be used by the slice sampler. The length of the vector should be equal to the total number of parameters. The default is 100 for all parameters (may not be appropriate for all situations).
#' @param slice.widths Vector of initial slice widths for parameters (\eqn{\beta} and \eqn{\lambda}, in this order) to be used by the slice sampler. The length of the vector should be equal to the total number of parameters. The default is 0.1 for all parameters (may not be appropriate for all situations). 
#' 
#' @inheritParams power.phm.fixed.a0
#' @inheritParams approximate.prior.beta
#' 
#' @details The proportional hazards model with piecewise constant hazard is implemented. 
#' We assume \eqn{\beta} is the regression coefficients. We assume the first column of the covariate matrix is the treatment indicator,
#' and the corresponding parameter is \eqn{\beta_1}. Here \eqn{a_0} is modeled as random with a normalized power prior. 
#' 
#' The normalized power prior for \eqn{\beta} is approximated by a weighted mixture of multivariate normal distributions provided in \code{prior.beta.mvn}. 
#' The user can use the \code{\link{approximate.prior.beta}} function to obtain samples of \eqn{\beta} from the normalized power prior, and use any mixture of multivariate normals to approximate 
#' the normalized power prior for \eqn{\beta}. By default, a single multivariate normal distribution is assumed.
#' 
#' Baseline hazard parameters for the 
#' current and historical data are NOT shared. The baseline hazards of the current data are denoted by \eqn{\lambda}. 
#' The baseline hazards of the historical data are denoted by \eqn{\lambda_0}. We assume Gamma initial priors for 
#' \eqn{\lambda} and \eqn{\lambda_0}.  
#' 
#' To perform sample size determination, we test the hypotheses
#' \deqn{H_0: \beta_1 \ge \delta} and \deqn{H_1: \beta_1 < \delta.}
#' 
#' The sampling prior for the treatment parameter can be generated from a normal distribution (see examples). 
#' For example, suppose one wants to compute the power for the hypotheses \eqn{H_0: \beta_1 \ge 0} and \eqn{H_1: \beta_1 < 0.} 
#' To approximate the sampling prior for \eqn{\beta_1}, one can simply sample from a normal distribution with negative mean, 
#' so that the mass of the prior falls in the alternative space. Conversely, to compute the type I error rate, one can 
#' sample from a normal distribution with positive mean, so that the mass of the prior falls in the null space.
#' 
#' The sampling prior for the other parameters (\eqn{\beta_2}, ..., \eqn{\beta_p} and \eqn{\lambda}) can be generated from the posterior based on the historical data. 
#' This can be achieved by the function \link{phm.fixed.a0}
#' with \code{current.data} set to \code{FALSE} (see the vignette).  
#' 
#' Posterior samples are obtained through slice sampling. 
#' The default lower limits are -100 for \eqn{\beta} and 0 for \eqn{\lambda}. The default upper limits 
#' for the parameters are 100. The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' If a sampling prior with support in the null space is used, the value returned is a Bayesian type I error rate.
#' If a sampling prior with support in the alternative space is used, the value returned is a Bayesian power.
#' 
#' @return Power or type I error is returned, depending on the sampling prior used.
#' The posterior probabilities of the alternative hypothesis are returned.  
#' The average posterior means of \eqn{\beta} and \eqn{\lambda} are also returned. 
#' @references Ibrahim, J. G., Chen, M.-H. and Sinha, D. (2001). Bayesian Survival Analysis. New York: Springer Science & Business Media.
#'
#' Psioda, M. A. and Ibrahim, J. G. (2019). Bayesian clinical trial design using historical data that inform the treatment effect. Biostatistics 20, 400–415.
#' 
#' Shen, Y., Psioda, M. A., and Joseph, J. G. (2023). BayesPPD: an R package for Bayesian sample size determination using the power and normalized power prior for generalized linear models. The R Journal, 14(4).
#' @seealso \code{\link{phm.random.a0}} and \code{\link{approximate.prior.beta}}
#' @examples
#' 
#' 
#' # Simulate two historical datasets
#' set.seed(1)
#' n <- 100
#' P <- 4
#' time1 <- round(rexp(n, rate=0.5),1)
#' event1 <- rep(1,n)
#' X1 <- matrix(rbinom(n*P,prob=0.5,size=1), ncol=P)
#' S1 <- c(rep(1,n/2),rep(2,n/2))
#' time2 <- round(rexp(n, rate=0.7),1)
#' event2 <- rep(1,n)
#' X2 <- matrix(rbinom(n*P,prob=0.5,size=1), ncol=P)
#' S2 <- c(rep(1,n/2),rep(2,n/2))
#' historical <- list(list(time=time1, event=event1, X=X1, S=S1),
#'                    list(time=time2, event=event2, X=X2, S=S2))
#' 
#' n.subjects <- 100
#' n.events <- 30
#' 
#' # We choose three intervals for the first stratum and two intervals for the second stratum
#' n.intervals <- c(3,2) 
#' 
#' # Generate sampling priors
#' 
#' # The null hypothesis here is H0: beta_1 >= 0. To calculate power,
#' # we can provide samples of beta_1 such that the mass of beta_1 < 0.
#' # To calculate type I error, we can provide samples of beta_1 such that
#' # the mass of beta_1 >= 0.
#' samp.prior.beta1 <- rnorm(100, mean=-1, sd=1)
#' # Here, mass is put on the alternative region, so power is calculated.
#' samp.prior.beta <- cbind(samp.prior.beta1, matrix(rnorm(100*(P-1)), 100, P-1))
#' 
#' # Point mass sampling priors are used for lambda
#' lambda_strat1 <- matrix(c(0.5, 0.5, 0.5), nrow=1)
#' lambda_strat2 <- matrix(c(0.7, 0.7), nrow=1)
#' samp.prior.lambda <- list(lambda_strat1, lambda_strat2)
#' 
#' 
#' nMC <- 50 # nMC should be larger in practice
#' nBI <- 50
#' N <- 5 # N should be larger in practice
#' 
#' result <- power.phm.random.a0(historical=historical, n.subjects=n.subjects, 
#'                               n.events=n.events, n.intervals=n.intervals, 
#'                               samp.prior.beta=samp.prior.beta, 
#'                               samp.prior.lambda=samp.prior.lambda,
#'                               prior.a0.shape1 = c(1,1), prior.a0.shape2 = c(1,1),
#'                               dist.enroll="Uniform", param.enroll=0.5,
#'                               nMC=nMC, nBI=nBI, delta=0, nullspace.ineq=">", N=N)
#' result$`power/type I error`
#' result$`average posterior mean of beta`
#' result$`average posterior mean of lambda`
#' 
#' 
#'
#' @export
#' @import dplyr tidyr
#' @importFrom stats cov runif rexp rbinom
power.phm.random.a0 <- function(historical, n.subjects, n.events, 
                               n.intervals, change.points=NULL, 
                               samp.prior.beta, samp.prior.lambda, # list of matrices
                               dist.enroll, param.enroll,
                               rand.prob=0.5, prob.drop=0, param.drop=0, 
                               dist.csr="Constant", param.csr=10000, min.follow.up=0, max.follow.up=10000,
                               prior.beta.mvn=NULL,
                               prior.a0.shape1=rep(1,10), prior.a0.shape2=rep(1,10), 
                               prior.lambda.hp1=rep(10^(-5),50), prior.lambda.hp2=rep(10^(-5),50), 
                               lower.limits=NULL, upper.limits=rep(100, 50), slice.widths=rep(0.1, 50),
                               nMC=10000, nBI=250, delta=0, nullspace.ineq=">", gamma=0.95, N=10000){
  
  # add zero and infinity to change.points
  change.points.new <- list()
  if(is.null(change.points)){
    change.points.new <- create_intervals_historical(historical,n.intervals)
  }else{
    for(i in 1:length(change.points)){
      l <- change.points[[i]]
      l1 <- unique(c(0, l, Inf))
      change.points.new[[i]] <- l1
    }
  }
  
  # if prior.beta.mvn is NULL, make its default value a single multivariate normal
  if(is.null(prior.beta.mvn)){
    prior.beta <- approximate.prior.beta(historical, n.intervals, change.points=change.points.new,
                                         prior.a0.shape1=prior.a0.shape1, prior.a0.shape2=prior.a0.shape2, 
                                         nMC=nMC, nBI=nBI)
    prior_beta_mu=colMeans(prior.beta)
    prior_beta_sigma=cov(prior.beta) # sigma is standard deviation
    prior.beta.mvn <- list(list(prior_beta_mu, prior_beta_sigma, 1))
  }


  # build matrix of covariates to sample from
  P <- ncol(historical[[1]]$X)
  x_sim <- matrix(NA, nrow=0, ncol=P-1)
  for(k in 1:length(historical)){
    dat <- historical[[k]]
    x_h <- as.matrix(dat$X[,-1])
    x_sim <- rbind(x_sim, x_h)
  }
  # build matrix of strata to sample from

  s_sim <- NULL
  for(k in 1:length(historical)){
    dat <- historical[[k]]
    s_h <- dat$S
    s_sim <- c(s_sim, s_h)
  }
  
  # repeat N times 
  
  posterior_probs <- NULL
  
  # save the posterior means of beta and lambda
  beta_mean <- matrix(NA, nrow=N, ncol=P)
  lambda_sum <- list()
  for(j in 1:length(n.intervals)){
    lambda_sum[[j]] <- rep(0,n.intervals[j])
  }

  for(l in 1:N){
    # print(l)
    # sample stratum variable
    ind <- sample(1:length(s_sim), n.subjects, replace = TRUE)
    s <- s_sim[ind]
    
    # sample beta
    ind <-  sample(1:nrow(samp.prior.beta), 1, replace = TRUE)
    beta <- samp.prior.beta[ind,]
    
    # sample lambda
    list_lambda <- list()
    for(j in 1:length(unique(s))){
      lams <- samp.prior.lambda[[j]]
      ind <-  sample(1:nrow(lams), 1, replace = TRUE)
      list_lambda <- c(list_lambda, list(lams[ind,]))
    }
    
    # sample enrollment times
    if(dist.enroll=="Uniform"){
      r <- runif(n.subjects,0,param.enroll)
    }else{
      r <- rexp(n.subjects,rate=param.enroll)
    }
    
    # sample covariates and treatment variable
    ind <-  sample(1:nrow(x_sim), n.subjects, replace = TRUE)
    x <- x_sim[ind,]
    z <- rbinom(n.subjects,size=1,rand.prob)
    x <- cbind(trt=z, x)
    
    
    phi <- exp(x%*%beta)
    
    # simulate event time ti
    eventtimes <- NULL

    for(i in 1:n.subjects){

      strat <- s[i]
      lambda_s <- list_lambda[[strat]]
      st <- change.points.new[[strat]]
      st_1 <- st[-length(st)]

      
      k <- 1
      
      while(TRUE){

        theta_k <- phi[i]*lambda_s[k]
        
        t_tilde <- rexp(1,rate=theta_k) + st_1[k]
        if(k > length(st_1) | t_tilde <= st[k+1]){
          break
        }
        k <- k + 1
      }
      eventtimes <- c(eventtimes, t_tilde)
    }
    
    # sample (administrative) censorship times
    if(dist.csr=="Uniform"){
      ctime <- runif(n.subjects,0,param.csr)
    }else if(dist.csr=="Constant"){
      ctime <- rep(param.csr, n.subjects)
    }else{
      ctime <- rexp(n.subjects,rate=param.csr)
    }
    y <- ifelse(eventtimes <= ctime, eventtimes, ctime)
    nu <- ifelse(eventtimes <= ctime, 1, 0)
    
    # simulate dropouts
    # simulate a dropout time; if dropout time < y, the person is a dropout; 
    # otherwise the person has an event or is administratively censored
    num_drops <- round(n.subjects * prob.drop, 0)
    
    if(num_drops > 0){
      drop_ind <- sample(1:n.subjects, size=num_drops, replace=FALSE)
      droptime <- runif(num_drops, 0, param.drop)
      obstime <- ifelse(y[drop_ind] <= droptime, y[drop_ind], droptime)
      bool_drop <- ifelse(y[drop_ind] <= droptime, 0, 1)
      nu[drop_ind][bool_drop==1] <- 0
      y[drop_ind] <- obstime
    }
    
    e <- r+y
    complete_data <- data.frame(X=x, S=s, enrtimes=r, eventtimes=eventtimes, ctime=ctime, y=y, nu=nu, t_elps=e)

    # create original data
    df_events <- complete_data[complete_data$nu==1,]
    df_events1 <- df_events[order(df_events$t_elps),]
    if(nrow(df_events1) < n.events){
      stoptime <- max.follow.up
    }else{
      stoptime <- df_events1[n.events, "t_elps"]
      if(stoptime < min.follow.up){
        stoptime <- min.follow.up
      }
    }

    finaldf <- complete_data[complete_data$enrtimes < stoptime,]
    finaldf$new_y <- ifelse(finaldf$t_elps > stoptime, stoptime-finaldf$enrtimes, finaldf$y)
    finaldf$new_nu <- ifelse(finaldf$t_elps > stoptime, 0, finaldf$nu)
    
    # create tables
    tables <- collapse_data(time=finaldf$new_y, event=finaldf$new_nu, X=finaldf[,1:ncol(x)], S=finaldf$S, 
                  historical=historical, n.intervals=n.intervals, change.points=change.points.new, dCurrent=TRUE)
    t1 <- tables[["curr_tables"]]
    t2 <- tables[["hist_tables"]]
    #print(t1)
    #print(t2)
    if(is.null(lower.limits)){
      lower.limits = c(rep(-100,P), rep(0,2*sum(n.intervals)))
    }
    samples <- phm_random_a0(t1, t2,
                             n.intervals, P, prior.beta.mvn, prior.lambda.hp1, prior.lambda.hp2, 
                             lower.limits, upper.limits, slice.widths, 
                             nMC, nBI)
    
    # compute probability of success
    beta1 <- samples$beta_samples[,1]
    
    if(nullspace.ineq == ">"){
      pr <- sum(beta1 < delta) / length(beta1)
    }else{
      pr <- sum(beta1 > delta) / length(beta1)
    }
    
    posterior_probs <- c(posterior_probs, pr)
    
    # average posterior mean of beta 
    beta_mean[l,] <- colMeans(samples$beta_samples)
    
    # sum of posterior mean of lambda
    list_lambda <- samples$lambda_samples
    list_lambda_mean <- list()
    for(j in 1:length(list_lambda)){
      list_lambda_mean[[j]] <- colMeans(list_lambda[[j]])
    }
    for(j in 1:length(list_lambda)){
      lambda_sum[[j]] <- lambda_sum[[j]]+list_lambda_mean[[j]]
    }

  }
  
  # average of posterior mean of lambda
  for(j in 1:length(n.intervals)){
    lambda_sum[[j]] <- lambda_sum[[j]]/N
  }
  
  power <- mean(posterior_probs >= gamma)
  
  
  return(list(#"simulated dataset"=simdf, 
      "power/type I error"=power, 
      "posterior probabilities"=posterior_probs,
      "average posterior mean of beta"=colMeans(beta_mean),
      "average posterior mean of lambda"=lambda_sum))


}




 
#' Model fitting for the proportional hazards model with piecewise constant hazard and random a0
#'
#' @description Model fitting using the normalized power prior for the proportional hazards model with piecewise constant hazard and random \eqn{a_0}
#' @param change.points List of vectors. Each vector in the list contains the change points for the baseline hazards for each stratum. The length of the vector should be equal to the total number of strata. 
#' By default, we assign the change points so that the same number of events are observed in all the intervals in the pooled current and historical data.  
#' 
#' @inheritParams power.phm.random.a0
#' @inheritParams phm.fixed.a0
#' 
#' @details The proportional hazards model with piecewise constant hazard is implemented. 
#' We assume \eqn{\beta} is the regression coefficients. We assume the first column of the covariate matrix is the treatment indicator,
#' and the corresponding parameter is \eqn{\beta_1}. Here \eqn{a_0} is modeled as random with a normalized power prior. 
#' 
#' The normalized power prior for \eqn{\beta} is approximated by a weighted mixture of multivariate normal distributions provided in \code{prior.beta.mvn}. 
#' The user can use the \code{\link{approximate.prior.beta}} function to obtain samples of \eqn{\beta} from the normalized power prior, and use any mixture of multivariate normals to approximate 
#' the normalized power prior for \eqn{\beta}. By default, a single multivariate normal distribution is assumed.
#' 
#' Posterior samples are obtained through slice sampling. 
#' The default lower limits are -100 for \eqn{\beta} and 0 for \eqn{\lambda}. The default upper limits 
#' for the parameters are 100. The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' 
#' @return Posterior samples of \eqn{\beta} and \eqn{\lambda} are returned.
#' @references Ibrahim, J. G., Chen, M.-H. and Sinha, D. (2001). Bayesian Survival Analysis. New York: Springer Science & Business Media.
#'
#' Psioda, M. A. and Ibrahim, J. G. (2019). Bayesian clinical trial design using historical data that inform the treatment effect. Biostatistics 20, 400–415.
#' 
#' Shen, Y., Psioda, M. A., and Joseph, J. G. (2023). BayesPPD: an R package for Bayesian sample size determination using the power and normalized power prior for generalized linear models. The R Journal, 14(4).
#' @seealso \code{\link{power.phm.random.a0}} and \code{\link{approximate.prior.beta}}
#' @examples 
#' 
#' 
#' 
#' set.seed(1)
#' # Simulate current data
#' n <- 50
#' P <- 4
#' time <- round(rexp(n, rate=0.5),1)
#' event <- rep(1,n)
#' X <- matrix(rbinom(n*P,prob=0.5,size=1), ncol=P)
#' S <- c(rep(1,n/2),rep(2,n/2))
#' 
#' # Simulate two historical datasets
#' n <- 100
#' time1 <- round(rexp(n, rate=0.5),1)
#' event1 <- rep(1,n)
#' X1 <- matrix(rbinom(n*P,prob=0.5,size=1), ncol=P)
#' S1 <- c(rep(1,n/2),rep(2,n/2))
#' time2 <- round(rexp(n, rate=0.7),1)
#' event2 <- rep(1,n)
#' X2 <- matrix(rbinom(n*P,prob=0.5,size=1), ncol=P)
#' S2 <- c(rep(1,n/2),rep(2,n/2))
#' historical <- list(list(time=time1, event=event1, X=X1, S=S1),
#'                    list(time=time2, event=event2, X=X2, S=S2))
#' 
#' # We choose three intervals for the first stratum and two intervals for the second stratum
#' n.intervals <- c(3,2) 
#' change.points <- list(c(1,2), 2)
#' 
#' # Get samples from the approximate normalized power prior for beta
#' nMC <- 100 # nMC should be larger in practice
#' nBI <- 50
#' prior.beta <- approximate.prior.beta(historical, n.intervals, change.points=change.points,
#'                                      prior.a0.shape1=c(1,1), prior.a0.shape2=c(1,1), 
#'                                      nMC=nMC, nBI=nBI)
#' prior_beta_mu=colMeans(prior.beta)
#' prior_beta_sigma=cov(prior.beta) 
#' 
#' # Aprroximate the discrete sames with a single multivariate normal with weight one
#' prior.beta.mvn <- list(list(prior_beta_mu, prior_beta_sigma, 1))
#' 
#' result <- phm.random.a0(time=time, event=event, X=X, S=S,
#'                         historical=historical, n.intervals=n.intervals, 
#'                         change.points=change.points, 
#'                         prior.beta.mvn=prior.beta.mvn,
#'                         nMC=nMC, nBI=nBI)
#' 
#' # posterior mean of beta
#' colMeans(result$beta_samples) 
#' # posterior mean of baseline hazards for stratum 1
#' colMeans(result$lambda_samples[[1]]) 
#' # posterior mean of baseline hazards for stratum 2
#' colMeans(result$lambda_samples[[2]]) 
#' 
#' 
#' 
#' @export
#' @import dplyr tidyr
#' @importFrom stats cov 
phm.random.a0 <- function(time, event, X, S, historical, n.intervals, change.points=NULL, prior.beta.mvn=NULL, 
                        prior.lambda.hp1=rep(10^(-5),50), prior.lambda.hp2=rep(10^(-5),50), 
                        prior.a0.shape1=rep(1,10), prior.a0.shape2=rep(1,10), 
                        lower.limits=NULL, upper.limits=rep(100, 50), slice.widths=rep(0.1, 50),
                        nMC=10000, nBI=250){
  
  # add zero and infinity to change.points
  change.points.new <- list()
  if(is.null(change.points)){
    change.points.new <- create_intervals(time, event, S, historical, n.intervals)
  }else{
    for(i in 1:length(change.points)){
      l <- change.points[[i]]
      l1 <- unique(c(0, l, Inf))
      change.points.new[[i]] <- l1
    }
  }
  # create tables
  tables <- collapse_data(time=time, event=event, X=X, S=S, 
                historical=historical, n.intervals=n.intervals, change.points=change.points.new, dCurrent=TRUE)
  t1 <- tables[["curr_tables"]]
  t2 <- tables[["hist_tables"]]
  
  P <- ncol(X)

  if(is.null(lower.limits)){
    lower.limits = c(rep(-100,P), rep(0,sum(n.intervals)))
  }
  
  # if prior.beta.mvn is NULL, make its default value a single multivariate normal
  if(is.null(prior.beta.mvn)){
    prior.beta <- approximate.prior.beta(historical, n.intervals, change.points=change.points.new,
                                         prior.a0.shape1=prior.a0.shape1, prior.a0.shape2=prior.a0.shape2, 
                                         nMC=nMC, nBI=nBI)
    prior_beta_mu=colMeans(prior.beta)
    prior_beta_sigma=cov(prior.beta) # sigma is standard deviation
    prior.beta.mvn <- list(list(prior_beta_mu, prior_beta_sigma, 1))
  }
  samples <- phm_random_a0(t1, t2,
                          n.intervals, P, prior.beta.mvn, prior.lambda.hp1, prior.lambda.hp2, 
                          lower.limits, upper.limits, slice.widths, 
                          nMC, nBI)
  
  return(samples)
}
    
