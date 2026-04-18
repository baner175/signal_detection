library(truncdist)
library(VGAM)
library(doSNOW)
library(parallel)
library(foreach)

#-------------------------------------------------------------------------------
## NEGATIVE LOG-LIKELIHOOD #####################################################
## with the benchmark model q_\alpha (when bkg data is not available) ##########
## or with proposal background g_\beta (when a bkg data is available) ##########

neg_loglik <- function(beta, data){
  q_i <- sapply(data, dtrunc, spec = 'pareto', a = l, b = u,
                scale = l, shape = beta)
  return(-sum(log(q_i)))
}

# Integrand in (||S||_{G_{\beta}})^2 as a function of x and \beta
normS_integrand <- function(x, beta){
  fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
               mean = mean_sig, sd = sd_sig)
  g <- dtrunc(x, spec = 'pareto', a = l, b = u,
              scale = l, shape = beta)
  return(((fs/g-1)^2)*g)
}

#-------------------------------------------------------------------------------
##################### FUNCTIONS FOR SIMULATIONS ################################

# power simulation function in the presence of a bkg data using Pareto type I proposal background:
sim_power_wbkg <- function(eta, n_phys,
                           r, nsims, seed = 12345,
                           signif.level = 0.05,
                           beta0 = NULL){
  ####################### DESCRIPTION OF ARGUMENTS #############################
  # eta: true value of eta for simulation
  # n_phys: Size of the physics data
  # r: bkg data to physics data raio
  # nsims: number of iterations
  # seed: seed for reproducibility
  # signif.level: significance level for the test
  # beta0: the known value of the parameter in g_\beta, estimated via MLE if NULL
  ##############################################################################
  
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f', eta))
  cat("\n--------------------------------------\n")
  
  B <- nsims; n_bkg <- r*n_phys
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  
  # parallelization for faster implementation
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  clusterExport(cl, ls(envir = environment()), envir = environment()) # exports the variables in the global environment of this script
  clusterExport(cl, ls(envir = parent.frame()), envir = parent.frame()) # exports the global variables of a script when this script is sourced
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # when beta0 is not provided and beta needs to be estimated via MLE
  if(is.null(beta0)){
    test_stat_eta <- foreach(i = 1:B, .combine = c,
                             .packages = c('truncdist', 'VGAM'),
                             .options.snow = opts) %dopar%
      {
        set.seed(seeds[i])
        # bkg-only sample:
        bkg_samp <- rtrunc(n_bkg, a = l, b = u, spec = 'gamma',
                           rate = bkg_rate, shape = bkg_shape)
        
        # physics-sample:
        s_samp <- rtrunc(n_phys, a = l, b = u, spec = 'norm',
                         mean = mean_sig, sd = sd_sig)
        b_samp <- rtrunc(n_phys, a = l, b = u, spec = 'gamma',
                         rate = bkg_rate, shape = bkg_shape)
        u_mask <- runif(n_phys)
        phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
        
        # obtaining MLE of beta:
        opt <- nlminb(start = 0.01,
                      objective = neg_loglik,
                      lower = 0, upper =10,
                      data = bkg_samp)
        beta_hat <- opt$par
        # ||S_\beta||_{G_\beta} evaluated at MLE of \beta
        norm_S <- integrate(normS_integrand, lower = l, upper = u,
                            beta = beta_hat)$value |> sqrt()
        # derivative of (||S_\beta||_{G_\beta})^2 w.r.t. \beta at MLE of \beta
        d_normS_sq <- -integrate(function(x){
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          g <- dtrunc(x, spec = 'pareto', a = l, b = u,
                      scale = l, shape = beta_hat)
          d_log_g <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          return(((fs^2)/g)*d_log_g)
        },l, u)$value
        # S0_{\hat\beta} evaluated on the physics sample
        S0_phys_vec <- sapply(phys_samp, 
                              function(x){
                                fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                             mean = mean_sig, sd = sd_sig)
                                g <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                            scale = l, shape = beta_hat)
                                return((fs/g-1)/(norm_S^2))
                              })
        # S0_{\hat\beta} evaluated on the bkg sample
        S0_bkg_vec <- sapply(bkg_samp,
                             function(x){
                               fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                            mean = mean_sig, sd = sd_sig)
                               g <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                           scale = l, shape = beta_hat)
                               return((fs/g-1)/(norm_S^2))
                             })
        # derivative of S0_{\beta} w.r.t. \beta at \hat\beta
        d_S0 <- function(x){
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          g <- dtrunc(x, spec = 'pareto', a = l, b = u,
                      scale = l, shape = beta_hat)
          d_log_g <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          num <- -((norm_S^2)*(fs/g)*d_log_g + (fs/g-1)*d_normS_sq)
          denom <- norm_S^4
          return(num/denom)
        }
        
        # computing numerator of the test statistic
        theta_0_hat <- mean(S0_phys_vec)
        delta_0_hat <- mean(S0_bkg_vec)
        eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat) # estimate of eta
        test_num <- sqrt(n_phys*n_bkg)*(eta_hat)
        
        # computing denominator of the test statistic
        J_hat <- -(-1/(beta_hat^2) - 
                     ((log(l)^2)*l^(-beta_hat) - (log(u)^2)*u^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)) -
                     ((log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)))^2)
        V_hat <- sapply(bkg_samp,
                        function(x){
                          val <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
                          return(val^2)
                        } ) |> mean()
        d_theta_hat <- sapply(phys_samp, d_S0) |> mean() # derivative of \hat\theta0_{\beta} w.r.t. \beta at \hat\beta
        d_delta_hat <- sapply(bkg_samp, d_S0) |> mean() # derivative of \hat\delta0_{\beta} w.r.t. \beta at \hat\beta
        d_theta_T <- 1/(1-delta_0_hat) # derivative of T w.r.t. its first component
        d_delta_T <- (theta_0_hat - 1)/((1-delta_0_hat)^2) # derivative of T w.r.t. its second component
        
        # d_log g_{\beta} w.r.t. \beta evaluated at bkg sample
        d_log_g_bkg <- sapply(bkg_samp, function(x){
          val <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          return(val)
        })
        cov_term <- mean(S0_bkg_vec*d_log_g_bkg)
        var_S0_F_hat <- mean(S0_phys_vec^2) - theta_0_hat^2
        var_S0_Fb_hat <- mean(S0_bkg_vec^2) - delta_0_hat^2
        denom1 <- n_bkg*var_S0_F_hat*(d_theta_T^2)
        denom2 <- n_phys*var_S0_Fb_hat*(d_delta_T^2)
        denom3 <- n_phys*(V_hat/J_hat^2)*(d_theta_T*d_theta_hat + d_delta_T*d_delta_hat)^2
        denom4 <- 2*(n_phys/J_hat)*d_delta_T*cov_term*(d_theta_T*d_theta_hat + d_delta_T*d_delta_hat)
        test_denom <- sqrt(denom1 + denom2 + denom3 + denom4) # denominator of the test statistic
        
        test_num/test_denom # test statistic
      }
  }else{ # when beta is fixed at beta0
    
    # ||S_\beta||_{G_\beta} at \beta = \beta0
    norm_S <- integrate(function(x) {
      fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                   mean = mean_sig, sd = sd_sig)
      g <- dtrunc(x, spec = 'pareto', a = l, b = u,
                  scale = l, shape = beta0)
      S_val <- (fs/g-1)
      return((S_val^2)*g)
    },l, u)$value |> sqrt()
    
    test_stat_eta <- foreach(i = 1:B, .combine = c,
                             .packages = c('truncdist', 'VGAM'),
                             .options.snow = opts) %dopar%
      {
        # background sample:
        bkg_samp <- rtrunc(n_bkg, a = l, b = u, spec = 'gamma',
                           rate = bkg_rate, shape = bkg_shape)
        
        # physics-sample:
        s_samp <- rtrunc(n_phys, a = l, b = u, spec = 'norm',
                         mean = mean_sig, sd = sd_sig)
        b_samp <- rtrunc(n_phys, a = l, b = u, spec = 'gamma',
                         rate = bkg_rate, shape = bkg_shape)
        u_mask <- runif(n_phys)
        phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
        
        # S0_{\beta0} evaluated on the physics sample
        S0_phys_vec <- sapply(phys_samp, 
                              function(x){
                                fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                             mean = mean_sig, sd = sd_sig)
                                g <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                            scale = l, shape = beta0)
                                S_val <- fs/g - 1
                                return((fs/g-1)/(norm_S^2))
                              })
        # S0_{\beta0} evaluated on the bkg sample
        S0_bkg_vec <- sapply(bkg_samp, 
                             function(x){
                               fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                            mean = mean_sig, sd = sd_sig)
                               g <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                           scale = l, shape = beta0)
                               S_val <- fs/g - 1
                               return((fs/g-1)/(norm_S^2))
                             })
        # calculating numerator of the test statistic:
        theta_0_hat <- mean(S0_phys_vec)
        delta_0_hat <- mean(S0_bkg_vec)
        eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat) # estimate of eta
        test_num <- sqrt(n_phys*n_bkg)*(eta_hat)
        
        # calculating denominator of the test statistic:
        sig_theta0_hat_sq <- mean(S0_phys_vec^2) - theta_0_hat^2
        sig_delta0_hat_sq <- mean(S0_bkg_vec^2) - delta_0_hat^2
        test_denom <- sqrt(
          n_bkg*sig_theta0_hat_sq/((1- delta_0_hat)^2) + 
            n_phys*sig_delta0_hat_sq*((theta_0_hat-1)^2)/((1-delta_0_hat)^4)
        ) # denominator of the test statistic
        
        test_num/test_denom # test statistic
      }
  }
  stopCluster(cl)
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE) # calculating p-values
  simulated_power <- mean(p_vals < signif.level) # proportion of rejections
  return(simulated_power)
}

# power simulation function in the presence of a bkg data using a uniform background:
sim_power_wbkg_unif <- function(eta, n_phys,
                                r, nsims, seed = 12345,
                                signif.level = 0.05){
  ####################### DESCRIPTION OF ARGUMENTS #############################
  # eta: true value of eta for simulation
  # n_phys: Size of the physics data
  # r: bkg data to physics data raio
  # nsims: number of iterations
  # seed: seed for reproducibility
  # signif.level: significance level for the test
  ##############################################################################
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f', eta))
  cat("\n--------------------------------------\n")
  
  B <- nsims; n_bkg <- r*n_phys
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  
  # parallelization for faster implementation
  cl <- makeCluster(8)
  clusterEvalQ(cl, {
    library(truncdist)
    library(VGAM)
  })
  registerDoSNOW(cl)
  clusterExport(cl, ls(envir = environment()), envir = environment()) # exports the variables in the global environment of this script
  clusterExport(cl, ls(envir = parent.frame()), envir = parent.frame()) # exports the global variables of a script when this script is sourced
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # ||S||_{G}:
  norm_S <- integrate(function(x) {
    fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                 mean = mean_sig, sd = sd_sig)
    g <- dunif(x, l, u)
    S_val <- (fs/g-1)
    return((S_val^2)*g)
  },l, u)$value |> sqrt()
  test_stat_eta <- foreach(i = 1:B, .combine = c,
                           .packages = c('truncdist', 'VGAM'),
                           .options.snow = opts) %dopar%
    {
      # bkg sample:
      bkg_samp <- rtrunc(n_bkg, a = l, b = u, spec = 'gamma',
                         rate = bkg_rate, shape = bkg_shape)
      
      # physics-sample:
      s_samp <- rtrunc(n_phys, a = l, b = u, spec = 'norm',
                       mean = mean_sig, sd = sd_sig)
      b_samp <- rtrunc(n_phys, a = l, b = u, spec = 'gamma',
                       rate = bkg_rate, shape = bkg_shape)
      u_mask <- runif(n_phys)
      phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
      
      # S0 evaluated on the physics sample
      S0_phys_vec <- sapply(phys_samp, 
                            function(x){
                              fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                           mean = mean_sig, sd = sd_sig)
                              g <- dunif(x, l, u)
                              S_val <- fs/g - 1
                              return((fs/g-1)/(norm_S^2))
                            })
      # S0 evaluated on the bkg sample
      S0_bkg_vec <- sapply(bkg_samp, 
                           function(x){
                             fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                          mean = mean_sig, sd = sd_sig)
                             g <- dunif(x, l, u)
                             S_val <- fs/g - 1
                             return((fs/g-1)/(norm_S^2))
                           })
      # calculating the numerator of the test statistic:
      theta_0_hat <- mean(S0_phys_vec)
      delta_0_hat <- mean(S0_bkg_vec)
      eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat) # estimate of eta
      test_num <- sqrt(n_phys*n_bkg)*(eta_hat)
      
      # calculating the denominator of the test statistic:
      sig_theta0_hat_sq <- mean(S0_phys_vec^2) - theta_0_hat^2
      sig_delta0_hat_sq <- mean(S0_bkg_vec^2) - delta_0_hat^2
      test_denom <- sqrt(
        n_bkg*sig_theta0_hat_sq/((1- delta_0_hat)^2) + 
          n_phys*sig_delta0_hat_sq*((theta_0_hat-1)^2)/((1-delta_0_hat)^4)
      )
      # test statistic:
      test_num/test_denom
    }
  stopCluster(cl)
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE) # calculating p-values
  simulated_power <- mean(p_vals < signif.level) # proportion of rejections
  return(simulated_power)
}

# power simulation function in the absence of a bkg data using Pareto type I density as the benchmark model:
sim_power_wobkg <- function(eta, n_phys, lambda,
                            mean1_in_g, mean2_in_g, sd_in_g,
                            nsims, seed = 12345,
                            signif.level = 0.05){
  ####################### DESCRIPTION OF ARGUMENTS #############################
  # eta: true value of eta for simulation
  # n_phys: Size of the physics data
  # mean1_in_g: location of the first Gaussian component in g
  # mean2_in_g: location of the second Gaussian component in g
  # sd_in_g: scale parameter for the Gaussian components in g
  # nsims: number of iterations
  # seed: seed for reproducibility
  # signif.level: significance level for the test
  ##############################################################################
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f and lambda = %.4f', eta, lambda))
  cat("\n--------------------------------------\n")
  
  B <- nsims
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  # parallelization for faster implementation
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  clusterExport(cl, ls(envir = environment()), envir = environment()) # exports the variables in the global environment of this script
  clusterExport(cl, ls(envir = parent.frame()), envir = parent.frame()) # exports the global variables of a script when this script is sourced
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  test_stat_eta <- foreach(i = 1:B, .combine = c,
                           .packages = c('truncdist', 'VGAM'),
                           .options.snow = opts) %dopar%
    {
      set.seed(seeds[i])
      
      # physics-sample:
      s_samp <- rtrunc(n_phys, a = l, b = u, spec = 'norm',
                       mean = mean_sig, sd = sd_sig)
      b_samp <- rtrunc(n_phys, a = l, b = u, spec = 'gamma',
                       rate = bkg_rate, shape = bkg_shape)
      u_mask <- runif(n_phys)
      phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
      
      # obtaining \hat\alpha using the benchmark model on the physics sample
      opt <- nlminb(start = 0.01,
                    objective = neg_loglik,
                    lower = 0, upper = 10,
                    data = phys_samp)
      alpha_hat <- opt$par
      # ||S_\beta||_{G_\beta} evaluated at \hat\beta = (\hat\alpha, \lambda)
      norm_S <- integrate(function(x) {
        fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                     mean = mean_sig, sd = sd_sig)
        q <- dtrunc(x, spec = 'pareto', a = l, b = u,
                    scale = l, shape = alpha_hat)
        phi1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean1_in_g, sd = sd_in_g)
        phi2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean2_in_g, sd = sd_in_g)
        g <- lambda*(phi1 + phi2) + (1-2*lambda)*q
        return(((fs/g-1)^2)*g)
      },l, u)$value |> sqrt()
      # derivative of (||S_{\beta}||_{G_\beta})^2 w.r.t. \alpha at \beta =  (\hat\alpha, \lambda)
      d_normS_sq <- -(1-2*lambda)*integrate(function(x){
        fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                     mean = mean_sig, sd = sd_sig)
        q <- dtrunc(x, spec = 'pareto', a = l, b = u,
                    scale = l, shape = alpha_hat)
        phi1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean1_in_g, sd = sd_in_g)
        phi2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean2_in_g, sd = sd_in_g)
        g <- lambda*(phi1 + phi2) + (1-2*lambda)*q
        d_log_q <- 1/alpha_hat - log(x) - (log(u)*u^(-alpha_hat) - log(l)*l^(-alpha_hat))/(l^(-alpha_hat) - u^(-alpha_hat))
        return(((fs/g)^2)*q*d_log_q)
      },l, u)$value
      
      # S0_{\hat\beta} evaluated on the physics sample
      S0_phys_vec <- sapply(phys_samp, 
                            function(x){
                              fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                           mean = mean_sig, sd = sd_sig)
                              q <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                          scale = l, shape = alpha_hat)
                              phi1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                             mean = mean1_in_g, sd = sd_in_g)
                              phi2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                             mean = mean2_in_g, sd = sd_in_g)
                              g <- lambda*(phi1 + phi2) + (1-2*lambda)*q
                              return((fs/g-1)/(norm_S^2))
                            })
      # derivative of the function S0 w.r.t. the parameter \alpha, evaluated at \hat\beta = (\hat\alpha, \lambda)
      d_S0 <- function(x){
        fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                     mean = mean_sig, sd = sd_sig)
        q <- dtrunc(x, spec = 'pareto', a = l, b = u,
                    scale = l, shape = alpha_hat)
        phi1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean1_in_g, sd = sd_in_g)
        phi2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean2_in_g, sd = sd_in_g)
        g <- lambda*(phi1 + phi1) + (1-2*lambda)*q
        d_log_q <- 1/alpha_hat - log(x) - (log(u)*u^(-alpha_hat) - log(l)*l^(-alpha_hat))/(l^(-alpha_hat) - u^(-alpha_hat))
        num <- -((norm_S^2)*(fs/(g^2))*q*(1-2*lambda)*d_log_q + (fs/g-1)*d_normS_sq)
        denom <- norm_S^4
        return(num/denom)
      }
      
      # Derivative of log q_\alpha at \hat\alpha evaluated on the physics sample
      d_log_q_vec <- sapply(phys_samp, function(x){
        1/alpha_hat - log(x) - (log(u)*u^(-alpha_hat) - log(l)*l^(-alpha_hat))/(l^(-alpha_hat) - u^(-alpha_hat))
      })
      # conservative estimate of eta, \hat\theta0:
      theta_0_hat <- mean(S0_phys_vec)
      
      #components for the denominator of the test statistic:
      J_hat <- -(-1/(alpha_hat^2) - 
                   ((log(l)^2)*l^(-alpha_hat) - (log(u)^2)*u^(-alpha_hat))/(l^(-alpha_hat) - u^(-alpha_hat)) -
                   ((log(u)*u^(-alpha_hat) - log(l)*l^(-alpha_hat))/(l^(-alpha_hat) - u^(-alpha_hat)))^2)
      V_hat <- sum(d_log_q_vec^2)/n_phys
      d_theta0_hat <- sapply(phys_samp, d_S0) |> mean() # derivative of \hat\theta0 at \hat\beta
      var_S0_F_hat <- mean(S0_phys_vec^2) - theta_0_hat^2
      sig_theta0_hat <- sqrt(
        var_S0_F_hat +  (d_theta0_hat^2)*V_hat/(J_hat^2) + 
          (2/J_hat)*d_theta0_hat*sum(S0_phys_vec*d_log_q_vec)/n_phys
      )
      # test statistic
      sqrt(n_phys)*theta_0_hat/sig_theta0_hat
    }
  stopCluster(cl)
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE) # calculating p-values
  simulated_power <- mean(p_vals < signif.level) # proportion of rejections
  return(simulated_power)
}

