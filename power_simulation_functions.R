library(truncdist)
library(VGAM)
library(doSNOW)
library(parallel)
library(foreach)

#-------------------------------------------------------------------------------
## background fitting functions ##

qb_bkg_model_unbinned <- function(beta, data){
  qb_i <- sapply(data, dtrunc, spec = 'pareto', a = l, b = u,
                 scale = l, shape = beta)
  return(-sum(log(qb_i)))
}

normS_integrand <- function(x, beta){
  fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
               mean = mean_sig, sd = sd_sig)
  qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
               scale = l, shape = beta)
  return(((fs/qb-1)^2)*qb)
}

#-------------------------------------------------------------------------------
## simulation functions ##

sim_power_wbkg <- function(eta, n_phys,
                                     r, nsims, seed = 12345,
                                     signif.level = 0.05,
                                     beta0 = NULL){
  
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f', eta))
  cat("\n--------------------------------------\n")
  
  B <- nsims; n_bkg <- r*n_phys
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  clusterExport(cl, ls(envir = environment()), envir = environment()) # exports the variables in the global environment of this script
  clusterExport(cl, ls(envir = parent.frame()), envir = parent.frame()) # exports the global variables of a script when this script is sourced
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
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
        
        opt <- nlminb(start = 0.01,
                      objective = qb_bkg_model_unbinned,
                      lower = 0, upper =10,
                      data = bkg_samp)
        beta_hat <- opt$par
        norm_S <- integrate(normS_integrand, lower = l, upper = u,
                            beta = beta_hat)$value |> sqrt()
        d_normS2 <- -integrate(function(x){
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          return(((fs^2)/qb)*d_log_qb)
        },l, u)$value
        
        S2_phys_vec <- sapply(phys_samp, 
                              function(x){
                                fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                             mean = mean_sig, sd = sd_sig)
                                qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                             scale = l, shape = beta_hat)
                                return((fs/qb-1)/(norm_S^2))
                              })
        S2_bkg_vec <- sapply(bkg_samp,
                             function(x){
                               fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                            mean = mean_sig, sd = sd_sig)
                               qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                            scale = l, shape = beta_hat)
                               return((fs/qb-1)/(norm_S^2))
                             })
        d_S2 <- function(x){
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          num <- -((norm_S^2)*(fs/qb)*d_log_qb + (fs/qb-1)*d_normS2)
          denom <- norm_S^4
          return(num/denom)
        }
        
        theta_0_hat <- mean(S2_phys_vec)
        delta_0_hat <- mean(S2_bkg_vec)
        
        J_hat <- -(-1/(beta_hat^2) - 
                     ((log(l)^2)*l^(-beta_hat) - (log(u)^2)*u^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)) -
                     ((log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)))^2)
        V_hat <- sapply(bkg_samp,
                        function(x){
                          val <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
                          return(val^2)
                        } ) |> mean()
        eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
        d_theta_hat <- sapply(phys_samp, d_S2) |> mean()
        d_delta_hat <- sapply(bkg_samp, d_S2) |> mean()
        d_theta_T <- 1/(1-delta_0_hat)
        d_delta_T <- (theta_0_hat - 1)/((1-delta_0_hat)^2)
        d_log_qb_bkg <- sapply(bkg_samp, function(x){
          val <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          return(val)
        })
        cov_term <- mean(S2_bkg_vec*d_log_qb_bkg)
        
        var_S2_F_hat <- mean(S2_phys_vec^2) - theta_0_hat^2
        var_S2_Fb_hat <- mean(S2_bkg_vec^2) - delta_0_hat^2
        
        
        test_num <- sqrt(n_phys*n_bkg)*(eta_hat)
        
        denom1 <- n_bkg*var_S2_F_hat*(d_theta_T^2)
        denom2 <- n_phys*var_S2_Fb_hat*(d_delta_T^2)
        denom3 <- n_phys*(V_hat/J_hat^2)*(d_theta_T*d_theta_hat + d_delta_T*d_delta_hat)^2
        denom4 <- 2*(n_phys/J_hat)*d_delta_T*cov_term*(d_theta_T*d_theta_hat + d_delta_T*d_delta_hat)
        
        test_denom <- sqrt(denom1 + denom2 + denom3 + denom4)
        
        test_num/test_denom
      }
  }else{
    norm_S <- integrate(function(x) {
      fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                   mean = mean_sig, sd = sd_sig)
      qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta0)
      S_val <- (fs/qb-1)
      return((S_val^2)*qb)
    },l, u)$value |> sqrt()
    test_stat_eta <- foreach(i = 1:B, .combine = c,
                             .packages = c('truncdist', 'VGAM'),
                             .options.snow = opts) %dopar%
      {
        bkg_samp <- rtrunc(n_bkg, a = l, b = u, spec = 'gamma',
                           rate = bkg_rate, shape = bkg_shape)
        
        # physics-sample:
        s_samp <- rtrunc(n_phys, a = l, b = u, spec = 'norm',
                         mean = mean_sig, sd = sd_sig)
        b_samp <- rtrunc(n_phys, a = l, b = u, spec = 'gamma',
                         rate = bkg_rate, shape = bkg_shape)
        u_mask <- runif(n_phys)
        phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
        
        S2_phys_vec <- sapply(phys_samp, 
                              function(x){
                                fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                             mean = mean_sig, sd = sd_sig)
                                gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                             scale = l, shape = beta0)
                                S_val <- fs/gb - 1
                                return((fs/gb-1)/(norm_S^2))
                              })
        S2_bkg_vec <- sapply(bkg_samp, 
                             function(x){
                               fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                            mean = mean_sig, sd = sd_sig)
                               gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                            scale = l, shape = beta0)
                               S_val <- fs/gb - 1
                               return((fs/gb-1)/(norm_S^2))
                             })
        theta_0_hat <- mean(S2_phys_vec)
        delta_0_hat <- mean(S2_bkg_vec)
        
        sig_theta0_hat_sq <- mean(S2_phys_vec^2) - theta_0_hat^2
        sig_delta0_hat_sq <- mean(S2_bkg_vec^2) - delta_0_hat^2
        
        eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
        
        test_num <- sqrt(n_phys*n_bkg)*(eta_hat)
        test_denom <- sqrt(
          n_bkg*sig_theta0_hat_sq/((1- delta_0_hat)^2) + 
            n_phys*sig_delta0_hat_sq*((theta_0_hat-1)^2)/((1-delta_0_hat)^4)
        )
        test_num/test_denom
      }
  }
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE)
  simulated_power <- mean(p_vals < signif.level)
  return(simulated_power)
}

sim_power_wbkg_unif <- function(eta, n_phys,
                                 r, nsims, seed = 12345,
                                 signif.level = 0.05){
  
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f', eta))
  cat("\n--------------------------------------\n")
  
  B <- nsims; n_bkg <- r*n_phys
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
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
  
  norm_S <- integrate(function(x) {
    fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                 mean = mean_sig, sd = sd_sig)
    qb <- dunif(x, l, u)
    S_val <- (fs/qb-1)
    return((S_val^2)*qb)
  },l, u)$value |> sqrt()
  test_stat_eta <- foreach(i = 1:B, .combine = c,
                           .packages = c('truncdist', 'VGAM'),
                           .options.snow = opts) %dopar%
    {
      bkg_samp <- rtrunc(n_bkg, a = l, b = u, spec = 'gamma',
                         rate = bkg_rate, shape = bkg_shape)
      
      # physics-sample:
      s_samp <- rtrunc(n_phys, a = l, b = u, spec = 'norm',
                       mean = mean_sig, sd = sd_sig)
      b_samp <- rtrunc(n_phys, a = l, b = u, spec = 'gamma',
                       rate = bkg_rate, shape = bkg_shape)
      u_mask <- runif(n_phys)
      phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
      
      S2_phys_vec <- sapply(phys_samp, 
                            function(x){
                              fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                           mean = mean_sig, sd = sd_sig)
                              gb <- dunif(x, l, u)
                              S_val <- fs/gb - 1
                              return((fs/gb-1)/(norm_S^2))
                            })
      S2_bkg_vec <- sapply(bkg_samp, 
                           function(x){
                             fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                          mean = mean_sig, sd = sd_sig)
                             gb <- dunif(x, l, u)
                             S_val <- fs/gb - 1
                             return((fs/gb-1)/(norm_S^2))
                           })
      theta_0_hat <- mean(S2_phys_vec)
      delta_0_hat <- mean(S2_bkg_vec)
      
      sig_theta0_hat_sq <- mean(S2_phys_vec^2) - theta_0_hat^2
      sig_delta0_hat_sq <- mean(S2_bkg_vec^2) - delta_0_hat^2
      
      eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
      
      test_num <- sqrt(n_phys*n_bkg)*(eta_hat)
      test_denom <- sqrt(
        n_bkg*sig_theta0_hat_sq/((1- delta_0_hat)^2) + 
          n_phys*sig_delta0_hat_sq*((theta_0_hat-1)^2)/((1-delta_0_hat)^4)
      )
      test_num/test_denom
    }
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE)
  simulated_power <- mean(p_vals < signif.level)
  return(simulated_power)
}

sim_power_wobkg <- function(eta, n_phys, lambda,
                            mean1_in_gb, mean2_in_gb,
                            sd_in_gb, nsims, seed = 12345,
                            signif.level = 0.05){
  
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f and lambda = %.4f', eta, lambda))
  cat("\n--------------------------------------\n")
  
  B <- nsims
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  
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
      
      opt <- nlminb(start = 0.01,
                    objective = qb_bkg_model_unbinned,
                    lower = 0, upper = 10,
                    data = phys_samp)
      beta_hat <- opt$par
      norm_S <- integrate(function(x) {
        fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                     mean = mean_sig, sd = sd_sig)
        qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                     scale = l, shape = beta_hat)
        fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                      mean = mean1_in_gb, sd = sd_in_gb)
        fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                      mean = mean2_in_gb, sd = sd_in_gb)
        gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
        return(((fs/gb-1)^2)*gb)
      },l, u)$value |> sqrt()
      d_normS2 <- -(1-2*lambda)*integrate(function(x){
        fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                     mean = mean_sig, sd = sd_sig)
        qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                     scale = l, shape = beta_hat)
        fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                      mean = mean1_in_gb, sd = sd_in_gb)
        fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                      mean = mean2_in_gb, sd = sd_in_gb)
        gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
        d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
        return(((fs/gb)^2)*qb*d_log_qb)
      },l, u)$value
      
      S2_phys_vec <- sapply(phys_samp, 
                            function(x){
                              fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                           mean = mean_sig, sd = sd_sig)
                              qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                           scale = l, shape = beta_hat)
                              fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                            mean = mean1_in_gb, sd = sd_in_gb)
                              fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                            mean = mean2_in_gb, sd = sd_in_gb)
                              gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
                              return((fs/gb-1)/(norm_S^2))
                            })
      d_S2 <- function(x){
        fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                     mean = mean_sig, sd = sd_sig)
        qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                     scale = l, shape = beta_hat)
        fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                      mean = mean1_in_gb, sd = sd_in_gb)
        fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                      mean = mean2_in_gb, sd = sd_in_gb)
        gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
        d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
        num <- -((norm_S^2)*(fs/(gb^2))*qb*(1-2*lambda)*d_log_qb + (fs/gb-1)*d_normS2)
        denom <- norm_S^4
        return(num/denom)
      }
      
      d_log_qb_vec <- sapply(phys_samp, function(x){
        1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
      })
      
      theta_0_hat <- mean(S2_phys_vec)
      
      J_hat <- -(-1/(beta_hat^2) - 
                   ((log(l)^2)*l^(-beta_hat) - (log(u)^2)*u^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)) -
                   ((log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)))^2)
      V_hat <- sum(d_log_qb_vec^2)/n_phys
      
      d_theta0_hat <- sapply(phys_samp, d_S2) |> mean()
      var_S2_F_hat <- mean(S2_phys_vec^2) - theta_0_hat^2
      sig_theta0_hat <- sqrt(
        var_S2_F_hat +  (d_theta0_hat^2)*V_hat/(J_hat^2) + 
          (2/J_hat)*d_theta0_hat*sum(S2_phys_vec*d_log_qb_vec)/n_phys
      )
      
      sqrt(n_phys)*theta_0_hat/sig_theta0_hat
    }
  stopCluster(cl)
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE)
  simulated_power <- mean(p_vals < signif.level)
  return(simulated_power)
}

