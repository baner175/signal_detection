rm(list = ls())
rm(list = ls())
library(truncdist)
library(VGAM)
library(doSNOW)
library(parallel)
library(foreach)

l <- 1; u <- 2

B <- 1e2 # In paper, we use 100,000
n_samp <- 5e3
eta_true <- 0; lambda <- 0
beta0 <- 4

################################################################
################ SIGNAL AND SIGNAL REGION ######################
################################################################

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02


# signal density
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb_true <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)
# true mixture density
f <- function(x) eta_true*fs(x)+(1-eta_true)*fb_true(x)

gb <- function(x, beta = beta0){
  lambda*fs(x) + (1-lambda)*dtrunc(x, spec = 'pareto', a = l, b = u,
                                   scale = l, shape = beta)
}

S <- function(x, beta = beta0) fs(x)/gb(x, beta) - 1
norm_S <- function(beta = beta0) integrate(function(x) (S(x, beta)^2)*gb(x,beta),
                                           l, u)$value |> sqrt()

neg_ll <- function(eta, data){
  fi <- eta*dtrunc(data, spec = 'norm', a = l, b = u,
                   mean = mean_sig, sd = sd_sig) + 
    (1-eta)*(
      lambda*dtrunc(data, spec = 'norm', a = l, b = u,
                    mean = mean_sig, sd = sd_sig) + 
        (1-lambda)*dtrunc(data, spec = 'pareto', a = l, b = u,
                          scale = l, shape = beta0)
    )
  return(-sum(log(fi)))
}

set.seed(12345)
seeds <- sample.int(.Machine$integer.max, B)
cl <- makeCluster(8)
registerDoSNOW(cl)
pb <- txtProgressBar(max = B, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start_time <- Sys.time()
test_stat_LRT <- foreach(i = 1:B, .combine = rbind,
                         .packages = c('truncdist', 'VGAM'),
                         .options.snow = opts) %dopar%
  {
    set.seed(seeds[i])
    
    # physics-sample:
    s_samp <- rtrunc(n_samp, a = l, b = u, spec = 'norm',
                     mean = mean_sig, sd = sd_sig)
    b_samp <- rtrunc(n_samp, a = l, b = u, spec = 'gamma',
                     rate = bkg_rate, shape = bkg_shape)
    u_mask <- runif(n_samp)
    phys_samp <- ifelse(u_mask <= eta_true, s_samp, b_samp)
    
    ll_res <- nlminb(start = 0.01,
                     objective = neg_ll,
                     lower = -Inf, upper = Inf,
                     data = phys_samp)
    eta_hat <- ll_res$par
    eta_hat_C <- (eta_hat>0)*eta_hat
    ll_eta_hat_C <- -neg_ll(eta = eta_hat_C, data = phys_samp)
    ll_0 <- -neg_ll(eta = 0, data = phys_samp)
    
    test_stat_eta_hat_C_0 <- 2*(ll_eta_hat_C - ll_0)
    c(test_stat_eta_hat_C_0)
  }
close(pb)

Sys.time() - start_time
stopCluster(cl)

df <- data.frame('test_stat_LRT' = test_stat_LRT)

file_name <- paste0('Numerical_example_Sec2.1/',
                    'LRT',
                    '_B(', B,
                    ')_beta0(', beta0,
                    ')_n_samp(', n_samp,
                    ')_eta(', eta_true,
                    ')_lambda(', lambda,
                    ').csv')
write.csv(df, file_name,
          row.names = FALSE)
