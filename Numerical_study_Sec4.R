rm(list = ls())
library(truncdist)
## Global parameters ##

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02

l <- 1; u <- 2

# signal density and CDF
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)
Fs <- function(x) ptrunc(x, a = l, b = u, spec = 'norm', mean = mean_sig, sd = sd_sig)

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb_true <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)


source('power_simulation_functions.R')

B <- 1e2; n_phys <- 250
eta <- 0
eps <- 1e-3
lambda_seq <- c(0.00, 0.01, 0.02, 0.03)

# Calculating signal region:
# Figuring out (mu_s -d, mu_s + d) that covers an area of 1-eps:
find_d <- function(d)
{
  pl <- Fs(mean_sig - d)
  pu <- Fs(mean_sig + d)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))

r <- sol$root

M_lower <- mean_sig - r
M_upper <- mean_sig + r

# PROPOSAL BACKGROUND DENSITY PARAMETERS:
mean1_in_gb <- 0.5*M_lower + 0.5*mean_sig
mean2_in_gb <- 0.4*M_upper + 0.6*mean_sig
sd_in_gb <- 4*sd_sig

neg_exp_log_lik <- function(beta, lambda)
{
  -integrate(function(x){
    fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                      a = l, b = u,
                      spec = 'norm')
    fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                       a = l, b = u,
                       spec = 'norm')
    qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                 scale = l, shape = beta)
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    
    return(log(gb)*fb_true(x))
  }, l, u)$value
}
neg_exp_log_lik <- Vectorize(neg_exp_log_lik)

# simulation results with known beta:
pow_seq <- c(); delta_seq <- c()
for(lam in lambda_seq){
  pow_val <- sim_power_wobkg(eta = eta, n_phys = n_phys,
                             lambda = lam,
                             mean1_in_gb = mean1_in_gb,
                             mean2_in_gb = mean2_in_gb,
                             sd_in_gb = sd_in_gb,
                             nsims = B, seed = 12345,
                             signif.level = 0.05)
  pow_seq <- c(pow_seq, pow_val)
  
  beta_star <- nlminb(start = 0.01, 
                      objective = neg_exp_log_lik,
                      lambda = lam,
                      lower = 1e-6,
                      upper = 10)$par
  normS <- integrate(function(x){
    fs <- fs(x)
    fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                      a = l, b = u,
                      spec = 'norm')
    fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                       a = l, b = u,
                       spec = 'norm')
    qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                 scale = l, shape = beta_star)
    gb <- lam*(fs_val1+fs_val2) + (1-2*lam)*qb
    S <- (fs/gb - 1)
    return((S^2)*gb)
  }, l, u)$value |> sqrt()
  
  delta_val <- integrate(function(x){
    fs <- fs(x)
    fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                      a = l, b = u,
                      spec = 'norm')
    fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                       a = l, b = u,
                       spec = 'norm')
    qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                 scale = l, shape = beta_star)
    gb <- lam*(fs_val1+fs_val2) + (1-2*lam)*qb
    S1 <- (fs/gb-1)/normS
    return(S1*fb_true(x))
  }, l, u)$value
  delta_seq <- c(delta_seq, delta_val)
  
}

power_res <- data.frame(lambda = lambda_seq,
                        sim_power = pow_seq,
                        delta = delta_seq)
colnames(power_res) <- c('lambda', 'power', 'delta')

file_name <- paste0('Numerical_example_Sec4/',
                    'WOBKG__',
                    'eta(',eta,')_',
                    'lambda(', paste0(lambda_seq, collapse = '_'),')_',
                    'B(',B,')_',
                    'n_phys(',n_phys,').csv')

write.csv(x = power_res, file_name,
          row.names = FALSE)
