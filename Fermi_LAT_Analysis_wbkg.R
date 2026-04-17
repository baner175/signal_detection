rm(list = ls())

################################################################################
################################################################################
##                                                                            ##
##     Code for reproducing the analysis on the simulated data                ##
##     from Fermi-LAT in the presence of a background-only sample             ##
##                                                                            ##
################################################################################
################################################################################

####################### Loading the required packages ##########################
library(VGAM)
library(truncdist)
################################################################################

# Loading the data and transforming into log-scale
phys_data <- read.table('Fermi_LAT_physics.txt', header = TRUE)
bkg_data <- read.table('Fermi_LAT_bkg_only.txt', header = TRUE)
y <- log(bkg_data$x)
x <- log(phys_data$x)
n <- length(x)
m <- length(y)

real_l <- 1; real_u <- 35 # Search region
l <- log(real_l); u <- log(real_u) # Search region on log-scale

# parameters for the signal density
mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)
# SIGNAL DENSITY:
fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd)*exp(x))
}

################################################################################
######## PARAMETRIC MODEL: TRUNCATED EXPONENTIAL ###############################
################################################################################

g_model <- function(beta){ # defining log-likelihood based g_beta
  gi <- dtrunc(y, spec = 'exp', rate = beta, a = l, b = u)
  return(-sum(log(gi)))
}
beta_hat <- nlminb(start = 0.01,
                   objective = g_model,
                   upper = 10, lower = 0)$par
# Defining proposal background g_beta at MLE beta_hat:
g <- function(x) dtrunc(x, spec = 'exp', rate = beta_hat, a = l, b = u)

norm_S <- integrate(function(t){
  fs <- fs(t)
  g <- g(t)
  return(((fs/g - 1)^2)*g)
}, l, u)$value |> sqrt() # ||S||_{G_\beta} at \beta = \hat{\beta}
S0_phys_vec <- sapply(x, function(t){ # Evaluating S_0 on the physics data
  fs <- fs(t)
  g <- g(t)
  S_val <- (fs/g - 1)
  return(S_val/(norm_S^2))
})
S0_bkg_vec <- sapply(y, function(t){ # Evaluating S_0 on the bkg data
  fs <- fs(t)
  g <- g(t)
  S_val <- (fs/g - 1)
  return(S_val/(norm_S^2))
})
theta0_hat <- mean(S0_phys_vec); delta0_hat <- mean(S0_bkg_vec)
eta_hat_exp <- (theta0_hat-delta0_hat)/(1-delta0_hat) # estimate of eta
test_num <- sqrt(m*n)*eta_hat_exp # numerator of the test statistic

###### Now we shall compute the components involved ###############
###### in the denominator of the test statistic ###################

# derivative of (||S||_G)^2 w.r.t. \beta
d_normS_sq <- -integrate(function(t){
  fs <- fs(t)
  g <- g(t)
  d_log_g <- (1/beta_hat) - t - 
    (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u))
  return((fs^2)*d_log_g/g)
}, l, u)$value

d_S0 <- function(t){ # derivative of S_0 w.r.t. \beta
  fs <- fs(t)
  g <- g(t)
  d_log_g <- (1/beta_hat) - t - 
    (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u))
  
  return(-((norm_S^2)*(fs/g)*g*d_log_g + (fs/g-1)*d_normS_sq)/(norm_S^4))
}
d_log_g_yi <- sapply(y, function(t){ # evaluating d_\beta log(g) on the bkg data
  return(
    (1/beta_hat) - t - 
      (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u)) 
  )
})
# double derivative of log g_\beta w.r.t. \beta:
d2_log_g_y <- -(1/(beta_hat^2)) - (((l^2)*exp(-beta_hat*l) - (u^2)*exp(-beta_hat*u))/(exp(-beta_hat*l) - exp(-beta_hat*u)) - 
                                     ((u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u)))^2)

d_theta0 <- sapply(x, d_S0) |> mean() # derivative of \hat\theta_0 w.r.t. \beta
d_delta0 <- sapply(y, d_S0) |> mean() # derivative of \hat\delta_0 w.r.t. \beta
d_theta_T <- 1/(1-delta0_hat) # derivative of T w.r.t. first component
d_delta_T <- (theta0_hat - 1)/((1-delta0_hat)^2) # derivative of T w.r.t. second component
# components for the denominator
cov_term <- mean(d_log_g_yi*S0_bkg_vec)
V_hat <- mean((d_log_g_yi)^2)
J_hat <- -d2_log_g_y
var_S0_F_hat <- mean(S0_phys_vec^2) - (theta0_hat^2)
var_S0_Fb_hat <- mean(S0_bkg_vec^2) - (delta0_hat^2)
denom1 <- m*(d_theta_T^2)*var_S0_F_hat
denom2 <- n*(d_delta_T^2)*var_S0_Fb_hat
denom3 <- n*(V_hat/(J_hat^2))*((d_theta_T*d_theta0 + d_delta_T*d_delta0)^2)
denom4 <- 2*(n/J_hat)*d_delta_T*cov_term*(d_theta_T*d_theta0 + d_delta_T*d_delta0)
test_denom <- sqrt(denom1 + denom2 + denom3 + denom4) # denominator of the test statistic
test_stat <- test_num/test_denom # test statistic
p_val_exp <- pnorm(test_stat, lower.tail = FALSE) # p-value
sig_hat_exp <- test_denom/(sqrt(m+n))
std_err_exp <- sig_hat_exp*sqrt((m+n)/(m*n)) # standard error
ci_95_exp <- eta_hat_exp + c(-1,1)*qnorm(0.975)*std_err_exp # 95% CI

################################################################################
######## PARAMETRIC MODEL: TAIL OF TRUNCATED GAUSSIAN  #########################
################################################################################
mu_in_g <- -1; sigma_factor_in_g <- 2 # so that g_\beta(x) \propto exp(-((x+1)^2)/(4*\beta))
# below we replicate the same code as above with necessary adjustments for the new background model
g_model <- function(beta){
  gi <- dtrunc(y, spec = 'norm',
               mean = mu_in_g,
               sd = sqrt(sigma_factor_in_g*beta),
               a = l, b = u)
  return(-sum(log(gi)))
}
beta_hat <- nlminb(start = 0.01,
                   objective = g_model,
                   upper = 10, lower = 0)$par
g <- function(t) dtrunc(t, spec = 'norm',
                        mean = mu_in_g,
                        sd = sqrt(sigma_factor_in_g*beta_hat),
                        a = l, b = u)
d_log_h <- function(t) ((t-mu_in_g)^2)/(2*sigma_factor_in_g*(beta_hat^2))
d2_log_h <- function(t) (-(t-mu_in_g)^2)/(sigma_factor_in_g*(beta_hat^3))
E_g_d_log_h <- integrate(function(t) d_log_h(t)*g(t), l, u)$value
d2_log_g_int_1 <- integrate(function(y){
  d2_log_h <- d2_log_h(y)
  d_log_h <- d_log_h(y)
  d2_h_by_h <- d2_log_h + d_log_h^2
  g <- g(y)
  return(d2_h_by_h*g)
}, l, u)$value
d2_log_g_int_2 <- integrate(function(y){
  d_log_h <- d_log_h(y)
  g <- g(y)
  return(d_log_h*g)
}, l, u)$value
d_log_g <- function(t) d_log_h(t) - E_g_d_log_h
d2_log_g <- function(t){
  d2_log_h <- d2_log_h(t)
  int_val <- d2_log_g_int_1 - (d2_log_g_int_2^2)
  return(d2_log_h - int_val)
}
norm_S <- integrate(function(t){
  fs <- fs(t)
  g <- g(t)
  return(((fs/g - 1)^2)*g)
}, l, u)$value |> sqrt()
S0_phys_vec <- sapply(x, function(t){
  fs <- fs(t)
  g <- g(t)
  S_val <- (fs/g - 1)
  return(S_val/(norm_S^2))
})
S0_bkg_vec <- sapply(y, function(t){
  fs <- fs(t)
  g <- g(t)
  S_val <- (fs/g - 1)
  return(S_val/(norm_S^2))
})
theta0_hat <- mean(S0_phys_vec); delta0_hat <- mean(S0_bkg_vec)
eta_hat_gt <- (theta0_hat-delta0_hat)/(1-delta0_hat) # estimate of eta 
test_num <- sqrt(m*n)*eta_hat_gt # numerator of the test statistic
d_normS_sq <- -integrate(function(t){
  fs <- fs(t)
  g <- g(t)
  d_log_g <- d_log_g(t)
  return((fs^2)*d_log_g/g)
}, l, u)$value
d_S0 <- function(t){
  fs <- fs(t)
  g <- g(t)
  d_log_g <- d_log_g(t)
  
  return(-((norm_S^2)*(fs/g)*g*d_log_g + (fs/g-1)*d_normS_sq)/(norm_S^4))
}
d_log_g_yi <- sapply(y, d_log_g)
d2_log_g_yi <- sapply(y, d2_log_g)
d_theta0 <- sapply(x, d_S0) |> mean()
d_delta0 <- sapply(y, d_S0) |> mean()
d_theta_T <- 1/(1-delta0_hat)
d_delta_T <- (theta0_hat - 1)/((1-delta0_hat)^2)
# calculating the denominator of the test statistic
cov_term <- mean(d_log_g_yi*S0_bkg_vec)
V_hat <- mean((d_log_g_yi)^2)
J_hat <- -mean(d2_log_g_yi)
var_S0_F_hat <- mean(S0_phys_vec^2) - (theta0_hat^2)
var_S0_Fb_hat <- mean(S0_bkg_vec^2) - (delta0_hat^2)
denom1 <- m*(d_theta_T^2)*var_S0_F_hat
denom2 <- n*(d_delta_T^2)*var_S0_Fb_hat
denom3 <- n*(V_hat/(J_hat^2))*((d_theta_T*d_theta0 + d_delta_T*d_delta0)^2)
denom4 <- 2*(n/J_hat)*d_delta_T*cov_term*(d_theta_T*d_theta0 + d_delta_T*d_delta0)
test_denom <- sqrt(denom1 + denom2 + denom3 + denom4) # denominator of the test statistic
test_stat <- test_num/test_denom # test statistic
p_val_gt <- pnorm(test_stat, lower.tail = FALSE) # p-value
sig_hat_gt <- test_denom/(sqrt(m+n))
std_err_gt <- sig_hat_gt*sqrt((m+n)/(m*n)) # std error
ci_95_gt <- eta_hat_gt + c(-1,1)*qnorm(0.975)*std_err_gt # p-value

################################################################################
##############  UNIFORM BACKGROUND (NO PARAMETERS INVOLVED)   ##################
################################################################################
g <- function(t) dunif(t, min = l, max = u) # uniform proposal background
norm_S <- integrate(function(t){
  fs <- fs(t)
  g <- g(t)
  return(((fs/g - 1)^2)*g)
}, l, u)$value |> sqrt()
# similar code as above for the uniform background
S0_phys_vec <- sapply(x, function(t){ 
  fs <- fs(t)
  g <- g(t)
  S_val <- (fs/g - 1)
  return(S_val/(norm_S^2))
})
S0_bkg_vec <- sapply(y, function(t){
  fs <- fs(t)
  g <- g(t)
  S_val <- (fs/g - 1)
  return(S_val/(norm_S^2))
})
theta0_hat <- mean(S0_phys_vec); delta0_hat <- mean(S0_bkg_vec)
eta_hat_unif <- (theta0_hat-delta0_hat)/(1-delta0_hat) # estimate of eta
test_num <- sqrt(m*n)*eta_hat_unif # numerator of the test statistic
# calculating the denominator of the test statistic
d_theta_T <- 1/(1-delta0_hat)
d_delta_T <- (theta0_hat - 1)/((1-delta0_hat)^2)
var_S0_F_hat <- mean(S0_phys_vec^2) - (theta0_hat^2)
var_S0_Fb_hat <- mean(S0_bkg_vec^2) - (delta0_hat^2)
denom1 <- m*(d_theta_T^2)*var_S0_F_hat
denom2 <- n*(d_delta_T^2)*var_S0_Fb_hat
test_denom <- sqrt(denom1 + denom2) # denominator of the test statistic
test_stat <- test_num/test_denom # test statistic
p_val_unif <- pnorm(test_stat, lower.tail = FALSE) # p-value
sig_hat_unif <- test_denom/(sqrt(m+n))
std_err_unif <- sig_hat_unif*sqrt((m+n)/(m*n)) # std error
ci_95_unif <- eta_hat_unif + c(-1,1)*qnorm(0.975)*std_err_unif # 95% CI

########### Table of results ###################################################
df_res <- data.frame('Proposal_bkg' = c('Exponential', 'Gaussian-tail', 'Uniform'),
                     'eta_hat' = c(eta_hat_exp, eta_hat_gt, eta_hat_unif),
                     'CI_Lower' = c(ci_95_exp[1], ci_95_gt[1], ci_95_unif[1]),
                     'CI_Upper' = c(ci_95_exp[2], ci_95_gt[2], ci_95_unif[2]),
                     'p_val' = c(p_val_exp, p_val_gt, p_val_unif))
caption <- paste0("Signal Search Results")
kable(df_res, format = "simple", digits = 10,
      booktabs = TRUE, escape = FALSE,
      caption = caption)
