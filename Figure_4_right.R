rm(list = ls())
library(truncdist)
library(VGAM)
library(latex2exp)

l <- 1; u <- 2 # search region

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02

# signal density
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)

# background model for the spurious signal method
beta0 <- 4
q <- function(x, beta = beta0){ # misspecified benchmark model for background
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta)
} 
g <- function(x, lambda) lambda*fs(x) + (1-lambda)*q(x) # misspecified background with spurious signal component

# function to calculate the compensator for each value of epsilon, i.e., size of spurious signal component
delta_vs_eps <- function(eps){
  S <- function(x) fs(x)/g(x, lambda = eps) - 1
  normS <- integrate(function(x) (S(x)^2)*g(x, lambda = eps), l, u)$value |> sqrt()
  
  S_dag <- function(x) S(x)/normS
  return(integrate(function(x) S_dag(x)*fb(x), l, u)$value)
}
delta_vs_eps <- Vectorize(delta_vs_eps)
eps_seq <- c(0, 0.005, 0.01) # candidate values for epsilon
delta_vals <- sapply(eps_seq, delta_vs_eps)
####################################################################
# LRT simulation parameters to load simulation results
B <- 1e5; n_samp <- 5e3; eta_true <- 0; beta0 <- 4
eps_seq <- c(0, 0.005, 0.01)

theo_CDF_0 <- function(x) (0.5 + 0.5*pchisq(x, df = 1))*(x>=0) # CDF of half chi square distribution

# generating the plot
op <- par(no.readonly = TRUE)
par(mgp = c(2, 0.8, 0))
curve(theo_CDF_0, from = 0, to = 20, 
      ylab = TeX('$P(X \\leq x)$'), ylim = c(0,1),
      lwd = 4, lty = 1, col = 'black',
      cex.lab = 2, cex.axis = 2)
abline(h = 0.5, v = 0, col = 'grey', lty = 2, lwd = 4)

mycols <- c('blue', 'orange', 'red')
my_lty <- c(2,6,4)
for(i in 1:length(eps_seq))
{
  file_name <- paste0('Numerical_example_Sec2.1/', 
                      'LRT_B(', B,
                      ')_beta0(', beta0,
                      ')_n_samp(', n_samp,
                      ')_eta(', eta_true,
                      ')_lambda(', eps_seq[i],').csv')
  df <- read.csv(file_name, header = TRUE)
  Fn_hat <- ecdf(df[,1]) # computing empirical CDF
  curve(Fn_hat, add = TRUE, lwd = 4, lty = my_lty[i], col = mycols[i])
}

legend(x = 2, y = 0.55,
       legend = c(
         TeX('$\\bar{\\chi}^2_{0 1}$'),
         TeX('$\\tilde{g}(\\epsilon = 0) \\equiv q; $'),
         TeX(sprintf('$\\tilde{g}(\\epsilon = %.4f); $', eps_seq[-1]))
       ),
       col = c('black', mycols),
       lty = c(1,my_lty),
       bty = 'n',
       lwd = 4,
       seg.len = 2,
       y.intersp = 1.5, cex = 2)
legend(x = 11.5, y = 0.53,
       legend = TeX(sprintf('$\\tilde{\\delta} = %.5f $', c(0, round(delta_vals, 3)))),
       bty = 'n',
       y.intersp = 1.5, cex = 2)
par(op)
