rm(list = ls())

library(truncdist)
library(VGAM)
library(doSNOW)
library(parallel)
library(foreach)
library(optparse)
library(latex2exp)
l <- 1; u <- 2

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

# beta0 <- 3.87
# or
beta0 <- 4
qb <- function(x, beta = beta0){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta)
}
gb <- function(x, lambda) lambda*fs(x) + (1-lambda)*qb(x)

delta_vs_eps <- function(eps){
  S <- function(x) fs(x)/gb(x, lambda = eps) - 1
  normS <- integrate(function(x) (S(x)^2)*gb(x, lambda = eps), l, u)$value |> sqrt()
  
  S1 <- function(x) S(x)/normS
  return(integrate(function(x) S1(x)*fb_true(x), l, u)$value)
}
delta_vs_eps <- Vectorize(delta_vs_eps)
curve(delta_vs_eps, 0, 0.02, lwd = 4)
abline(h = 0, lty = 2, col = 'red', lwd = 4)


eps_seq <- c(0, 0.005, 0.01)
delta_vals <- sapply(eps_seq, delta_vs_eps)

op <- par(no.readonly = TRUE)
par(mgp = c(2.5, 0.8, 0))
curve(fb_true, from = l, to = u, 
      col = 'black', lwd = 4, lty = 1, ylab = 'Density',
      cex.lab = 2, cex.axis = 2)
my_cols <- c('blue', 'orange', 'red')
my_lty <- c(2,6,4)

# curve(gb, l, u, col = 'red')
for(i in 1:length(eps_seq)){
  curve(gb(x, lambda = eps_seq[i]), col = my_cols[i],
        lwd = 4, lty = my_lty[i],
        add = TRUE)
}
legend(x = 1.2, y = 4,
       legend = c(
         TeX('$f_b(x);$'),
         TeX('$q(x);$'),
         TeX(sprintf('$\\tilde{g}(x, \\epsilon = %.4f);$', eps_seq[-1]))
       ),
       col = c('black', my_cols),
       lty = c(1, my_lty),
       bty = 'n',
       lwd = 4,
       y.intersp = 1.5, cex = 2)
legend(x = 1.62, y = 4,
       legend = TeX(sprintf('$\\tilde{\\delta} = %.5f $', c(0, round(delta_vals, 3)))),
       bty = 'n',
       y.intersp = 1.5, cex = 2)
par(op)
