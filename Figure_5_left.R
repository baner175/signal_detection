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
gb <- function(x, beta){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta)
}

neg_E_fb_log_gb <- function(beta)
{
  return(
    -integrate(function(x){
      gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta)
      fb <- dtrunc(x, a = l, b = u, spec = 'gamma',
                   rate = bkg_rate, shape = bkg_shape)
      return(log(gb)*fb)
    }, l, u)$value 
  )
}

beta_star <- nlminb(0.01, neg_E_fb_log_gb)$par

op <- par(no.readonly = TRUE)
par(mgp = c(2.5, 0.8, 0))

curve(fb_true, from = l, to = u, lwd = 4, lty = 1, col = 'black',
      xlab = 'x', ylab = 'Density', cex.axis = 2, cex.lab = 2)
curve(gb(x, beta = beta_star), col = 'blue', lwd = 4, lty = 2, add = TRUE)
curve(gb(x, beta = 2), col = 'brown', lwd = 4, lty = 4, add = TRUE)
curve(dunif(x, l, u), col = 'red', lwd = 4, lty = 3, add = TRUE)


legend('top',
       legend = c(
         TeX('$f_b(x)$'),
         TeX('$g_{\\beta^*}(x) \\propto x^{-(\\beta^*+1)}$'),
         TeX('$g(x) \\propto x^{-3}$'),
         TeX('$g(x) \\propto 1$')
       ),
       bty = 'n', lty = c(1,2,4,3), 
       col = c('black', 'blue', 'brown', 'red'), lwd = 4,
       cex = 2, y.intersp = 1.3)
par(op)
