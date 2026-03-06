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
eta_true <- 0.02

eps <- 1e-3; beta <- 4

# signal density
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)
Fs <- function(x) ptrunc(x, a = l, b = u, spec = 'norm', mean = mean_sig, sd = sd_sig)

# Figuring out (mu_s -d, mu_s + d) that covers an area of 1-eps:
find_d <- function(d)
{
  pl <- Fs(mean_sig-d)
  pu <- Fs(mean_sig+d)
  return(pu-pl-1+eps)
}
sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))
r <- sol$root
M_lower <- mean_sig - r; M_upper <- mean_sig + r

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb_true <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)
f_true <- function(x) eta_true*fs(x)+(1-eta_true)*fb_true(x)
qb <- function(x){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta0)
}

mean1_in_gb <- 0.5*M_lower + 0.5*mean_sig; sd_in_gb <- 4*sd_sig
mean2_in_gb <- 0.4*M_upper + 0.6*mean_sig

gb <- function(x, lambda) {
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    a = l, b = u,
                    spec = 'norm')
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     a = l, b = u,
                     spec = 'norm')
  qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
               scale = l, shape = beta)
  
  return(lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb)
}
op <- par(no.readonly = TRUE)
par(mar = c(5, 5, 4, 2))
curve(f_true, from = 1.1, to = 1.6, col = 'orange', lwd = 5,
      ylab = 'Density', xlab = 'x', lty = 2,
      cex.axis = 2, cex.lab = 2)
curve(fb_true, add = TRUE, col = 'black', lwd = 4, lty = 1,
      ylab = '', xlab = 'x')
curve(gb(x, 0), add = TRUE, col = 'blue', lwd = 4, lty = 2,
      ylab = '', xlab = 'x')
curve(gb(x, 0.03), add = TRUE, col = 'purple', lwd = 4, lty = 4,
      ylab = '', xlab = 'x')

legend(x = 1.26, y = 2.6, legend = TeX('$g_{\\beta}(x; \\lambda = 0.03)$'), 
       col = 'purple', lwd = 4, lty = 4,
       bty = 'n', cex = 2,
       seg.len = 2.2)
legend(x = 1.26, y = 2.4, legend = TeX('$q(x) \\propto x^{-5}$'), 
       col = 'blue', lwd = 4, lty = 2,
       bty = 'n', cex = 2,
       seg.len = 2.2)
legend(x = 1.26, y = 2.2, legend = 'f(x)', 
       col = 'orange', lwd = 4, lty = 5,
       bty = 'n', cex = 2,
       seg.len = 2.2)
legend(x = 1.26, y = 2, legend = TeX('$f_b(x)$'), 
       col = 'black', lwd = 4, lty = 1,
       bty = 'n', cex = 2,
       seg.len = 2.2)
abline(v = c(M_lower, M_upper), col = ggplot2::alpha('grey', 0.5), lwd = 4, lty = 2)
rect(
  xleft   = M_lower,
  xright  = M_upper,
  ybottom = 0,
  ytop    = 3,
  col     = rgb(0, 1, 0, 0.15),
  border  = NA
)
par(op)