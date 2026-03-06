rm(list = ls())
library(truncdist)
library(VGAM)
library(doSNOW)
library(latex2exp)
####################################################################
# LRT simulation:
B <- 1e5; n_samp <- 5e3; eta_true <- 0; beta0 <- 4
eps_seq <- c(0, 0.005, 0.01)

theo_CDF_0 <- function(x) (0.5 + 0.5*pchisq(x, df = 1))*(x>=0)
op <- par(no.readonly = TRUE)
par(mgp = c(2.5, 0.8, 0))
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
  Fn_hat <- ecdf(df[,1])
  curve(Fn_hat, add = TRUE, lwd = 4, lty = my_lty[i], col = mycols[i])
}

legend(x = 7, y = 0.55,
       legend = c(
         TeX('$\\bar{\\chi}^2_{0 1}$'),
         TeX('$q(x)$'),
         TeX(sprintf('$\\tilde{g}(x, \\epsilon = %.4f)$', eps_seq[-1]))
       ),
       col = c('black', mycols),
       lty = c(1,my_lty),
       bty = 'n',
       lwd = 4,
       y.intersp = 1.5, cex = 2)
par(op)