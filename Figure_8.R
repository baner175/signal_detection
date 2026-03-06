rm(list = ls())

library(latex2exp)

# WITH BKG ONLY SAMPLE ##########################################

# beta estimated:
B <- 1e5
n_seq <- c(50, 250, 5e2, 1e3, 2e3)
lambda_seq <- c(0.00, 0.01, 0.02, 0.03)

eta <- 0.

transp <- 0.8
op <- par(no.readonly = TRUE)
par(mgp = c(2.5, 0.8, 0))
if(eta == 0){
  plot(x = 500, y = 0, 
       xlim = c(50, 2e3),
       ylim = c(0.0, 0.25), type = 'n',
       ylab = 'P(Type I error)',
       xlab = 'Size of the physics data',
       cex.lab = 2, cex.axis = 2,
       xaxt = 'n',
       yaxt = 'n',
       cex.lab = 2.2, cex.axis = 2.2, cex.main = 2.2,
       main = TeX(sprintf('$\\eta = %f$', eta))
  )
  ax_marks <- seq(0.0, 0.25, 0.05)
  axis(1, at = n_seq,
       labels = as.character(n_seq),
       cex.axis = 1.8)
  axis(2, at = ax_marks,
       labels = as.character(ax_marks),
       cex.axis = 1.8)
  abline(h = ax_marks, lty = 2, lwd = 4,
         col = 'grey')
}else{
  plot(x = 500, y = 0, 
       xlim = c(50, 2e3),
       ylim = c(0,1), 
       yaxt = 'n',
       xaxt = 'n',
       type = 'n',
       ylab = 'Power',
       xlab = 'Size of the physics data',
       cex.lab = 2.2, cex.axis = 2.2, cex.main = 2.2,
       main = TeX(sprintf('$\\eta = %f$', eta))
  )
  ax_marks <- seq(0.05, 1, 0.1)
  axis(1, at = n_seq,
       labels = as.character(n_seq),
       cex.axis = 1.8)
  axis(2, at = ax_marks,
       labels = as.character(ax_marks),
       cex.axis = 1.8)
  abline(h = ax_marks, lty = 2, lwd = 4,
         col = 'grey')
}

pow_res <- c()

for(n in n_seq){
  file_name <- paste0('Numerical_example_Sec4/',
                      'WOBKG__',
                      'eta(',eta,')_',
                      'lambda(', paste0(lambda_seq, collapse = '_'),')_',
                      'B(',B,')_',
                      'n_phys(', n,').txt')
  pow <- read.table(file_name, header = TRUE)[,2]
  pow_res <- cbind(pow_res, pow)
}

delta_seq <- read.table(file_name, header = TRUE)[,3]

# lambda = 0:
lines(x = n_seq, y = pow_res[1,], lty = 2, lwd = 4,
      col = ggplot2::alpha('blue', transp))
points(x = n_seq, y = pow_res[1,], cex = 2,
       col = ggplot2::alpha('blue', transp), lwd = 4, pch = 16)


# lambda = 0.01:
lines(x = n_seq, y = pow_res[2,], lty = 6, lwd = 4,
      col = ggplot2::alpha('cyan3', transp))
points(x = n_seq, y = pow_res[2,], cex = 2,
       col = ggplot2::alpha('cyan3', transp), lwd = 4, pch = 17)

# lambda = 0.02:
lines(x = n_seq, y = pow_res[3,], lty = 5, lwd = 4,
      col = ggplot2::alpha('orange', transp))
points(x = n_seq, y = pow_res[3,], cex = 2,
       col = ggplot2::alpha('orange', transp), lwd = 4, pch = 18)

# lambda = 0.03:
lines(x = n_seq, y = pow_res[4,], lty = 4, lwd = 4,
      col = ggplot2::alpha('purple', transp))
points(x = n_seq, y = pow_res[4,], cex = 2,
       col = ggplot2::alpha('purple', transp), lwd = 4, pch = 19)

pdf("power_curve_legend_wobkg.pdf", width = 20, height = 2)
par(mar = c(0, 0, 0, 0))
plot.new()
labs <- Map(function(lam, del) {
  bquote(atop(lambda == .(sprintf("%.2f", lam)),
              (delta[beta^"*"] == .(sprintf("%.3f", del)))))
}, lambda_seq, delta_seq)
tw <- strwidth(TeX("$\\lambda = 0.00$\\n$(\\delta_{\\beta^*} = -0.000)$")) * 2
legend(
  "bottom",
  legend = labs,
  col = c("blue", "cyan3", "orange", "purple"),
  lty = c(2, 6, 5, 4),
  pch = 16:19,
  lwd = 4,
  cex = 2,
  bty = "n",
  ncol = 4,
  x.intersp = 0.5,
  seg.len = 2.5,
  text.width = tw
)

dev.off()

par(op)
