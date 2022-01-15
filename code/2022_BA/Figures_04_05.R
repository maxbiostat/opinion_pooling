library(ggplot2)
library(ggthemes)
library(Ternary)

source("one_right_many_wrong_beta_aux_new.r")

simus <- read.csv("saved_data/oracle_beta/oracle_Beta_estX=0.1.csv")

par_names <- c(
  `alpha[1]` = "alpha[0]",
  `alpha[2]` = "alpha[1]",
  `alpha[3]` = "alpha[2]",
  `10` = "n==10",
  `100` = "n==100",
  `10000` = "n==10000",
  `theta` = "theta",
  `astar` = 'a^"*"',
  `bstar` = 'b^"*"'
)

pp <- ggplot(data = simus[grep("alpha", simus$par),],
             aes(x = xbar, y = mean, colour = method, fill = method)) +
  geom_line() + 
  scale_x_continuous(expression(bar(y[n]))) +
  scale_y_continuous("Posterior weight") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .4) +
  facet_grid(par~n, scales = "free_y",
             labeller = as_labeller(par_names,
                                    label_parsed)) +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  theme_bw(base_size = 20)

pdf("../../plots/alpha_posterior_oracle_flexible.pdf",
    paper = "a4r",
    width = 10, height = 8)
pp
dev.off()


beta_linear_kernel <- function(x, s1, s2, w, log = FALSE){
  require(matrixStats)
  K <- length(s1)
  if(length(s2) != K || length(w) != K){
    stop("Hyperparameters and weights not the same length")
  } 
  lps <- sapply(1:K, function(i){
    log(w[i]) + dbeta(x = x, shape1 = s1[i], shape2 = s2[i], log = TRUE)
  })
  ans <- matrixStats::logSumExp(lps)
  if(!log) ans <- exp(ans)
  return(ans)
}

beta_linear_mixture <- function(x, s1, s2, w, log = FALSE){
  ans <- sapply(x,
                function(x) beta_linear_kernel(x, s1, s2, w, log))
  return(ans)
}

produce_plots <- function(sampleMean, n, cvr, cvw, est_X){
  
  fName <- paste0("saved_data/oracle_beta/oracleSimulationBeta_",
                  cvr, "_", cvw,
                  "_x_bar=", sampleMean,
                  "_n=", n,
                  "_X=", est_X[1], ".RData")
  
  load(file = fName)
  
  K <- length(hyperpars$av)
  the.colours <- RColorBrewer::brewer.pal(K, name = "Dark2")
  
  ## Figure 1: Densities
  fn1 <- paste0("../../plots/oracle_beta/densities_",
                cvr, "_", cvw,
                "_x_bar=", sampleMean,
                "_n=", n,
                "_X=", est_X[1], ".pdf")
  
  pdf(file = fn1, paper = "a4")
  
  plot_densities(hyperpars, lg = FALSE, main = "")
  
  curve(dbeta(x,
              shape1 = hpar.bma[1],
              shape2 = hpar.bma[2]),
        lwd = 3, lty = 2, col = "blue", add = TRUE)
  
  curve(beta_linear_mixture(x, s1 = hyperpars$av,
                            s2 = hyperpars$bv,
                            w = alphas.bma),
        lwd = 3, lty = 3, col = "purple", add = TRUE)
  
  legend(x = "topright", 
         legend = c("expert_0", "expert_1", "expert_2",
                    "log-linear BMA", "Linear BMA"),
         col = c(the.colours, "blue", "purple"),
         lty = c(1:K, 2, 3),
         bty = 'n')
  
  dev.off()
  
  ## Figure 2: Ternary plot, Dirichlet
  
  fn2 <- paste0("../../plots/oracle_beta/ternary_Dirichlet_",
                cvr, "_", cvw,
                "_x_bar=", sampleMean,
                "_n=", n,
                "_X=", est_X[1], ".pdf")
  
  
  ttt <- TernaryDensity(posterior.alphas.Dirichlet,
                        resolution = 30L)
  
  pdf(file = fn2, paper = "a4")
  
  par(mar = rep(0.75, 4))
  Ternary::TernaryPlot(atip = "Expert 0", 
              btip = "Expert 1",
              ctip = "Expert 2",
              axis.labels = seq(0, 1, length.out = 11),
              main = paste0("Dirichlet posterior (", est_X[1], ")"))
  Ternary::ColourTernary(TernaryDensity(posterior.alphas.Dirichlet,
                               resolution = 30L))
  Ternary::TernaryPoints(posterior.alphas.Dirichlet, col = 'red', pch = '.')
  Ternary::TernaryPoints(alphas.bma, col = 'black', pch = 17, cex = 2)
  Ternary::TernaryPoints(colMeans(posterior.alphas.Dirichlet), col = 'black',
                pch = 16, cex = 2)
  Ternary::TernaryDensityContour(posterior.alphas.Dirichlet, resolution = 30L)
  
  dev.off()
  
  ## Figure 3: Ternary plot, logistic-Normal
  
  fn3 <- paste0("../../plots/oracle_beta/ternary_logisticNormal_",
                cvr, "_", cvw,
                "_x_bar=", sampleMean,
                "_n=", n,
                "_X=", est_X[1], ".pdf")
  
  pdf(file = fn3, paper = "a4")
  
  par(mar = rep(0.75, 4))
  TernaryPlot(axis.labels = seq(0, 1, length.out = 11),
              atip = "Expert 0", 
              btip = "Expert 1",
              ctip = "Expert 2",
              main = paste0("Logistic-normal posterior (", est_X[1], ")"))
  ColourTernary(TernaryDensity(posterior.alphas.logisticNormal,
                               resolution = 30L))
  TernaryPoints(posterior.alphas.logisticNormal, col = 'red', pch = '.')
  TernaryPoints(alphas.bma, col = 'black', pch = 17, cex = 2)
  TernaryPoints(colMeans(posterior.alphas.logisticNormal), col = 'black',
                pch = 16, cex = 2)
  TernaryDensityContour(posterior.alphas.logisticNormal, resolution = 30L)
  dev.off()
}

xbs <- 1:9/10
sizes <- c(10, 100, 1E4)

grid <- expand.grid(xbar = xbs, n = sizes)
J <- nrow(grid)

lapply(1:J, function(j){
  produce_plots(sampleMean = grid[j, 1],
               n = grid[j, 2],
               cvr = .2,
               cvw = .2,
               est_X = rep(1/10, 3))})
