# devtools::install_github("hyunjimoon/SBC")
library(SBC)
library(LearnBayes)
library(logPoolR)
source("one_right_many_wrong_beta_aux_new.r")

cvr <- 0.2
cvw <- 0.2
pars <- get_parameter_vectors_mean_cv(cv_correct = cvr,
                                      cv_wrong = cvw)
###

# set.seed(666)
compute <- TRUE
est_X <- 1/10
n <- 100

load(paste0("saved_data/BetaDataSets_Dirichlet_prior_X=", est_X,
            "_n=", n, ".RData"))

load(
  paste0("saved_data/ResultsBetaDirichlet_prior_X=", est_X, "_n=", n, ".RData")
)

###############
get_marginal_liks <- function(dt){
  K <- dt$K
  
  marginal.likelihoods <- densities.true <- rep(NA, K)
  
  for (k in 1:K){
    marginal.likelihoods[k] <- logPoolR:::ml_beta(yi = dt$Y,
                                       ni = dt$N,
                                       a = dt$a[k],
                                       b = dt$b[k])
  }
  alphas.bma <- marginal.likelihoods/sum(marginal.likelihoods) 
  return(alphas.bma)
}

All.alphas <- results_Dirichlet$stats[grep("alpha",
                                           results_Dirichlet$stats$parameter), ]
SSAlphas <- split(All.alphas, All.alphas$dataset_id)

SSAlphas.with.ranks <- lapply(1:length(SSAlphas), function(i){
  out <- SSAlphas[[i]]
  out$bma <- get_marginal_liks(datasets_Dirichlet$generated[[i]])
  out$rbma <- rank(out$bma)
  out$rtrue <- rank(out$simulated_value)
  out$rmean <- rank(out$mean)
  out$rmedian <- rank(out$median)
  return(out)
})
agree <- do.call(rbind, lapply(1:length(SSAlphas), function(i){
  apply(SSAlphas.with.ranks[[i]][, c("rmean", "rmedian", "rbma")], 2,
        function(x){
          identical(x, SSAlphas.with.ranks[[i]]$rtrue, )
        })
}))
colMeans(agree)

agree.biggest <- do.call(rbind, lapply(1:length(SSAlphas), function(i){
  apply(SSAlphas.with.ranks[[i]][, c("rmean", "rmedian", "rbma")], 2,
        function(x){
          identical(x[3], SSAlphas.with.ranks[[i]]$rtrue[3], )
        })
}))
colMeans(agree.biggest)

mse <- do.call(rbind, lapply(1:length(SSAlphas), function(i){
  apply(SSAlphas.with.ranks[[i]][, c("rmean", "rmedian", "rbma")], 2,
        function(x){
          (x-SSAlphas.with.ranks[[i]]$rtrue)^2
        })
}))
colMeans(mse)

j <- 213
SSAlphas.with.ranks[[j]][, c("parameter", "simulated_value",
                             "mean", "median", "bma")]
SSAlphas.with.ranks[[j]][, c("rtrue", "rmean", "rmedian",
                             "rbma")]
###############

dplyr::filter(results_Dirichlet$stats, parameter == "alpha[1]")
dplyr::filter(results_Dirichlet$stats, parameter == "alpha[2]")
dplyr::filter(results_Dirichlet$stats, parameter == "alpha[3]")
dplyr::filter(results_Dirichlet$stats, parameter == "theta")


plot(simulated_value ~mean, 
     data = dplyr::filter(results_Dirichlet$stats, parameter == "alpha[1]"))

plot(simulated_value ~mean, 
     data = dplyr::filter(results_Dirichlet$stats, parameter == "astar"))

plot(simulated_value ~mean, 
     data = dplyr::filter(results_Dirichlet$stats, parameter == "bstar"))

pcdf <- SBC::plot_ecdf_diff(results_Dirichlet)
phist <- SBC::plot_rank_hist(results_Dirichlet)

par_names <- c(
  `alpha[1]` = "alpha[0]",
  `alpha[2]` = "alpha[1]",
  `alpha[3]` = "alpha[2]",
  `theta` = "theta",
  `astar` = 'a^"*"',
  `bstar` = 'b^"*"'
)

cdfplot <- pcdf  +
  theme_bw(base_size = 20) +
  facet_wrap(parameter ~ .,
             labeller = as_labeller(par_names, label_parsed))
cdfplot

ggsave(
  filename = paste0("../../plots/BetaDirichlet_prior_X=",
                    est_X, "_n=", n, "_ecdf.pdf") ,
  plot = cdfplot,
  width = 12, height = 6, units = "in",
  device = "pdf"
)



histplot <- phist +
  theme_bw(base_size = 20) +
  facet_wrap(parameter ~ .,
             labeller = as_labeller(par_names, label_parsed))
histplot

ggsave(
  filename = paste0("../../plots/BetaDirichlet_prior_X=",
                    est_X, "_n=", n, "_hist.pdf") ,
  plot = histplot,
  width = 12, height = 6, units = "in",
  device = "pdf"
)
