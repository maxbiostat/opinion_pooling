# Leo Bastos & Luiz Max Carvalho (2019)
# This example was taken from Malta et al. (2010)

source("pooling_aux.r")

meta <- read.csv("../data/meta_analysis_Malta_2010.csv")

meta$SampleSize

K <- nrow(meta)
av <- meta$HIV + 1
bv <- meta$SampleSize - meta$HIV + 1

mv <- meta$HIV/meta$SampleSize
vv <- mv*(1-mv)/meta$SampleSize
sv <- sqrt(vv)

cbind(
  t(apply(cbind(av, bv), 1, stat_beta) ),
  t(apply(cbind(mv, sv), 1, stat_gauss) )
) 

# Individual entropies
entropies <- rep(NA, K)
for(k in 1:K) entropies[k] <- entropy_gauss(mv[k], vv[k])

entropies


############
PaperMSMGauss.tbl <- data.frame(mean.prior = rep(NA, 6), lower.prior = NA, 
                             upper.prior = NA)
rownames(PaperMSMGauss.tbl) <- c("equal_weights", "maximum_entropy", "minimum_KL",
                             "hierarchical_Dirichlet", "hierarchical_LogisticNormal", "Sample_size")

AlphasMSMGauss.tbl <- data.frame(matrix(NA, nrow = 3, ncol = length(av)))
rownames(AlphasMSMGauss.tbl) <- c("maximum_entropy", "minimum_KL", "Sample_size")
colnames(AlphasMSMGauss.tbl) <- paste("alpha_", 0:(K-1), sep = "")


library(ggplot2)

phi.grid <- seq(0, 1, length.out = 1000)
study.densities <- vector(K, mode = "list")
for(k in 1:K){
  study.densities[[k]] <- data.frame(phi = phi.grid,
                                      dens = dnorm(phi.grid, mean = mv[k], sd = sv[k] ),
                                      study = paste("study_", k-1, sep = ""))
  
}
study.densities.df <- do.call(rbind, study.densities)
study.densities.df$distribution <- "Gaussian"
write.csv(study.densities.df, file =  "../data/output/MSM_Gaussian_expert_densities.csv", row.names = FALSE)

study_priors <- ggplot(study.densities.df, aes(x = phi, y = dens,
                                                 linetype = study, colour = study)) + 
  geom_line(size = 2) +
  scale_x_continuous(expression(phi), expand = c(0, 0), limits = c(0, .4)) +
  scale_y_continuous(expression(f[i](phi)), expand = c(0, 0)) +
  theme_bw(base_size = 20)

study_priors
ggsave(study_priors, filename = "../plots/study_densities_MSMGaussian.pdf")

###### Equal weights

alphaEqual <- rep(1/K, K)

ab.Equal.star <- pool_par_gauss(alphaEqual, mv, vv)
# Prior
(PaperMSMGauss.tbl[1, 1:3] <- stat_gauss(ab.Equal.star))

####### Maximum entropy

## WARNING: For the Gaussian case, we do not need to optimise, if you don't believe the maths, just run the code below
# N <- 1000 ## could increase to, say, 10000 in order to make sure, but it's fine
# ent.many.startingPoints <- matrix(rnorm(n = (K-1)*N, mean = 0, sd = 100), ncol = K-1, nrow = N)
# many.ents <- lapply(1:N, function(i) {
#   optim(ent.many.startingPoints[i, ], optentgauss_inv, mp = mv, vp = vv)
# })
# optimised.ents <- unlist(lapply(many.ents, function(x) x$value))
# 
# hist(optimised.ents)
# abline(v = optimised.ents[which.min(optimised.ents)], lty = 2, lwd = 2)
# 
# alphaMaxEnt.opt <- alpha_01(many.ents[[which.min(optimised.ents)]]$par)

## Maximum entropy "analytical" solution,
alphaMaxEnt.opt  <- rep(0, K)
alphaMaxEnt.opt[which.max(vv)] <- 1
round(alphaMaxEnt.opt, 2)

( AlphasMSMGauss.tbl[1, ] <- alphaMaxEnt.opt )

ab.MaxEnt.star <- pool_par_gauss(alphaMaxEnt.opt, mv, vv)

# Prior
(PaperMSMGauss.tbl[2, 1:3] <- stat_gauss(ab.MaxEnt.star))


####### Minimum KL

N <- 1000 ## could increase to, say, 10000 in order to make sure, but it's fine
kl.many.startingPoints <- matrix(rnorm(n = (K-1)*N, mean = 0, sd = 100), ncol = K-1, nrow = N)
many.kls <- lapply(1:N, function(i) {
  optim(kl.many.startingPoints[i, ], optklgauss_inv, mp = mv, vp = vv, type = "fp")
})
optimised.kls <- unlist(lapply(many.kls, function(x) x$value))

hist(optimised.kls)
abline(v = optimised.kls[which.min(optimised.kls)], lty = 2, lwd = 2)

alphaKL.opt <- alpha_01(many.kls[[which.min(optimised.kls)]]$par)

round(AlphasMSMGauss.tbl[2, ] <- alphaKL.opt, 2)

ab.KL.star <- pool_par_gauss(alphaKL.opt, mv, vv)

# Prior
(PaperMSMGauss.tbl[3, 1:3] <- stat_gauss(ab.KL.star))
  
####### Hierarchical priors
require("LearnBayes")

M <- 100000
X <- c(1, 1, 1, 1, 1, 1)/10
alpha.MC.dirichlet <- rdirichlet(M, X)
alpha.MC.logisticNormal <- rlogisticnorm(N = M,
                              m = digamma(X)-digamma(X[K]),
                              Sigma = constructSigma(X))

apply(alpha.MC.dirichlet, 2, mean)
apply(alpha.MC.logisticNormal, 2, mean)

apply(alpha.MC.dirichlet, 2, sd)
apply(alpha.MC.logisticNormal, 2, sd)

gauss.par.dirichlet <- apply(alpha.MC.dirichlet, 1, function(w) pool_par_gauss(w, mv, vv))
gauss.par.logisticNormal <- apply(alpha.MC.logisticNormal, 1, function(w) pool_par_gauss(w, mv, vv))

phi.par.dirichlet <- apply(gauss.par.dirichlet, 2, function(x) rnorm(1, x[1], x[2]))
phi.par.logisticNormal <- apply(gauss.par.logisticNormal, 2, function(x) rnorm(1, x[1], x[2]))
# Prior
PaperMSMGauss.tbl[4, 1] <- mean(phi.par.dirichlet)
PaperMSMGauss.tbl[4, 2:3] <- quantile(phi.par.dirichlet, c(.025, .975))

PaperMSMGauss.tbl[5, 1] <- mean(phi.par.logisticNormal)
PaperMSMGauss.tbl[5, 2:3] <- quantile(phi.par.logisticNormal, c(.025, .975))


####### Using sample sizes

alphas.sampleSize <- meta$SampleSize/sum(meta$SampleSize)

( AlphasMSMGauss.tbl[3, ] <- alphas.sampleSize )

ab.sampleSize <- pool_par_gauss(alphas.sampleSize, mv, vv)

# Prior
( PaperMSMGauss.tbl[6, 1:3] <- stat_gauss(ab.sampleSize) )


#### Finally, tables!

round(PaperMSMGauss.tbl, 3)
round(AlphasMSMGauss.tbl, 3)

round(PaperMSMGauss.tbl, 2)
round(AlphasMSMGauss.tbl, 2)
write.csv(round(PaperMSMGauss.tbl, 3), file = "../data/output/MSM_Gaussian_stat.csv", row.names = TRUE)
write.csv(round(AlphasMSMGauss.tbl, 3), file = "../data/output/MSM_Gaussian_weights.csv", row.names = TRUE)

####### Plotting

posterior_studies <- data.frame(
  alpha = as.numeric(c(AlphasMSMGauss.tbl[1, ], AlphasMSMGauss.tbl[2, ], AlphasMSMGauss.tbl[3, ])),
  lwr = rep(NA, 18),
  upr = rep(NA, 18),
  study = rep(paste("study_", 0:(K-1), sep = ""), 3),
  method = rep(c("maximum_entropy", "minimum_KL", "Sample_size"), each = K)
)

####
radar_alphas <- ggplot(data = posterior_studies,
       aes(x = study, y = alpha, group = method, colour = method, fill = method)) +
  geom_point() +
  geom_polygon(alpha = 0.4) +
  theme_bw(base_size = 16) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                     breaks = number_ticks(10)) + 
  coord_radar() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )
radar_alphas

ggsave(plot = radar_alphas, filename = "../plots/alphas_radar_MSMGaussian.pdf")


#############
# Now  let's look at marginal likelihoods for the pooled priors

pars <- list(equal_weights = ab.Equal.star,
             maximum_entropy = ab.MaxEnt.star,
             minimum_KL = ab.KL.star,
             sample_Size = ab.sampleSize)
pars

apply(AlphasMSMGauss.tbl, 1, get_ratio)

J <- length(pars)
posterior.densities.list <- vector(J, mode ="list")
for (j in 1:J){
  posterior.densities.list[[j]] <- data.frame(
    phi = phi.grid,
    dens = dnorm(phi.grid, mean = pars[[j]][1], sd = pars[[j]][2]),
    method = names(pars)[j]
  )
}

posterior.densities.df <- do.call(rbind, posterior.densities.list)
posterior.densities.df$distribution <- "Gaussian"

write.csv(posterior.densities.df, "../data/output/MSM_Gaussian_densities.csv", row.names = FALSE)

method_posteriors <- ggplot(posterior.densities.df, aes(x = phi, y = dens,
                                                        linetype = method, colour = method)) + 
  geom_line(size = 2) +
  scale_x_continuous(expression(phi), expand = c(0, 0), limits = c(0, .4)) +
  scale_y_continuous(expression(pi(phi)), expand = c(0, 0)) +
  geom_vline(xintercept = sum(meta$HIV)/sum(meta$SampleSize), linetype = "dashed", size = 1.2) +
  theme_bw(base_size = 16)

method_posteriors
ggsave(method_posteriors, filename = "../plots/method_posterior_densities_MSMGaussian.pdf")
