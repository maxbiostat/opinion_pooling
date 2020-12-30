# This example was taken from Malta et al. (2010)

source("pooling_aux.r")

meta <- read.csv("../data/meta_analysis_Malta_2010.csv")

meta$SampleSize

av <- meta$HIV + 1
bv <- meta$SampleSize - meta$HIV + 1
K <- nrow(meta)

beta_mean(av, bv)
beta_mode(av, bv)
beta_sd(av, bv)
beta_sd(av, bv)^2
data.frame(meta$Study, round ( t( apply(cbind(av, bv), 1, stat_beta)  ), 3))

# Individual entropies
entropies <- rep(NA, K)
for(k in 1:K) entropies[k] <- entropy_beta(av[k], bv[k])

entropies

## Entropy surface (for a future dominance analysis) 

library(fields)
ES <- entropy_surface_beta(av, bv)
export <- TRUE
if(export){
  pdf("../plots/entropy_surface_MSMBetaExample.pdf")  
}
image.plot(ES$as, ES$bs, ES$M,
           xlab = expression(a), ylab = expression(b), horizontal = TRUE,
           cex.lab = 1.5, cex.axis = 1.5, axis.args = list(font = 2),
           legend.cex = 1.5,
           legend.lab = expression(H[pi]), main = "Entropy Beta distribution", font = 2)
if(export){
  dev.off()  
}

############
PaperMSMBeta.tbl <- data.frame(mean.prior = rep(NA, 6), lower.prior = NA, 
                             upper.prior = NA)
rownames(PaperMSMBeta.tbl) <- c("equal_weights", "maximum_entropy", "minimum_KL",
                             "hierarchical_Dirichlet", "hierarchical_LogisticNormal", "Sample_size")

AlphasMSMBeta.tbl <- data.frame(matrix(NA, nrow = 3, ncol = length(av)))
rownames(AlphasMSMBeta.tbl) <- c("maximum_entropy", "minimum_KL", "Sample_size")
colnames(AlphasMSMBeta.tbl) <- paste("alpha_", 0:(K-1), sep = "")


library(ggplot2)

phi.grid <- seq(0, 1, length.out = 1000)
study.densities <- vector(K, mode = "list")
for(k in 1:K){
  study.densities[[k]] <- data.frame(phi = phi.grid,
                                      dens = dbeta(phi.grid, shape1 = av[k], shape2 = bv[k]),
                                      study = paste("study_", k-1, sep = ""))
  
}
study.densities.df <- do.call(rbind, study.densities)
study.densities.df$distribution <- "Beta"
write.csv(study.densities.df, file =  "../data/output/MSM_Beta_expert_densities.csv", row.names = FALSE)

study_priors <- ggplot(study.densities.df, aes(x = phi, y = dens,
                                                 linetype = study, colour = study)) + 
  geom_line(size = 2) +
  # scale_linetype_manual(values = c("twodash", "dotted", "longdash", "solid"))+
  scale_x_continuous(expression(phi), expand = c(0, 0), limits = c(0, .4)) +
  scale_y_continuous(expression(f[i](phi)), expand = c(0, 0)) +
  theme_bw(base_size = 20)

study_priors
ggsave(study_priors, filename = "../plots/study_densities_MSMBeta.pdf")

###### Equal weights

alphaEqual <- rep(1/K, K)

ab.Equal.star <- pool_par(alphaEqual, av, bv)
# Prior
(PaperMSMBeta.tbl[1, 1:3] <- stat_beta(ab.Equal.star))

####### Maximum entropy

N <- 1000 ## could increase to, say, 10000 in order to make sure, but it's fine
ent.many.startingPoints <- matrix(rnorm(n = (K-1)*N, mean = 0, sd = 100), ncol = K-1, nrow = N)
many.ents <- lapply(1:N, function(i) {
  optim(ent.many.startingPoints[i, ], optentbeta_inv, ap = av, bp = bv)
})
optimised.ents <- unlist(lapply(many.ents, function(x) x$value))

hist(optimised.ents)
abline(v = optimised.ents[which.min(optimised.ents)], lty = 2, lwd = 2)

alphaMaxEnt.opt <- alpha_01(many.ents[[which.min(optimised.ents)]]$par)

round(alphaMaxEnt.opt, 2)

( AlphasMSMBeta.tbl[1, ] <- alphaMaxEnt.opt )

ab.MaxEnt.star <- pool_par(alphaMaxEnt.opt, av, bv)

# Prior
(PaperMSMBeta.tbl[2, 1:3] <- stat_beta(ab.MaxEnt.star))


####### Minimum KL

N <- 1000 ## could increase to, say, 10000 in order to make sure, but it's fine
kl.many.startingPoints <- matrix(rnorm(n = (K-1)*N, mean = 0, sd = 100), ncol = K-1, nrow = N)
many.kls <- lapply(1:N, function(i) {
  optim(kl.many.startingPoints[i, ], optklbeta_inv, ap = av, bp = bv, type = "fp")
})
optimised.kls <- unlist(lapply(many.kls, function(x) x$value))

hist(optimised.kls)
abline(v = optimised.kls[which.min(optimised.kls)], lty = 2, lwd = 2)

alphaKL.opt <- alpha_01(many.kls[[which.min(optimised.kls)]]$par)

round(AlphasMSMBeta.tbl[2, ] <- alphaKL.opt, 2)

ab.KL.star <- pool_par(alphaKL.opt, av, bv)

# Prior
(PaperMSMBeta.tbl[3, 1:3] <- stat_beta(ab.KL.star))
  
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

beta.par.dirichlet <- alpha.MC.dirichlet %*% cbind(av, bv)
beta.par.logisticNormal <- alpha.MC.logisticNormal %*% cbind(av, bv)

phi.par.dirichlet <- apply(beta.par.dirichlet, 1, function(x) rbeta(1, x[1], x[2]))
phi.par.logisticNormal <- apply(beta.par.logisticNormal, 1, function(x) rbeta(1, x[1], x[2]))
# Prior
PaperMSMBeta.tbl[4, 1] <- mean(phi.par.dirichlet)
PaperMSMBeta.tbl[4, 2:3] <- quantile(phi.par.dirichlet, c(.025, .975))

PaperMSMBeta.tbl[5, 1] <- mean(phi.par.logisticNormal)
PaperMSMBeta.tbl[5, 2:3] <- quantile(phi.par.logisticNormal, c(.025, .975))



####### Using sample sizes

alphas.sampleSize <- meta$SampleSize/sum(meta$SampleSize)

(AlphasMSMBeta.tbl[3, ] <- alphas.sampleSize)

ab.sampleSize <- pool_par(alphas.sampleSize, av, bv)

# Prior
(PaperMSMBeta.tbl[6, 1:3] <- stat_beta(ab.sampleSize))


#### Finally, tables!

round(PaperMSMBeta.tbl, 3)
round(AlphasMSMBeta.tbl, 3)

round(PaperMSMBeta.tbl, 2)
round(AlphasMSMBeta.tbl, 2)

write.csv(round(PaperMSMBeta.tbl, 3), file = "../data/output/MSM_Beta_stat.csv", row.names = TRUE)
write.csv(round(AlphasMSMBeta.tbl, 3), file = "../data/output/MSM_Beta_weights.csv", row.names = TRUE)

####### Plotting

posterior_studies <- data.frame(
  alpha = as.numeric(c(AlphasMSMBeta.tbl[1, ], AlphasMSMBeta.tbl[2, ], AlphasMSMBeta.tbl[3, ])),
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

ggsave(plot = radar_alphas, filename = "../plots/alphas_radar_MSMBeta.pdf")


#############
# Now  let's look at marginal likelihoods for the pooled priors

pars <- list(equal_weights = ab.Equal.star,
             maximum_entropy = ab.MaxEnt.star,
             minimum_KL = ab.KL.star,
             sample_Size = ab.sampleSize)

lapply(pars, function(x) c(beta_mean(a = x[1], b = x[2]), beta_sd(a = x[1], b = x[2]) ) )

apply(AlphasMSMBeta.tbl, 1, get_ratio)

J <- length(pars)
method.densities.list <- vector(J, mode ="list")
for (j in 1:J){
  method.densities.list[[j]] <- data.frame(
    phi = phi.grid,
    dens = dbeta(phi.grid, shape1 = pars[[j]][1], shape2 = pars[[j]][2]),
    method = names(pars)[j]
  )
}

method.densities.df <- do.call(rbind, method.densities.list)
method.densities.df$distribution <- "Beta"

write.csv(method.densities.df, "../data/output/MSM_Beta_densities.csv", row.names = FALSE)

method_priors <- ggplot(posterior.densities.df, aes(x = phi, y = dens,
                                                        linetype = method, colour = method)) + 
  geom_line(size = 2) +
  scale_x_continuous(expression(phi), expand = c(0, 0), limits = c(0, .4)) +
  scale_y_continuous(expression(pi(phi)), expand = c(0, 0)) +
  geom_vline(xintercept = sum(meta$HIV)/sum(meta$SampleSize), linetype = "dashed", size = 1.2) +
  theme_bw(base_size = 16)

method_priors
ggsave(method_priors, filename = "../plots/method_prior_densities_MSMBeta.pdf")