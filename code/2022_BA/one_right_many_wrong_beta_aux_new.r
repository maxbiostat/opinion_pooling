get_parameter_vectors_mean_cv <- function(cv_correct, cv_wrong = .1){
  pars_0 <- logPoolR:::elicit_beta_mean_cv(m0 = .1, cv = cv_wrong)
  pars_1 <- logPoolR:::elicit_beta_mean_cv(m0 = .5, cv = cv_correct) # big shot
  pars_2 <- logPoolR:::elicit_beta_mean_cv(m0 = .9, cv = cv_wrong)
  
  av <- unlist( c(pars_0[1], pars_1[1], pars_2[1]) )
  bv <- unlist( c(pars_0[2], pars_1[2], pars_2[2]) )
  return(list(
    av = as.numeric(av),
    bv = as.numeric(bv)
  ))
}
#
plot_densities <- function(pars, lg = TRUE, add = FALSE, main = ""){
  av <- pars$av
  bv <- pars$bv
  K <- length(av)
  the.colours <- RColorBrewer::brewer.pal(K, name = "Dark2")
  
  curve(dbeta(x, av[1], bv[1]),
        ylab = "Density",
        xlab = expression(theta),
        main = main,
        col = the.colours[1],
        lwd = 2, add = add)
  for(k in 2:K){
    curve(dbeta(x, av[k], bv[k]),
          lwd = 2, col = the.colours[k],
          lty = k, add = TRUE)
  }
  if(lg){
    if(add){
      legend(x = "topright",
             col = c(the.colours, "grey50"),
             lwd = 2, lty = 1:K, bty = 'n',
             cex = 0.85,
             legend = c(paste0("Expert_", (1:K)-1),
                        "Induced pooled prior") )
    }else{
      legend(x = "topright", col = the.colours,
             lwd = 2, lty = 1:K, bty = 'n',
             cex = 0.85,
             legend = paste0("Expert_", (1:K)-1) )
    }
  }
}
