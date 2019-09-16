fake <-  FALSE
if(fake){
  ## values compatible with boarding school data
  true.N <- 763
  true.gamma <- .5
  true.R0 <- 0.00218/.5 * true.N
  true.beta <- true.gamma * true.R0
  true.sigma <- .2
  
  library(deSolve)
  sir <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS <- -beta * S * I
      dI <- beta * S * I - gamma * I
      dR <- gamma * I
      
      return(list(c(dS, dI, dR)))
    })
  }
  true.s0 <- .9998
  iniI <- 1-true.s0
  init    <- c(S = true.s0, I = iniI, R = 0.0)
  parameters <- c(beta = true.beta, gamma = true.gamma)
  times      <- seq(0, 14, length.out = 15)
  sol <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters)) 
  plot(I ~time, sol)
  forfit.sol <- subset(sol, time > 0)
  N <- nrow(forfit.sol)
  noisy_I <- rlnorm(N, mean = log(forfit.sol$I), sd = true.sigma)
  iniTime <- 0
  
  hamper <- FALSE
  if(hamper){
    noisy_I <- ifelse(noisy_I >= 1, 1, noisy_I)
  }
}else{
  library(outbreaks)
  # data("zika_girardot_2015")
  # Ndata <- 102225 ## pop size for the 'zika_girardot_2015' data set
  # sol <- zika_girardot_2015
  # sol$time <- as.numeric(zika_girardot_2015$date-min(zika_girardot_2015$date))
  # sol$I <- sol$cases
  # forfit.sol <- subset(sol, time > 0)
  # noisy_I <- forfit.sol$I/Ndata
  # iniTime <- 0
  # iniI <- 1/Ndata
  #$$$$$$$$$$$$$$$$$$
  data(influenza_england_1978_school)
  Ndata <- 763
  sol <- influenza_england_1978_school
  sol$time <- as.numeric(sol$date-min(sol$date)) + 2
  sol$I <- sol$in_bed
  forfit.sol <- sol
  noisy_I <- forfit.sol$I/Ndata
  iniTime <- 0
  iniI <- 1/Ndata
  #$$$$$$$$$$$$$$$$$$$
  # data("zika_sanandres_2015")
  # Ndata <- 54513 ## pop size for the 'zika_sanandres_2015' data set
  # sol <- zika_sanandres_2015
  # sol$time <- as.numeric(zika_sanandres_2015$date-min(zika_sanandres_2015$date))
  # sol$I <- sol$cases
  #$$$$$$$$$$$$$$$$$$$
}

explore <- FALSE
if(explore){
  expose_stan_functions(SIR_code)
  qtilde_phi <- function(x, mub, sdb, mug, sdg, mur, sdr, alpha){
    ld <- qtilde_phi_lpdf(y = x, mb = mub, sb = sdb, mg = mug, sg = sdg, mr = mur, sr = sdr, alpha = alpha)
    return(exp(ld))
  }
  qtilde_phi <- Vectorize(qtilde_phi)
  
  curve(qtilde_phi(x,
                   mub = epi.data$mu_beta, sdb = epi.data$sigma_beta,
                   mug = epi.data$mu_gamma, sdg = epi.data$sigma_gamma,
                   mur = epi.data$mu_r0, sdr = epi.data$sigma_r0,
                   alpha = 0.5), 0 , 5, ylab = "Density", xlab = expression(R[0]), lwd = 2)
  curve(qtilde_phi(x,
                   mub = epi.data$mu_beta, sdb = epi.data$sigma_beta,
                   mug = epi.data$mu_gamma, sdg = epi.data$sigma_gamma,
                   mur = epi.data$mu_r0, sdr = epi.data$sigma_r0,
                   alpha = 0), 0 , 5, ylab = "Density", xlab = expression(R[0]), lwd = 2, col = 2, add = TRUE)
  curve(qtilde_phi(x,
                   mub = epi.data$mu_beta, sdb = epi.data$sigma_beta,
                   mug = epi.data$mu_gamma, sdg = epi.data$sigma_gamma,
                   mur = epi.data$mu_r0, sdr = epi.data$sigma_r0,
                   alpha = 1), 0 , 5, ylab = "Density", xlab = expression(R[0]), lwd = 2, col = 3, add = TRUE)  
}