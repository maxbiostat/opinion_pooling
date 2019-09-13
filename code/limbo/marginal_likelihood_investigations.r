ml.beta <- function(yi, ni, a, b){ # Equation (9) in Raftery et al (2007) 
  ( gamma(ni + 1)/{gamma(ni - yi + 1) * gamma(yi + 1)} ) *
    ( gamma(a + b)/gamma(a + b + ni)  ) *
    ( gamma(a+yi)/gamma(a) ) * ( gamma(b + ni - yi)/gamma(b)) 
}
##
# parameters
a0 <- 18.1 ; b0 <- .995
a1 <- 3.44 ; b1 <- .860 
a2 <- 8.32 ; b2 <- .924
a3 <- 1.98 ; b3 <- .848

av <- c(a0, a1, a2, a3)
bv <- c(b0, b1, b2, b3)
K <- length(av)
# fix y and n
y <- 9
n <- 10

# create the grid
pa <- seq(min(av), max(av), length.out = 50)
pb <- seq(min(bv), max(bv), length.out = 50)
grid <- expand.grid(a = pa , b = pb)


results <- apply(grid, 1, function(x) ml.beta(yi = y, ni = n, a = x[1], b = x[2]))
mat <- matrix(results, ncol = length(pa))
library(fields)
image.plot(pa, pb, mat, xlab = expression(a), ylab = expression(b))

marglikes <- rep(NA, K)
for (k in 1:K){ marglikes[k] <- ml.beta(yi = y, ni = n, a = av[k], b = bv[k]) }

marglikes
max(results)

grid[which(results == max(results)), ]

## These results indicate it is possible to obtain a pooled prior with
## higher integrated likelihood than any of the original priors.