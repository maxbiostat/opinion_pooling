# Just a little something for us to begin studying how to do a PSA for \pi(\theta)
# imagine an equal-parameter prior with X_i = x for a n-dimensional simplex

dirichlet_var <- function(x, n){
  X <- rep(x, n)
  N <- n*x
  return( (x*(N-x))/ (N^2*(N+1)))
}

xs <- seq(.1, 10, .1)
# Stuff for the Jeffreys' prior
which(xs == .5)
dirichlet_var(1/2, 4)
dirichlet_var(1/2, 8)

svg("../plots/Dirichlet_var.svg")
plot(xs, dirichlet_var(xs, 2), type = "l", lwd = 3, ylab = expression(Var(theta)), xlab = expression(X[i]))
lines(xs, dirichlet_var(xs, 4), type = "l", col = 2, lwd = 3)
lines(xs, dirichlet_var(xs, 5), type = "l", col = 3, lwd = 3)
lines(xs, dirichlet_var(xs, 10), type = "l", col = 4, lwd = 3)
legend(x = "top",legend = c("K=2", "K=4", "K=5", "K=10"),
       title = "Number of priors (K)", lwd = 2, col = 1:4, bty="n")
dev.off()