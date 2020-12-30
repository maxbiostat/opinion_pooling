M <- 1e6
X1 <- rnorm(M)
X2 <- rbeta(M, shape1 = 2, shape2 = 3)
X3 <- rgamma(M, shape = 1, rate = 1)

phi <- function(x) floor(x)

Y1 <- phi(X1)
Y2 <- phi(X2)
Y3 <- phi(X3)

barplot(table(Y1)/M)
barplot(table(Y2)/M)
barplot(table(Y3)/M)
