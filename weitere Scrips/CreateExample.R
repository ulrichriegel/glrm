v <- c(101, 105, 107, 102, 115, 140, 142, 135, 133, 140, 145, 146, 180, 178, 155)
m <- c(0.1, 0.2, 0.15, 0.1, 0.05, 0.04, 0.03, 0.03, 0.02, 0.01, 0.005, 0.005, 0.006, -0.005, 0.002)
s_sq <- (m * v[10] * 0.1)^2 / v[10]
epsilon <- 0.02 * rep(1, 14)
K <- 1:10
NumberOfSimulations <- 5
Sim <- SimulateClaims(v, m, s_sq, epsilon, NumberOfSimulations, K)
