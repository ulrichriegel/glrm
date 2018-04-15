rm(list = ls())



v <- Volumes
C <- Cumulative_Triangle
S <- Decumulated_Triangle
v <- v/1000000
C <- C/1000000
S <- S/1000000
n <- ncol(C)
S <- matrix(0, ncol = n, nrow = n)
for (k in n:2) {
  S[1:(n-k+1),k] <- C[1:(n-k+1),k]-C[1:(n-k+1),k-1]
}
S[,1] <- C[,1]


#Optionen:
OptionRequirePositiveWeights <- F
#OptionRequirePositiveWeights <- T
eps <- 0.05
c <- rep(1,n-1)
epsilon <- c * eps
epsilon_sq <- epsilon^2
K <- 1:(n-1)
K <- 1:2
K <- c(1,3)

lower_bound_eps <- 0.001
upper_bound_eps <- 0.1
lower_bound_cv_S <- rep(0.001, n) # lower_bound_cv_S[k] = lower bound for cv(S_{i,k})
upper_bound_cv_S <- rep(100, n) # lower_bound_cv_S[k] = lower bound for cv(S_{i,k})



# apply additive method
ResultsAdditive <- LR_Method(S,v)
m <- ResultsAdditive$m_hat
s_squared <- ResultsAdditive$s_hat_sq
UltimatesAdditiveMethod <- apply(ResultsAdditive$Predictions, 1, 'sum')
W <- ResultsAdditive$W_opt

# alphas for calculation of MSE (coefficients for predictions)
# Alpha[[i]] = Reserve of Accident Year n
# Alpha[[n+1]] = Total Reserve
Alpha <- list()
Temp <- matrix(NA, ncol = n, nrow = n)
for (i in 1:n) {
  Alpha[[i]] <- (row(Temp) == i & row(Temp) + col(Temp) > n+1 ) + 0 # add 0 to convert logical to numeric
}
Alpha$TotalReserve <- (row(Temp) + col(Temp) > n+1 ) + 0 # Alpha[[n+1]]
Temp <- NULL
NrOfAlpha <- n+1


StandardError_Additive <- sqrt( MSE_LR(S,W,Alpha[[NrOfAlpha]],rep(0.0,n-1), m, s_squared, v)$MSE)
CovMatrix_S_hat_Additive_eps0 <- MSE_LR(S,W,Alpha[[NrOfAlpha]],rep(0.0,n-1), m, s_squared, v)$CovMatrix_S_hat

# Check and reduce K if necessary (if not all incremental claims are positive)
K <- CheckK(S,K)
# Matrix J[[i]] contains the set $J_i$
J <- CreateJ(S, K, UseAllEstimators)
# Matrix containing the set $J$, N[i] = $#J_i$
J_all <- J[[1]]
N <- nrow(J[[1]])
for (i in 2:n) {
  J_all <- rbind(J_all,J[[i]])
  N <- c(N,nrow(J[[i]]))
}

#ResultsGLR <- Generalized_LR_Method(S, v, Precision = 0.0000000001, Iterations = 101, eps_start = 0.05)
#ResultsGLR <- Generalized_LR_Method(S, v, setS = c(1,2,3,4,5,6,7))
#ResultsGLR <- Generalized_LR_Method(S, v, Iterations = 50, eps_ext = 0.0315, K=K)
ResultsGLR <- Generalized_LR_Method(S, v, Iterations = 50, K=K)
#ResultsGLR <- Generalized_LR_Method(S, v, lower_bound_eps = 0.05, upper_bound_eps = 0.05, Iterations = 200)
#ResultsGLR <- Generalized_LR_Method(S, v, Iterations = 200, eps_start = 0.0315)
#ResultsGLR <- Generalized_LR_Method(S, v, setS=c(1,2,3,4,5,6,7), RequirePositiveWeights = OptionRequirePositiveWeights)

g_hat_infty <- ResultsGLR$g_hat_infty
h_hat_infty <- ResultsGLR$h_hat_infty
eps_hat_infty <- ResultsGLR$eps_hat_infty
s_hat_sq_infty <- ResultsGLR$s_hat_sq_infty
m_hat_infty <- ResultsGLR$m_hat_infty
v_hat_infty <- ResultsGLR$v_hat_infty


W <- OptimalWeights_LR_Method(v, m, s_squared, eps_hat_infty * c)
ResultsAdditiveOpt <- LR_Method(S, v, W)

MSE_LR_epshatinfty <- MSE_LR(S,W,Alpha[[NrOfAlpha]],eps_hat_infty * c, m_hat_infty, s_hat_sq_infty, v_hat_infty)
CovMatrix_S_hat_Additive_eps_hat_infty <- MSE_LR_epshatinfty$CovMatrix_S_hat
L <- MSE_LR_epshatinfty$L
seLRepshatinfty <- sqrt(MSE_LR_epshatinfty$MSE)

#MSE_GLR_temp <- MSE_GLR(S, Alpha[[NrOfAlpha]], K,  g_hat_infty, h_hat_infty, eps_hat_infty * c, m_hat_infty, s_hat_sq_infty, v_hat_infty, CovDotMatrixName = "CovdotMatrix - Bsp ASTIN Bulletin.txt")
MSE_GLR_temp <- MSE_GLR(S, Alpha[[NrOfAlpha]], K,  g_hat_infty, h_hat_infty, eps_hat_infty * c, m_hat_infty, s_hat_sq_infty, v_hat_infty, blnUseRcpp = T)

J_s <- c(J_all)
K_s <- K
Number_of_Est_s <- numeric(8)
g_s <- numeric()
h_s <- numeric()

for (i in 1:8) {
  Number_of_Est_s[i] <- nrow(J[[i]])
  g_s <- c(g_s,g_hat_infty[[i]])
  h_s <- c(h_s,h_hat_infty[[i]])
}
Index_S_hat_infty1_s <- c(7,3)
Index_S_hat_infty2_s <- c(7,3)
epsilon_s <- c*eps_hat_infty
m_s <- m_hat_infty
s_squared_s <- s_hat_sq_infty
v_s <- v_hat_infty

rslt <- Cov_S_hat_infty1_S_hat_infty2_cpp(J_s, K_s, g_s, h_s, Index_S_hat_infty1_s, Index_S_hat_infty2_s, epsilon_s, m_s, s_squared_s, v_s, n, Number_of_Est_s)
rslt2 <- Cov_S_hat_infty1_S_hat_infty2(J,  K,  g_hat_infty, h_hat_infty, Index_S_hat_infty1_s, Index_S_hat_infty2_s, epsilon_s, m_hat_infty, s_hat_sq_infty, v_hat_infty)


Index_S <- matrix(c(2,2,2,2), ncol=2, byrow = T)
Index_r <- matrix(c(2,4,2,1,4,1,1,3,2,2,4,1), ncol=3, byrow = T)

Index_S_s <- c(t(Index_S))
Index_r_s <- c(t(Index_r))

rslt_cpp <- Test_cpp(J_s, K_s, g_s, h_s, Index_S_s, Index_r_s, epsilon_s, m_s, s_squared_s, v_s, n, Number_of_Est_s)
rslt_R <- Test_R(J,  K,  g_hat_infty, h_hat_infty, Index_S, Index_r, epsilon_s, m_hat_infty, s_hat_sq_infty, v_hat_infty)
rslt_cpp
rslt_R
rslt_cpp - rslt_R


sum(ResultsGLR$Predictions * Alpha[[NrOfAlpha]])
sum(ResultsAdditive$Predictions * Alpha[[NrOfAlpha]])


