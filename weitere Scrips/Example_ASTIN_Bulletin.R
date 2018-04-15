rm(list = ls())

library(stats)
library(readxl)
library(ChainLadder)

source('LR.R')
source('GLR.R')
source('MSE_GLR.R')

Create_Table_v_S <- T                # Table 1
Create_Table_Iteration_Eps <- F      # Table 2
Create_Table_Parameters <- T         # Table 3 
Create_Table_W_opt <- T              # Table 4
Create_Table_W_conventional <- T     # Table 5
Create_Table_Predictions <- T        # Table 6
Create_Table_Predictions_GLR <- T    # Table 7
Create_Table_MSE_GLR <- T            # Table 8

# Options:
NrOfIterations <- 5                # for ASTIN Bulletin: 200

# Bsp aus ASTIN Bulletin einlesen (v = Volumina, C = kumulierte Zahlungen)
v <- read_excel("./Examples/Example1/v.xlsx", sheet = 1)
C <- read_excel("./Examples/Example1/C.xlsx", sheet = 1)

v <- as.matrix(v)
C <- as.matrix(C)
v <- v/1000000
C <- C/1000000
n <- ncol(C)
S <- matrix(0, ncol = n, nrow = n)
for (k in n:2) {
  S[1:(n-k+1),k] <- C[1:(n-k+1),k]-C[1:(n-k+1),k-1]
}
S[,1] <- C[,1]

c <- rep(1,n-1)


if (Create_Table_v_S) {
  Table_v_S <- round(cbind(v *1000, S * 1000), digits = 0)
  View(Table_v_S)
}

# apply additive method
ResultsAdditive <- LR_Method(S,v)
m_hat <- ResultsAdditive$m_hat
s_hat_sq <- ResultsAdditive$s_hat_sq
UltimatesAdditiveMethod <- apply(ResultsAdditive$Predictions, 1, 'sum')
W_conventional <- ResultsAdditive$W_opt
if (Create_Table_W_conventional) {
  Table_W_conventional <- W_conventional
  View(Table_W_conventional)  
}

if(Create_Table_Iteration_Eps) {
  Table_Iteration_Eps <- NULL
  for (eps_start in c(0.005, 0.02, 0.05, 0.1)) {
    ResultsGLR <- Generalized_LR_Method(S, v, Iterations = NrOfIterations, eps_start = eps_start)  
    Table_Iteration_Eps <- c(Table_Iteration_Eps, ResultsGLR$eps_hat[c(1,2,3,4,6,11,21,51,101)])
  }
  Table_Iteration_Eps <- matrix(Table_Iteration_Eps, nrow = 4, byrow = T)
  View(Table_Iteration_Eps)
}

ResultsGLR <- Generalized_LR_Method(S, v, Iterations = NrOfIterations)
K <- ResultsGLR$K
J <- ResultsGLR$J
J_all <- ResultsGLR$J_all

g_hat_infty <- ResultsGLR$g_hat_infty
h_hat_infty <- ResultsGLR$h_hat_infty
eps_hat_infty <- ResultsGLR$eps_hat_infty
s_hat_sq_infty <- ResultsGLR$s_hat_sq_infty
m_hat_infty <- ResultsGLR$m_hat_infty
v_hat_infty <- ResultsGLR$v_hat_infty

W_opt <- OptimalWeights_LR_Method(v, m_hat_infty, s_hat_sq_infty, eps_hat_infty * c)
if (Create_Table_W_opt) {
  Table_W_opt <- W_opt 
  View(Table_W_opt)
}

ResultsAdditiveOpt <- LR_Method(S, v, W_opt)

if (Create_Table_Parameters) {
  Table_Parameters <- c( round(v * 1000, digits = 0), s_hat_sq)
  eps_fix <- c(0.01, 0.0306, 0.05)
  for (nu in 1:length(eps_fix)) {
    Temp <- Generalized_LR_Method(S, v, Iterations = NrOfIterations, eps_ext = eps_fix[nu], eps_start = eps_fix[nu])
    Table_Parameters <- c(Table_Parameters, round(Temp$v_hat_infty * 1000, digits = 0), Temp$s_hat_sq_infty)
  }
  Table_Parameters <- matrix(Table_Parameters, nrow = n)
  View(Table_Parameters)
  rm(Temp)
}

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

if (Create_Table_Predictions) {
  Table_Predictions <- NULL
  for (NrOfAlpha in 1:(n+1)) {
    # Predictions conventional additive with canonical weights
    temp <- sum(ResultsAdditive$Predictions * Alpha[[NrOfAlpha]])
    Table_Predictions <- c(Table_Predictions, temp)
    # Predictions conventional additive with optimal weights
    temp <- sum(ResultsAdditiveOpt$Predictions * Alpha[[NrOfAlpha]])
    Table_Predictions <- c(Table_Predictions, temp)
    # conventional formula for standard error of loss ratio method
    temp <- sqrt( MSE_LR(S, W_conventional, Alpha[[NrOfAlpha]], rep(0.0,n-1), m_hat, s_hat_sq, v)$MSE)
    Table_Predictions <- c(Table_Predictions, temp)
    # conventional formula for standard error of loss ratio method
    temp <- sqrt( MSE_LR(S, W_opt, Alpha[[NrOfAlpha]], rep(0.0,n-1), m_hat, s_hat_sq, v)$MSE)
    Table_Predictions <- c(Table_Predictions, temp)
    # formula for standard error of loss ratio method using the estimators from GLR
    temp <- sqrt( MSE_LR(S, W_conventional, Alpha[[NrOfAlpha]], eps_hat_infty * c, m_hat_infty, s_hat_sq_infty, v_hat_infty)$MSE)
    Table_Predictions <- c(Table_Predictions, temp)
    # formula for standard error of loss ratio method using the estimators from GLR
    temp <- sqrt( MSE_LR(S, W_opt, Alpha[[NrOfAlpha]], eps_hat_infty * c, m_hat_infty, s_hat_sq_infty, v_hat_infty)$MSE)
    Table_Predictions <- c(Table_Predictions, temp)
  }
  Table_Predictions <- matrix(Table_Predictions, nrow = n+1, byrow = T)
  View(Table_Predictions)
}

Cum <- incr2cum(as.triangle(S))
Cum[row(Cum)+col(Cum)>n+1] <- NA 
ResultsChainLadder <- MackChainLadder(Cum, est.sigma="Mack")
ResultsChainLadder
CLFullIncrementalTriangle <- cum2incr(ResultsChainLadder$FullTriangle)

if (Create_Table_Predictions_GLR) {
  Table_Predictions_GLR <- numeric()
  temp1 <- Generalized_LR_Method(S,v,eps_start = 0.01, eps_ext = 0.01, Iterations = NrOfIterations)
  temp2 <- Generalized_LR_Method(S,v, Iterations = NrOfIterations)
  temp3 <- Generalized_LR_Method(S,v,eps_start = 0.05, eps_ext = 0.05, Iterations = NrOfIterations)
  
  for (NrOfAlpha in 1:(n+1)) {
    Table_Predictions_GLR <- c(Table_Predictions_GLR, sum(ResultsAdditive$Predictions * Alpha[[NrOfAlpha]]))
    Table_Predictions_GLR <- c(Table_Predictions_GLR, sum(temp1$Predictions * Alpha[[NrOfAlpha]]))
    Table_Predictions_GLR <- c(Table_Predictions_GLR, sum(temp2$Predictions * Alpha[[NrOfAlpha]]))
    Table_Predictions_GLR <- c(Table_Predictions_GLR, sum(temp3$Predictions * Alpha[[NrOfAlpha]]))
    Table_Predictions_GLR <- c(Table_Predictions_GLR, sum(CLFullIncrementalTriangle * Alpha[[NrOfAlpha]]))
  }
  Table_Predictions_GLR <- matrix(Table_Predictions_GLR, nrow = n+1, byrow = T)
  View(Table_Predictions_GLR)
}




if (Create_Table_MSE_GLR) {
  Table_MSE_GLR <- numeric()
  for (NrOfAlpha in 1:(n+1)) {
    MSE_GLR_temp <- MSE_GLR(S, Alpha[[NrOfAlpha]], J,  K,  g_hat_infty, h_hat_infty, eps_hat_infty * c, m_hat_infty, s_hat_sq_infty, v_hat_infty, CovDotMatrixName = "./Examples/Example1/CovdotMatrix 3,06% 200 Iterations.txt")
    Table_MSE_GLR <- c(Table_MSE_GLR,sqrt(MSE_GLR_temp$MSE))
    Table_MSE_GLR <- c(Table_MSE_GLR,sqrt(MSE_GLR_temp$RandomError))
    Table_MSE_GLR <- c(Table_MSE_GLR,sqrt(MSE_GLR_temp$EstimationError))
    Table_MSE_GLR <- c(Table_MSE_GLR,MSE_GLR_temp$Bias)
  }
  Table_MSE_GLR <- matrix(Table_MSE_GLR, ncol = 4, byrow = T)
  Table_MSE_GLR <- cbind(Table_MSE_GLR, c(ResultsChainLadder$Mack.S.E[,n],ResultsChainLadder$Total.Mack.S.E))
  View(Table_MSE_GLR)
}

