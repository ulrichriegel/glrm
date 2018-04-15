rm(list = ls())

start <- Sys.time()

library(stats)
library(readxl)
source('Simulation.R')
source('LR.R')
source('GLR.R')
source('MSE_GLR.R')

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

#Optionen:
c <- rep(1,n-1)
K <- 1:(n-1)
#K <- 1:2


ResultsGLR <- Generalized_LR_Method(S, v)
g_hat_infty <- ResultsGLR$g_hat_infty
h_hat_infty <- ResultsGLR$h_hat_infty
eps_hat_infty <- ResultsGLR$eps_hat_infty
s_hat_sq_infty <- ResultsGLR$s_hat_sq_infty
m_hat_infty <- ResultsGLR$m_hat_infty
v_hat_infty <- ResultsGLR$v_hat_infty




NumberOfSimulations <- 1000


#SimulatedData <- SimulateClaims(v, m, s_squared, epsilon, NumberOfSimulations)
SimulatedData <- SimulateClaims(v_hat_infty, m_hat_infty, s_hat_sq_infty, eps_hat_infty * c, NumberOfSimulations)
Predictions_GLR <- array(NA, dim = c(n,n,NumberOfSimulations))
Predictions_LR <- array(NA, dim = c(n,n,NumberOfSimulations))
eps_hat_GLR <- rep(NA, NumberOfSimulations)
m_hat_GLR <- matrix(NA, n, NumberOfSimulations)
v_hat_GLR <- matrix(NA, n, NumberOfSimulations)
s_hat_sq_GLR <- matrix(NA, n, NumberOfSimulations)
m_hat_LR <- matrix(NA, n, NumberOfSimulations)
s_hat_sq_LR <- matrix(NA, n, NumberOfSimulations)

for (SimNr in 1:NumberOfSimulations) {
  GLR <- Generalized_LR_Method_fixed_g_h(SimulatedData$Triangle[,,SimNr], v, K, g_hat_infty, h_hat_infty, m_hat_infty, s_hat_sq_infty)
  Predictions_GLR[,,SimNr] <- GLR$Predictions
  #eps_hat_GLR[SimNr] <- GLR$eps_hat_infty
  #m_hat_GLR[,SimNr] <- GLR$m_hat_infty
  #v_hat_GLR[,SimNr] <- GLR$v_hat_infty
  #s_hat_sq_GLR[,SimNr] <- GLR$s_hat_sq_infty
  LR <- LR_Method(SimulatedData$Triangle[,,SimNr],SimulatedData$volumes[,SimNr])
  Predictions_LR[,,SimNr] <- LR$Predictions
  m_hat_LR[,SimNr] <- LR$m_hat 
  s_hat_sq_LR[,SimNr] <- LR$s_hat_sq
  print(SimNr)
}

Ultimates_GLR <- apply(Predictions_GLR, c(1,3), 'sum')
Ultimates_LR <- apply(Predictions_LR, c(1,3), 'sum')
Ultimates_Sim <- apply(SimulatedData$Triangle, c(1,3), 'sum')

qa_GLR <- (Ultimates_Sim - Ultimates_GLR)^2
qa_LR <- (Ultimates_Sim - Ultimates_LR)^2
sqrt(apply(qa_GLR,1,'sum')/NumberOfSimulations)
sqrt(apply(qa_LR,1,'sum')/NumberOfSimulations)


#SimulatedData$volumes
#SimulatedData$Triangle
#SimulatedData$m
#SimulatedData$s_squared
#SimulatedData$epsilon

#versuch <- AdditiveMethod(SimulatedData$Triangle[,,1],SimulatedData$volumes[,1])
#versuch$Predictions

ende <- Sys.time()

ende - start




