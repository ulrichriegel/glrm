# #' This function provides the conventional loss ratio method
# #' @param S decumulated claims triangle (n x n Matrix)
# #' @param v vector containing the volume estimates
# #' @param W_used matrix containing the weights to be used (n x n Matrix)
# #' @export
LR_Method <- function(S, v, W_used = NULL) {
  ###########################################################################################################
  # This function provides the conventional loss ratio method
  ###########################################################################################################
  #
  # Required Parameters:
  # --------------------
  # S = incremental triange
  # v = vector of volumes
  #
  # Optional Parameters:
  # --------------------
  # W_used = if this parameter is not NULL then these weights are used
  # (otherwise W_opt is used - inversely proportional to the volumes)
  #
  # Output:
  # -------
  # Predictions = Full triangle including predictions
  # m_hat = vector containing the estimators for the m_k
  # s_hat_sq = vector containing the estimators for s_k^2
  # W_used = weights that have been used to calculate the predictions
  # W_opt = optimal weights in the additive model (inversely proportional to the volumes)

  n <- length(v)
  m <- rep(NA, n)
  s_squared <- rep(NA, n)

  W_opt <- matrix(0,ncol = n, nrow = n)
  for (j in 1:n) {
    W_opt[1:(n-j+1),j] <- v[1:(n-j+1)]/sum(v[1:(n-j+1)])
  }
  if (is.null(W_used)) {
    W_used <- W_opt
  }

  for (k in 1:n) {
    m[k] <- W_used[1:(n-k+1),k] %*% (S[1:(n-k+1),k]/v[1:(n-k+1)])
  }
  for (k in 1:(n-1)) {
    s_squared[k] <- sum((v[1:(n-k+1)]*(S[1:(n-k+1),k]/v[1:(n-k+1)]-m[k])^2))/(n-k)
  }
  # s_squared[n] <- min(s_squared[n-1]^2/s_squared[n-2],s_squared[n-1])
  s_squared[n] <- min(s_squared[n-1]^2/s_squared[n-2],s_squared[n-2])
  S_incl_predictions <- S
  for (k in 2:n) {
    S_incl_predictions[(n-k+2):n,k] <- v[(n-k+2):n] * m[k]
  }

  rslt <- list(Predictions = S_incl_predictions, m_hat = m, s_hat_sq = s_squared, W_used = W_used, W_opt = W_opt)
  return(rslt)
}


# #' This function calculates the mean squared error of the conventional loss ratio method in the extended additive model
# #' @param S decumulated claims triangle (n x n Matrix)
# #' @param W matrix containing the weights that have been used in the loss ratio method (n x n Matrix)
# #' @param Alpha matrix containing the weights alpha
# #' @param epsilon estimator for the parameter epsilon
# #' @param m vector of dimension n containing estimators for the parameters m_k
# #' @param s_squared vector of dimension n containing estimators for the s^2_k
# #' @param v vector containing the volume estimates
# #' @export

MSE_LR <- function (S, W, Alpha, epsilon, m, s_squared, v) {
  ###########################################################################################################
  # This function calculates the mean squared error of the loss ratio method in the extended additive model
  ###########################################################################################################
  #
  # Required Parameters:
  # --------------------
  # S = incremental Triangle
  # W: n x n matrix containing the weights W[i,j] for estimators S_ij/v_i in the estimation of m[k]
  # Alpha: n x n matrix containing Alpha[i,j] for i+j > n+1
  # m = vector with estimators for m_k
  # epsilon = vector with epsilons
  # Note: the MSE of the LR method in the conventional additive model is obtained with epsilon = rep(0,n-1)
  # v = vector with volume estimates
  #
  # Optional Parameters:
  # --------------------
  #
  # Output:
  # -------
  # MSE = mean squared error for the given Alpha (in the extended additive model)
  # RandomError = parameter variance
  # EstimationError = squared parameter estimation error
  # CovMatrix_S_hat = Covariance matrix of the predictors in the extended additive model
  # L = matrix containing the used order of the set L = {(i,k)|i+k>n+1}:
  #   - 1st column: number of the element
  #   - 2nd & 3rd column: element of L

  n <- length(v)

  #ensure that weights sum up to 1 for each column:
  for (k in 1:n) {
    W[1:(n-k+1),k] <- W[1:(n-k+1),k]/sum(W[1:(n-k+1),k])
  }

  # calculate prod1plusepssq[i,j] = (1+epsilon[i]^2) * ... * (1+epsilon[j]^2)
  epsilon_squared <- epsilon^2
  prod1plusepssq <- matrix(1,ncol=n,nrow=n)
  for (i in 1:(n-1)) {
    for (j in (i:(n-1))) {
      prod1plusepssq[i,j] <- prod(1+epsilon_squared[i:j])

    }
  }

  # ensure that Alpha[i,k] = 0 for i+k <= n+1
  Alpha <- (row(Alpha) + col(Alpha) > n+1) * Alpha

  # index matrix containing the set L in columns 2&3:
  L <- matrix(c(1:(n*(n-1)/2),row(S)[col(S)+row(S)>n+1],col(S)[col(S)+row(S)>n+1]),ncol=3)
  NumberRowsL <- nrow(L)

  # Random error:
  RandomError <- sum(Alpha^2 * (v %*% t(s_squared)))

  # Estimation error
  EstimationError <- 0
  CovMatrix_S_hat <- matrix(NA, nrow = NumberRowsL, ncol = NumberRowsL)
  for (p in 1:NumberRowsL) {
    i <- L[p,2]
    k <- L[p,3]
    for (q in 1:NumberRowsL) {
      j <- L[q,2]
      l <- L[q,3]
      Temp <- 0
      for (nu in 1:(n-k+1)) {
        for (kappa in 1:(n-l+1)) {
          Temp <- Temp + W[nu,k] * m[k] * W[kappa,l] * m[l] * (prod1plusepssq[max(nu,kappa),min(i,j)-1]-1)
        }
      }
      if (k==l) {
        for (nu in 1:(n-k+1)) {
          Temp <- Temp + W[nu,k]^2 * s_squared[k] / v[nu] * prod1plusepssq[nu,min(i,j)-1]
        }
      }
      EstimationError <- EstimationError + Temp * Alpha[i,k] * Alpha[j,l] * v[i] * v[j]
      CovMatrix_S_hat[p,q] <- Temp * v[i] * v[j]
    }
  }

  MSE <- RandomError + EstimationError
  return(list(MSE = MSE, RandomError = RandomError, EstimationError = EstimationError, CovMatrix_S_hat = CovMatrix_S_hat, L = L))
}




OptimalWeights_LR_Method <- function(v, m, s_squared, epsilon) {
  n <- length(m)
  OptWeights_Matrix <- matrix(0,ncol = n, nrow = n)

  epsilon_sqared <- epsilon^2

  prod1plusepssq <- rep(1,n)
  for (nu in 1:(n-1)) {
    prod1plusepssq[nu] <- prod(1+epsilon_sqared[nu:(n-1)])
  }

  for (k in 1:(n-1)) {
    # Kovarianzmatrix \Sigma_k
    CovMatrix <- matrix(0, ncol = n-k+1, nrow = n-k+1)
    for (i in 1:(n-k+1-1)) {
      for (j in (i+1):(n-k+1)) {
        CovMatrix[i,j] <- m[k]^2 * (prod1plusepssq[max(i,j)]-1)
      }
    }
    CovMatrix <- CovMatrix + t(CovMatrix)
    for (i in 1:(n-k+1)) {
      CovMatrix[i,i] <- (m[k]^2 + s_squared[k]/v[i]) * prod1plusepssq[i] - m[k]^2
    }
    OptWeights_Matrix[1:(n-k+1),k] <- CalculateOptimalWeights(CovMatrix)
  }
  OptWeights_Matrix[1,n] <- 1
  return(OptWeights_Matrix)
}



