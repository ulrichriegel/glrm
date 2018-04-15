MSE_GLR <- function (S, Alpha, K,  g, h, epsilon, m, s_squared, v, CovdotMatrix = NULL, UseAllEstimators = FALSE, blnUseRcpp = T) {
  ###########################################################################################################
  # This function calculates MSE^\bullet for the generalized loss ratio method
  # This function calls the other functions of this file.
  ###########################################################################################################
  #
  # Required Parameters:
  # --------------------
  # S = incremental triange
  # Alpha: n x n matrix containing Alpha[i,j] for i+j > n+1
  # K = subset of {1,...,n-1} such that S_{i,k} is lognormally distributed for all (i,k) with k\in K
  # g = list containing the used weight vectors (typically \hatg_i^{(\infty)})
  # h = list containing the used weight vectors (typically \hath_i^{(\infty)})
  # epsilon
  # m = vector containing the used estimators for m_k (typically \hatm_k^{(\infty)})
  # s_squared
  # v = vector containing the used estimators for v_i (typically \hatv_k^{(\infty)})
  #
  # Optional Parameters:
  # --------------------
  # CovDotMatrixName: covariance matrix of the predictions S_{i,k}^{(\infty)}: Cov^\bullet(S_{i,k}^{(\infty)},S_{j,l}^{(\infty)});
  # if this is NULL then the covariance matrix is calculated using the function Create_Cov_Matrix_S_infty;
  # this is computationally very intensive (might take several hours)
  #
  # Output:
  # -------
  # MSE = mean squared error for the given Alpha
  # ProcessVariance = Process variance
  # PredictionVariance = Variance of the prediction
  # Bias = bias of the prediction (not squared)
  # CovdotMatrix used covariance matrix of the predictions
  # L = matrix containing the used order of the set L = {(i,k)|i+k>n+1}:
  #   - 1st column: number of the element
  #   - 2nd & 3rd column: element of L

  #print(Sys.time())

  n <- length(v)

  # Check and reduce K if necessary (if not all incremental claims are positive)
  K <- CheckK(S, K)
  # Matrix J[[i]] contains the set $J_i$
  J <- CreateJ(S, K, UseAllEstimators)
  # J[[i]] = used list of indices for the estimators of v_i^{-1}; for \hatr_{i,(j,k)}
  #         - J[[i]][,1] = i, i.e. estimator for v_i^{-1}
  #         - J[[i]][,2] = j
  #         - J[[i]][,3] = k


  # index matrix containing the set L in columns 2&3:
  L <- matrix(c(1:(n*(n-1)/2),row(S)[col(S)+row(S)>n+1],col(S)[col(S)+row(S)>n+1]),ncol=3)
  NumberRowsL <- nrow(L)

  # ensure that Alpha[i,k] = 0 for i+k <= n+1
  Alpha <- (row(Alpha) + col(Alpha) > n+1) * Alpha

  # Process Variance:
  ProcessVariance <- sum(Alpha^2 * (v %*% t(s_squared)))

  Alpha_Vec <- rep(NA,NumberRowsL)
  for (i in 1:NumberRowsL) {
    Alpha_Vec[i] <- Alpha[L[i,2],L[i,3]]
  }

  if (is.null(CovdotMatrix)) {
    CovdotMatrix <- Create_Cov_Matrix_S_infty(J,  K,  g, h, epsilon, m, s_squared, v, blnUseRcpp)$CovMatrix
  }

  PredictionVariance <- t(Alpha_Vec) %*% CovdotMatrix %*% Alpha_Vec

  Bias_Vec <- Create_Edot_Bias_Vector_S_infty(J,  K,  g, h, epsilon, m, s_squared, v)$Bias
  Bias <- Bias_Vec %*% Alpha_Vec
  Bias_Squared <- Bias^2

  MSE <- ProcessVariance + PredictionVariance + Bias_Squared

  return(list(MSE = MSE, ProcessVariance = ProcessVariance, PredictionVariance = PredictionVariance, Bias = Bias, CovdotMatrix = CovdotMatrix, L = L))
}










Create_Edot_Bias_Vector_S_infty <- function(J,  K,  g, h, epsilon, m, s_squared, v) {

  # index matrix containing the set L in columns 2&3:
  n <- length(v)

  M <- matrix(0, nrow = n, ncol = n)
  L <- matrix(c(1:(n*(n-1)/2),row(M)[col(M)+row(M)>n+1],col(M)[col(M)+row(M)>n+1]),ncol=3)
  NumberRowsL <- nrow(L)

  Edot <- rep(NA, NumberRowsL)
  Bias <- rep(NA, NumberRowsL)


  for (mu in 1:NumberRowsL) {
    Index_S_hat_infty <- L[mu,2:3]
    Edot[mu] <- Calculate_Edot_S_infty(J,  K,  g, h, Index_S_hat_infty, epsilon, m, s_squared, v)
    Bias[mu] <- Edot[mu] - v[Index_S_hat_infty[1]] * m[Index_S_hat_infty[2]]
  }
  return(list(Edot = Edot, Bias = Bias, L = L))
}



Calculate_Edot_S_infty <- function(J,  K,  g, h, Index_S_hat_infty, epsilon, m, s_squared, v) {

  n <- length(v)

  # calculate gammas
  # set gam[i,l] <- 1 if $l\not\in K$
  gam <- matrix(1, ncol = n, nrow = n)
  for (i in 1:n) {
    for (k in K) {
      gam[i,k] <- 1+ s_squared[k]/(v[i]*m[k]^2)
    }
  }

  # calculate squared epsilons
  epsilon_sqared <- epsilon^2

  # calculate \prod_{\nu=i}^{n-1} (1+\varepsilon_\nu^2)
  # set to 1 if i=n
  prod1plusepssq <- rep(1,n)
  for (i in 1:(n-1)) {
    prod1plusepssq[i] <- prod(1+epsilon_sqared[i:(n-1)])
  }

  i<- Index_S_hat_infty[1]
  k<- Index_S_hat_infty[2]

  Edot <- 0

  for (mu in 1:(n-k+1)) {
    Temp <- 0
    Index_S <- c(mu,k)
    Index_prodv <- mu
    Temp <- Temp + 2 * E_S_by_vhat(J, g, Index_S, Index_prodv, prod1plusepssq, gam, K, m, s_squared, v)

    Index_prodv <- c(mu,i)
    Temp <- Temp - v[i] * E_S_by_vhat(J, g, Index_S, Index_prodv, prod1plusepssq, gam, K, m, s_squared, v)


    Edot <- Edot + h[[k]][mu] * Temp
  }
  Edot <- v[i] * Edot
  return(Edot)

}



E_S_by_vhat <- function(J, g, Index_S, Index_prodv, prod1plusepssq, gam, K, m, s_squared, v) {
  # Calculates E(S_{i,k}/prod_hatv)
  # prod_hatv is either of the form hatv_i or hatv_i*hatv_j
  # Index_S:  contains (i,k)
  # Index_prodv:  if prod_hatv = hatv_i then Index_v = i, if prod_hatv = hatv_i*hatv_j then Index_v = (i,j),
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # gam[i,k] must be 1 for k \not\in K
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Index_S <- matrix(Index_S, nrow = 1)

  E_S_by_vhat <- 0

  if (length(Index_prodv) == 1) {
    for (j in 1:nrow(J[[Index_prodv[1]]])) {
      Index_r <- matrix(J[[Index_prodv[1]]][j,],ncol = 3)
      Temp <- Calculate_E_Sr(Index_S, Index_r, prod1plusepssq, gam, K, m, s_squared, v)
      E_S_by_vhat <- E_S_by_vhat + g[[Index_prodv[1]]][j] * Temp
    }
  } else if (length(Index_prodv) == 2) {
    for (j1 in 1:nrow(J[[Index_prodv[1]]])) {
      for (j2 in 1:nrow(J[[Index_prodv[2]]])) {
        Index_r <- matrix(c(J[[Index_prodv[1]]][j1,],J[[Index_prodv[2]]][j2,]),ncol = 3,byrow = T)
        Temp <- Calculate_E_Sr(Index_S, Index_r, prod1plusepssq, gam, K, m, s_squared, v)
        E_S_by_vhat <- E_S_by_vhat + g[[Index_prodv[1]]][j1] * g[[Index_prodv[2]]][j2] * Temp
      }
    }
  } else {
    stop("prod_hatv has more than 2 factors")
  }
  return(E_S_by_vhat)
}



Create_Cov_Matrix_S_infty <- function(J,  K,  g, h, epsilon, m, s_squared, v, blnUseRcpp = T) {

  n <- length(v)

  # index matrix containing the set L in columns 2&3:
  M <- matrix(0, nrow = n, ncol = n)
  L <- matrix(c(1:(n*(n-1)/2),row(M)[col(M)+row(M)>n+1],col(M)[col(M)+row(M)>n+1]),ncol=3)
  NumberRowsL <- nrow(L)

  CovMatrix <- matrix(0, nrow = NumberRowsL, ncol = NumberRowsL)

  starttime <- Sys.time()
  #print(c("Number of Rows:", NumberRowsL))
  print("Calculating covariance matrix of GLRM predictions")

  #blnUseRcpp <- T # use C++ function for better performance




  if (blnUseRcpp) {
    J_all <- J[[1]]
    for (i in 2:n) {
      J_all <- rbind(J_all,J[[i]])
    }
    J_s <- c(J_all)
    K_s <- K
    Number_of_Est_s <- numeric(n)
    g_s <- numeric()
    h_s <- numeric()

    for (i in 1:n) {
      Number_of_Est_s[i] <- nrow(J[[i]])
      g_s <- c(g_s,g[[i]])
      h_s <- c(h_s,h[[i]])
    }
    epsilon_s <- epsilon
    m_s <- m
    s_squared_s <- s_squared
    v_s <- v

    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    progress <- 0
    step_pb <- 1 / (NumberRowsL * (NumberRowsL + 1) / 2)
    setTxtProgressBar(pb, progress)
    for (mu in 1:NumberRowsL) {
      for (nu in mu:NumberRowsL) {
        Index_S_hat_infty1 <- c(L[mu,2:3])
        Index_S_hat_infty2 <- c(L[nu,2:3])
        CovMatrix[mu,nu] <- Cov_S_hat_infty1_S_hat_infty2_cpp(J_s, K_s, g_s, h_s, Index_S_hat_infty1, Index_S_hat_infty2, epsilon_s, m_s, s_squared_s, v_s, n, Number_of_Est_s)
        #print(c("Nr.", mu, nu, "S", Index_S_hat_infty1, "S", Index_S_hat_infty2, "Covariance", CovMatrix[mu,nu]))
        #write.table(CovMatrix, file = "CovDotMatrix.csv", sep = ";", row.names = FALSE)
        progress <- progress + step_pb
        setTxtProgressBar(pb, progress)
      }
    }
    close(pb)
  } else {
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    progress <- 0
    step_pb <- 1 / (NumberRowsL * (NumberRowsL + 1) / 2)
    setTxtProgressBar(pb, progress)

    for (mu in 1:NumberRowsL) {
      for (nu in mu:NumberRowsL) {

        Index_S_hat_infty1 <- L[mu,2:3]
        Index_S_hat_infty2 <- L[nu,2:3]
        CovMatrix[mu,nu] <- Cov_S_hat_infty1_S_hat_infty2(J,  K,  g, h, Index_S_hat_infty1, Index_S_hat_infty2, epsilon, m, s_squared, v)
        # finished <- Sys.time()
        # print(c("Nr.", mu, nu, "S", Index_S_hat_infty1, "S", Index_S_hat_infty2, "Covariance:", CovMatrix[mu,nu]))
        # start <- finished
        write.table(CovMatrix, file = "CovdotMatrix.csv", sep = ";", row.names = FALSE)

        progress <- progress + step_pb
        setTxtProgressBar(pb, progress)
      }
    }

    close(pb)

  }
  D <- diag(diag(CovMatrix))
  CovMatrix <- CovMatrix + t(CovMatrix) - D
  return(list(CovMatrix = CovMatrix, L = L))
}



Cov_S_hat_infty1_S_hat_infty2 <- function(J,  K,  g, h, Index_S_hat_infty1, Index_S_hat_infty2, epsilon, m, s_squared, v) {
  n <- length(v)

  # calculate gammas
  # set gam[i,l] <- 1 if $l\not\in K$
  gam <- matrix(1, ncol = n, nrow = n)
  for (i in 1:n) {
    for (k in K) {
      gam[i,k] <- 1+ s_squared[k]/(v[i]*m[k]^2)
    }
  }

  # calculate squared epsilons
  epsilon_sqared <- epsilon^2

  # calculate \prod_{\nu=i}^{n-1} (1+\varepsilon_\nu^2)
  # set to 1 if i=n
  prod1plusepssq <- rep(1,n)
  for (i in 1:(n-1)) {
    prod1plusepssq[i] <- prod(1+epsilon_sqared[i:(n-1)])
  }

  i<- Index_S_hat_infty1[1]
  k<- Index_S_hat_infty1[2]
  j<- Index_S_hat_infty2[1]
  l<- Index_S_hat_infty2[2]

  Covariance <- 0

  for (mu in 1:(n-k+1)) {
    for (nu in 1:(n-l+1)) {
      Temp <- 0
      Index_S1 <- c(mu,k)
      Index_prodv1 <- mu
      Index_S2 <- c(nu,l)
      Index_prodv2 <- nu
      Temp <- Temp + 4 * Cov_S_by_vhat_S_by_vhat(J, g, Index_S1, Index_prodv1, Index_S2, Index_prodv2, prod1plusepssq, gam, K, m, s_squared, v)
      #print(Temp)

      Index_prodv1 <- c(mu,i)
      Index_prodv2 <- nu
      Temp <- Temp - 2 * v[i] * Cov_S_by_vhat_S_by_vhat(J, g, Index_S1, Index_prodv1, Index_S2, Index_prodv2, prod1plusepssq, gam, K, m, s_squared, v)
      #print(Temp)

      Index_prodv1 <- mu
      Index_prodv2 <- c(nu,j)
      Temp <- Temp - 2 * v[j] * Cov_S_by_vhat_S_by_vhat(J, g, Index_S1, Index_prodv1, Index_S2, Index_prodv2, prod1plusepssq, gam, K, m, s_squared, v)
      #print(Temp)

      Index_prodv1 <- c(mu,i)
      Index_prodv2 <- c(nu,j)
      Temp <- Temp + v[i] * v[j] * Cov_S_by_vhat_S_by_vhat(J, g, Index_S1, Index_prodv1, Index_S2, Index_prodv2, prod1plusepssq, gam, K, m, s_squared, v)
      #print(Temp)

      Covariance <- Covariance + h[[k]][mu] * h[[l]][nu] * Temp
      #print(Covariance)
    }
  }
  Covariance <- v[i] * v[j] * Covariance
  return(Covariance)
}




Cov_S_by_vhat_S_by_vhat <- function(J, g, Index_S1, Index_prodv1, Index_S2, Index_prodv2, prod1plusepssq, gam, K, m, s_squared, v) {
  # Calculates the Cov(S_{i1,k1}/prod_hatv1,S_{i2,k2}/prod_hatv2)
  # prod_hatv1 is either of the form hatv_i or hatv_i*hatv_j
  # prod_hatv2 is either of the form hatv_i or hatv_i*hatv_j
  # Index_S1:  contains (i1,k1)
  # Index_S2:  contains (i2,k2)
  # Index_prodv1:  if prod_hatv1 = hatv_i then Index_v1 = i, if prod_hatv1 = hatv_i*hatv_j then Index_v1 = (i,j),
  # Index_prodv2:  if prod_hatv2 = hatv_i then Index_v2 = i, if prod_hatv2 = hatv_i*hatv_j then Index_v2 = (i,j),
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # gam[i,k] must be 1 for k \not\in K
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Covariance <- 0
  Counter <- 0

  Index_S1 <- matrix(Index_S1, nrow = 1)
  Index_S2 <- matrix(Index_S2, nrow = 1)
  Index_S12 <- rbind(Index_S1, Index_S2)


  if (length(Index_prodv1) == 1 && length(Index_prodv2) == 1) {
    for (j1 in 1:nrow(J[[Index_prodv1[1]]])) {
      Index_r1 <- matrix(J[[Index_prodv1[1]]][j1,],ncol = 3)
      for (k1 in 1:nrow(J[[Index_prodv2[1]]])) {
        Index_r2 <- matrix(J[[Index_prodv2[1]]][k1,],ncol = 3)
        Index_r12 <- rbind(Index_r1, Index_r2)
        #E_1 <- Calculate_E_Sr_bc(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v)
        #E_2 <- Calculate_E_Sr_bc(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v)
        #E_12 <- Calculate_E_Sr_bc(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v)
        E_1 <- Calculate_E_Sr(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v)
        E_2 <- Calculate_E_Sr(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v)
        E_12 <- Calculate_E_Sr(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v)
        Counter <- Counter + 1
        if (Counter < 0) {
        #print ("E1")
      #print("E2")
        #print(E_2)
          #print("E12")
          #print(E_12)
          #print("g1")
          #print(g[[Index_prodv1[1]]][j1])
          #print("g2")
          #print(g[[Index_prodv2[1]]][k1])
        }
        Covariance <- Covariance + g[[Index_prodv1[1]]][j1] * g[[Index_prodv2[1]]][k1] * (E_12 - E_1 * E_2)
        #print(c(Counter, Covariance))
      }
    }
  } else if (length(Index_prodv1) == 2 && length(Index_prodv2) == 2) {
    for (j1 in 1:nrow(J[[Index_prodv1[1]]])) {
      for (j2 in 1:nrow(J[[Index_prodv1[2]]])) {
        Index_r1 <- matrix(c(J[[Index_prodv1[1]]][j1,],J[[Index_prodv1[2]]][j2,]),ncol = 3,byrow = T)
        for (k1 in 1:nrow(J[[Index_prodv2[1]]])) {
          for (k2 in 1:nrow(J[[Index_prodv2[2]]])) {
            Index_r2 <- matrix(c(J[[Index_prodv2[1]]][k1,],J[[Index_prodv2[2]]][k2,]), ncol = 3, byrow = T)
            Index_r12 <- rbind(Index_r1, Index_r2)
            #E_1 <- Calculate_E_Sr_bc(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v)
            #E_2 <- Calculate_E_Sr_bc(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v)
            #E_12 <- Calculate_E_Sr_bc(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v)
            E_1 <- Calculate_E_Sr(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v)
            E_2 <- Calculate_E_Sr(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v)
            E_12 <- Calculate_E_Sr(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v)
            Covariance <- Covariance + g[[Index_prodv1[1]]][j1] * g[[Index_prodv1[2]]][j2] * g[[Index_prodv2[1]]][k1] * g[[Index_prodv2[2]]][k2] * (E_12 - E_1 * E_2)
          }
        }
      }
    }

  } else {
    if (length(Index_prodv1) == 1) {
      Temp <- Index_S1
      Index_S1 <- Index_S2
      Index_S2 <- Temp
      Temp <- Index_prodv1
      Index_prodv1 <- Index_prodv2
      Index_prodv2 <- Temp
    }
    for (j1 in 1:nrow(J[[Index_prodv1[1]]])) {
      for (j2 in 1:nrow(J[[Index_prodv1[2]]])) {
        Index_r1 <- matrix(c(J[[Index_prodv1[1]]][j1,],J[[Index_prodv1[2]]][j2,]),ncol = 3,byrow = T)
        for (k1 in 1:nrow(J[[Index_prodv2[1]]])) {
          Index_r2 <- matrix(J[[Index_prodv2[1]]][k1,], ncol = 3, byrow = T)
          Index_r12 <- rbind(Index_r1, Index_r2)
          #E_1 <- Calculate_E_Sr_bc(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v)
          #E_2 <- Calculate_E_Sr_bc(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v)
          #E_12 <- Calculate_E_Sr_bc(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v)
          E_1 <- Calculate_E_Sr(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v)
          E_2 <- Calculate_E_Sr(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v)
          E_12 <- Calculate_E_Sr(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v)
          Covariance <- Covariance + g[[Index_prodv1[1]]][j1] * g[[Index_prodv1[2]]][j2] * g[[Index_prodv2[1]]][k1] * (E_12 - E_1 * E_2)
        }
      }
    }
  }
  return(Covariance)
}






Calculate_E_Sr <- function(Index_S, Index_r, prod1plusepssq, gam, K, m, s_squared, v) {
  # Index_S:  if S appears l times then lx2 matrix; columns contain the two indices of S
  #           the function works with l=1 and l=2
  # Index_r:  if index is in J^k then matrix with k rows and 3 columns;
  #           for the estimator $\hatr_{i,(j,k)}$ the columns contain i,j and k
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # gam[i,k] must be 1 for k \not\in K
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  n <- length(v)
  chi <- matrix(0,ncol = n,nrow = n)
  xi <- matrix(0,ncol = n,nrow = n)

  Count_r <- nrow(Index_r)
  Count_S <- nrow(Index_S)

  if (Count_S <1 || Count_S >2) {
    stop("Number of S not admissible")
  }

  for (nu in 1:Count_r) {
    i <- Index_r[nu,1]
    j <- Index_r[nu,2]
    k <- Index_r[nu,3]
    if (j!=0) {
      chi[j,k] <- chi[j,k] + 1
      chi[i,k] <- chi[i,k] - 1
      xi[i,k] <- xi[i,k] - 1
    }
  }

  for (nu in 1: Count_S) {
    i <- Index_S[nu,1]
    k <- Index_S[nu,2]
    if (k %in% K) {
      chi[i,k] <- chi[i,k] + 1
    }
  }


  rho <- 0.5* chi * (chi - 1)

  #print(chi)


  Temp <- prod(gam^(rho+xi))


  pi <- ifelse(Index_r[,2]==0,Index_r[,1],Index_r[,2])

  # Calculate beta
  beta <- 1
  if (Count_r > 1) {
    i_vec <- sort(pi)
    for (kappa in 2:Count_r) {
      beta <- beta * prod1plusepssq[i_vec[kappa]]^(kappa-1)
    }
  }

  #print(beta)
  #print(gam)

  tau <- Index_r[,1]

  E_Sr <- beta / prod(v[tau]) * prod(v[Index_S[,1]]) * Temp

  if (Count_S > 1) {
    if ((Index_S[1,1] == Index_S[2,1]) && (Index_S[1,2] == Index_S[2,2]) && (Index_S[1,2] %in% (1:n)[-K])) {
      i <- Index_S[1,1]
      l <- Index_S[1,2]
      E_Sr <- E_Sr * (m[l]^2 + s_squared[l] / v[i])
    } else {
      E_Sr <- E_Sr * prod(m[Index_S[,2]])
      #print(E_Sr)
    }
  } else {
    E_Sr <- E_Sr * prod(m[Index_S[,2]])
  }
  return(E_Sr)
}























Test_R <- function(J,  K,  g, h, Index_S, Index_r, epsilon, m, s_squared, v) {
  n <- length(v)

  # calculate gammas
  # set gam[i,l] <- 1 if $l\not\in K$
  gam <- matrix(1, ncol = n, nrow = n)
  for (i in 1:n) {
    for (k in K) {
      gam[i,k] <- 1+ s_squared[k]/(v[i]*m[k]^2)
    }
  }

  # calculate squared epsilons
  epsilon_sqared <- epsilon^2

  # calculate \prod_{\nu=i}^{n-1} (1+\varepsilon_\nu^2)
  # set to 1 if i=n
  prod1plusepssq <- rep(1,n)
  for (i in 1:(n-1)) {
    prod1plusepssq[i] <- prod(1+epsilon_sqared[i:(n-1)])
  }



  return(Calculate_E_Sr(Index_S, Index_r, prod1plusepssq, gam, K, m, s_squared, v))


}


