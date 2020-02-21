# #' This function provides the generalized loss ratio method
# #' @param S decumulated claims triangle (n x n Matrix)
# #' @param v vector containing the volume estimates
# #' @param K set containing the columns with lognormally distributed increments
# #' @export
Generalized_LR_Method <- function(S, v, K = NULL, c = NULL, eps_start = 0.05, lower_bound_eps = 0.001, lower_bound_cv_S = 0.0001, eps_ext = NULL, adjust_v_hat_infty_for_small_eps_ext = F, s_sq_ext = NULL, upper_bound_eps = 0.1, upper_bound_cv_S = 100, Iterations = 50, setS = NULL, Relativities_s_sq = NULL, RequirePositiveWeights = FALSE, UseAllEstimators = TRUE, UseRcpp = TRUE, ShowProgressBar = TRUE, Bias_Correction_m_hat_infty = TRUE) {
  ###########################################################################################################
  # This function provides the generalized loss ratio method
  ###########################################################################################################
  #
  # Required Parameters:
  # --------------------
  # S = incremental triange
  # v = vector of volume estimates
  #
  # Optional Parameters:
  # --------------------
  # K = subset of {1,...,n-1} such that S_{i,k} is lognormally distributed for all (i,k) with k\in K;
  #       if K = NULL then the maximal K is used
  # c = (n-1)-dimensional vector such that eps with epsilon[i] = c[i] * eps for all i=1,..n;
  #       if c = NULL then c = (1,...,1) is used
  # eps_start = start value for eps in the recursion
  # lower_bound_eps = lower bound for eps in the recursion
  # upper_bound_eps = upper bound for eps in the recursion
  # eps_ext = external estimator for eps
  # s_sq_ext = exernal estimators for s^2, must be NULL (i.e. no external estimators given) or a vector of length n
  # lower_bound_cv_s = n-dimensional vector containing lower bounds for the coefficients of variantiosn, i.e.
  #       lower_bound_cv_s[k] = lower bound for the coefficient of variance of S_{i,k} in the recursion
  # upper_bound_cv_s = n-dimensional vector containing upper bounds for the coefficients of variantiosn, i.e.
  #       upper_bound_cv_s[k] = upper bound for the coefficient of variance of S_{i,k} in the recursion
  # Iterations = number of iterations in the recursion (default = 50)
  # setS = set of indices k for which we assume s_k^2 = s^2 * Relativities_s_sq with a given vector Relativities_s_sq;
  #         if NULL then no relativities are precsribed
  # Relativities_s_sq = vector containing the relativities;
  #         if NULL and setS != 0 then the relativities of \shat_k^2 (from the conventional LR method) are used
  # RequirePositiveWeights = Boolean; if TRUE than weights for estimators are forced to be positive (default = FALSE)
  # Bias_Correction_m_hat_infty = Boolean, it TRUE then m_hat_infty and h_hat_infty are multiplied with a factor to correct the bias of m_hat_infty
  #         (with notation of the paper: if FALSE, then \tildem_k and weights d are used, if TRUE then \hatm_k and weights h are used)
  #
  # Output:
  # -------
  # Predictions = Full triangle including predictions
  # v_hat_infty = vector containing the revised volumes estimates
  # m_hat_infty = vector containing the estimators for m_k
  # s_hat_sq_infty = vector containing the estimators for s_k^2
  # eps_hat_infty = estimator for eps (scalar); estimator for epsilon_i = c_i * eps_hat_infty
  # g_hat_infty = list containing the weight vectors \hatg_i^{(\infty)}
  # h_hat_infty = list containing the weight vectors \hath_i^{(\infty)}
  # v_hat = matrix containing the revised volume estimates after each iteration
  # m_hat = matrix containing the estimators for m_k after each iteration
  # s_hat_sq = matrix containing the estimators for s_k^2 after each iteration
  # eps_hat = vector containing the estimator for eps after each iteration
  # K = used set K
  # J[[i]] = used list of indices for the estimators of v_i^{-1}; for \hatr_{i,(j,k)}
  #         - J[[i]][,1] = i, i.e. estimator for v_i^{-1}
  #         - J[[i]][,2] = j
  #         - J[[i]][,3] = k
  # J_all = all J[[i]] pasted together in one matrix
  # Precision_eps_hat = estimated deviation from \hatvarepsilon^{(\infty)} due to limited number of iterations
  # Precision_s_hat_sq = estimated deviation from (\hats^{(\infty)})^2 due to limited number of iterations

  n <- length(v)
  if (is.null(K)) {
    K <- 1:(n-1)
  }
  if (is.null(c)) {
    c <- rep(1,n-1)
  }

  if (!is.null(eps_ext)) {
    lower_bound_eps <- eps_ext
    upper_bound_eps <- eps_ext
    eps_start <- eps_ext
  }


  # calculate m for conventional LR method
  m_start <- numeric()
  for (k in 1:n) {
    m_start <- c(m_start, sum(S[1:(n-k+1),k])/sum(v[1:(n-k+1)]))
  }
  # calculate s_squared for conventional LR method
  s_squared_start <- numeric()
  for (k in 1:(n-1)) {
    s_squared_start <- c(s_squared_start, sum((v[1:(n-k+1)]*(S[1:(n-k+1),k]/v[1:(n-k+1)]-m_start[k])^2))/(n-k))
  }
  if (s_squared_start[n-2] > 0) {
    s_squared_start <- c(s_squared_start, min(s_squared_start[n-1]^2/s_squared_start[n-2],s_squared_start[n-2]))
  } else {
    s_squared_start <- c(s_squared_start, 0)
  }


  if (!is.null(s_sq_ext)) {
    lower_bound_cv_S <- 0
    upper_bound_cv_S <- Inf
    s_squared_start <- s_sq_ext
  }

  # Check and reduce K if necessary (if not all incremental claims are positive)
  K <- CheckK(S, K)
  # Matrix J[[i]] contains the set $J_i$
  J <- CreateJ(S, K, UseAllEstimators)
  # Matrix containing the set $J$, N[i] = $#J_i$
  J_all <- J[[1]]
  N <- nrow(J[[1]])
  for (i in 2:n) {
    J_all <- rbind(J_all,J[[i]])
    N <- c(N,nrow(J[[i]]))
  }


  v_hat <- matrix(NA, ncol = Iterations, nrow = n)
  m_hat <- matrix(NA, ncol = Iterations, nrow = n)
  s_hat_sq <- matrix(NA, ncol = Iterations, nrow = n)
  eps_hat <- rep(NA, Iterations)

  v_hat[,1] <- v
  m_hat[,1] <- m_start
  s_hat_sq[,1] <- s_squared_start
  eps_hat[1] <- eps_start


  if (ShowProgressBar) {
    print(paste0("Recursion for parameter estimation (",Iterations," Iterations)"))
    pb <- utils::txtProgressBar(min = 0, max = 1, style = 3)
    progress <- 1 / Iterations
    utils::setTxtProgressBar(pb, progress)
  }

  Scale_Difference_v_v_hat <- 1



  for (nu in 2:Iterations) {
    Param <- ParameterEstimation(K, J, J_all, N, eps_hat[nu-1], c, S, m_hat[,nu-1], s_hat_sq[,nu-1], v_hat[,nu-1], v, lower_bound_eps, lower_bound_cv_S, upper_bound_eps, upper_bound_cv_S, setS = setS, Relativities_s_sq = Relativities_s_sq, s_squared_start=s_squared_start, s_sq_ext = s_sq_ext, RequirePositiveWeights = RequirePositiveWeights, UseRcpp = UseRcpp, Scale_Difference_v_v_hat = Scale_Difference_v_v_hat, Bias_Correction_m_hat_infty = Bias_Correction_m_hat_infty)
    if (adjust_v_hat_infty_for_small_eps_ext) {
      Scale_Difference_v_v_hat <- Param$Ratio_eps_cap
    }
    v_hat[,nu] <- Param$v_new
    m_hat[,nu] <- Param$m_new
    s_hat_sq[,nu] <- Param$s_sq_new
    eps_hat[nu] <- Param$eps_new
    h_hat_infty <- Param$h_new
    g_hat_infty <- Param$g_new
    if (ShowProgressBar) {
      progress <- nu / Iterations
      utils::setTxtProgressBar(pb, progress)
    }
  }
  if (ShowProgressBar) {close(pb)}

  if (eps_hat[Iterations]-eps_hat[Iterations-1] == 0) {
    Precision_eps_hat <- 0
  } else {
    a <-abs((eps_hat[Iterations]-eps_hat[Iterations-1])/(eps_hat[Iterations-1]-eps_hat[Iterations-2]))
    if (a<1) {
      Precision_eps_hat <-abs(eps_hat[Iterations]-eps_hat[Iterations-1])/(1-a)
    } else {
      Precision_eps_hat <- 1/0
    }
  }

  a <- abs(s_hat_sq[,Iterations]-s_hat_sq[,Iterations-1]) / abs(s_hat_sq[,Iterations-1]-s_hat_sq[,Iterations-2])
  Precision_s_hat_sq <- abs(s_hat_sq[,Iterations]-s_hat_sq[,Iterations-1])/(1-a)
  Precision_s_hat_sq <- ifelse(s_hat_sq[,Iterations]-s_hat_sq[,Iterations-1]==0, 0, Precision_s_hat_sq)

  v_hat_infty <- v_hat[,Iterations]
  m_hat_infty <- m_hat[,Iterations]
  s_hat_sq_infty <- s_hat_sq[,Iterations]
  eps_hat_infty <- eps_hat[Iterations]

  Predictions <- S
  for (k in 2:n) {
    for (i in (n-k+2):n) {
      Predictions[i,k] <- v_hat_infty[i] * m_hat_infty[k]
    }
  }



  Rslt <- list(Predictions = Predictions, v_hat_infty = v_hat_infty, m_hat_infty = m_hat_infty, s_hat_sq_infty = s_hat_sq_infty, eps_hat_infty = eps_hat_infty, g_hat_infty = g_hat_infty, h_hat_infty = h_hat_infty, v_hat = v_hat, m_hat = m_hat, s_hat_sq = s_hat_sq, eps_hat = eps_hat, K = K, J = J, J_all = J_all, Precision_eps_hat = Precision_eps_hat, Precision_s_hat_sq = Precision_s_hat_sq)
}



ParameterEstimation <- function(K, J, J_all, N, eps, c, S, m, s_squared, v, v_hat_original, lower_bound_eps = 0.001, lower_bound_cv_S = 0.001, upper_bound_eps = 0.1, upper_bound_cv_S = 100, setS = NULL, Relativities_s_sq = NULL, s_squared_start, s_sq_ext, RequirePositiveWeights = F, UseRcpp = T, Scale_Difference_v_v_hat = 1, Bias_Correction_m_hat_infty = T) {
  ###########################################################################################################
  # This function provides iteration for the recursive parameter estimation
  # The function is used by the function Generalized_LR_Method
  ###########################################################################################################
  #
  # Required Parameters:
  # --------------------
  # K = subset of {1,...,n-1} such that S_{i,k} is lognormally distributed for all (i,k) with k\in K
  # J[[i]] = used list of indices for the estimators of v_i^{-1}; for \hatr_{i,(j,k)}
  #         - J[[i]][,1] = i, i.e. estimator for v_i^{-1}
  #         - J[[i]][,2] = j
  #         - J[[i]][,3] = k
  # J_all = all J[[i]] pasted together in one matrix
  # N [i] = number of rows of J[[i]] (= N_i in the paper)
  # eps = estimator for eps before the iteration
  # c = (n-1)-dimensional vector such that eps with epsilon[i] = c[i] * eps for all i=1,..n
  # S = incremental triange
  # m = vector containing the estimators for m_k before the iteration
  # s_squared = vector containing the estimators for s_k^2 before the iteration
  # v = vector containing the estimators for v stemming from the previous iteration
  # v_hat_original = vector of volume estimates \hatv_i
  #
  # Optional Parameters:
  # --------------------
  # lower_bound_eps = lower bound for eps in the recursion
  # lower_bound_cv_s = n-dimensional vector containing lower bounds for the coefficients of variantiosn, i.e.
  #       lower_bound_cv_s[k] = lower bound for the coefficient of variance of S_{i,k} in the recursion
  # upper_bound_eps = upper bound for eps in the recursion
  # upper_bound_cv_s = n-dimensional vector containing upper bounds for the coefficients of variantiosn, i.e.
  #       upper_bound_cv_s[k] = upper bound for the coefficient of variance of S_{i,k} in the recursion
  # setS = set of indices k for which we assume s_k^2 = s^2 * Relativities_s_sq with a given vector Relativities_s_sq;
  #         if NULL then no relativities are precsribed
  # Relativities_s_sq = vector containing the relativities;
  #         if NULL and setS != NULL then the relativities of \shat_k^2 (from the conventional LR method) are used
  # s_squared_start
  # s_sq_ext = exernal estimators for s^2, must be NULL (i.e. no external estimators given) or a vector of length n
  # RequirePositiveWeights = Boolean; if TRUE than weights for estimators are forced to be positive (default = FALSE)
  # Bias_Correction_m_hat_infty = Boolean, it TRUE then m_new and h_new are multiplied with a factor to correct the bias of m_new
  #         (with notation of the paper: if FALSE, then \tildem_m and weights d are used, if TRUE then \hatm_k and weights h are used)
  #
  # Output:
  # -------
  # eps_new = estimator for eps after the iteration
  # g_new = list containing weight vectors \hatg_i used in this iteration
  # v_new = vector containing the revised volume estimates the iteration
  # h_new = list containing weight vectors \hath_i used in this iteration
  # m_new = vector containing the estimators for m_k after the iteration
  # s_sq_new = vector containing the estimators for s_k^2 after the iteration

  epsilon <- eps * c
  n <- length(v_hat_original)

  # calculate gammas
  gam <- matrix(0, ncol = n, nrow = n)
  for (i in 1:n) {
    for (k in K) {
      gam[i,k] <- 1+ s_squared[k]/(v[i]*m[k]^2)
    }
  }




  if (UseRcpp) {
    J_s <- c(J_all)
    CovR_vector <- CreateR_cpp(J_s, N, K, epsilon, m, s_squared, v, n)
    Dimensions <- numeric(n-1)
    k = 0
    for (i in 1:(n-1)) {
      k <- k +  N[i] * N[i]
      Dimensions[i] = k
    }
    CovR <- list()
    CovR[[1]] <- matrix(CovR_vector[1:Dimensions[1]], nrow = N[1])
    for (i in 2:(n-1)) {
      CovR[[i]] <- matrix(CovR_vector[(Dimensions[i-1]+1):Dimensions[i]], nrow = N[i])
    }
  } else {
    CovR <- CreateR_fast(J, K, epsilon, m, s_squared, v)
  }





  # Calculate estimators r_hat
  r_hat <- list()
  for (i in 1:n) {
    r_hat[[i]] <- numeric()
    for (nu in 1:N[i]) {
      j <- J[[i]][nu,2]
      k <- J[[i]][nu,3]
      if (j==0) {
        Temp <- 1/v_hat_original[i]
      } else {
        Temp <- 1/v_hat_original[j] * S[j,k] / S[i,k] / gam[i,k]
      }
      r_hat[[i]] <- c(r_hat[[i]], Temp)
    }
  }


  g <- list()
  for (i in 1:(n-1)) {
    g[[i]] <- CalculateOptimalWeights(CovR[[i]], RequirePositiveWeights = RequirePositiveWeights)
  }
  g[[n]] <- 1

  v_hat <- numeric()
  for (i in 1:n) {
    Temp <- t(g[[i]]) %*% r_hat[[i]]
    v_hat <- c(v_hat, 1/Temp)
  }

  if (Scale_Difference_v_v_hat < 1) {
    v_hat <- 1 / (Scale_Difference_v_v_hat / v_hat + (1 - Scale_Difference_v_v_hat) / v_hat_original)
  }

  m_hat_per_cell <- S / v_hat
  m_hat_vec <- list()
  for (k in 1:n) {
    m_hat_vec[[k]] <- m_hat_per_cell[1:(n-k+1),k]
  }

  a_per_cell <- matrix(0,ncol = n, nrow = n)
  for (l in K) {                    # for $l\not\in K$: $a_{i,l}=0$
    for (i in 1:(min(n-l+1,n-1))) { # $a_{1,n}$ = 0 =>  1:(min(n-l+1,n-1))) instead of 1:(n-l+1)
      a_per_cell[i,l] <- sum(g[[i]][J[[i]][,3]==l]) * (1/gam[i,l] - 1)
    }
  }

  a <- list()
  b <- list()
  for (l in 1:n) {
    a[[l]] <- a_per_cell[1:(n-l+1),l]
    b[[l]] <- m[l] * a[[l]]
  }

  if (UseRcpp) {
    # J_all <- J[[1]]
    # for (i in 2:n) {
    #   J_all <- rbind(J_all,J[[i]])
    # }
    J_s <- c(J_all)
    g_s <- numeric()
    for (i in 1:n) {
      g_s <- c(g_s,g[[i]])
    }
    CovM_vector <- CreateM_cpp(J_s, N, K, g_s, epsilon, m, s_squared, v, n)
    Dimensions <- numeric(n)
    k <- 0
    for (l in 1:n) {
      k <- k+(n-l+1)^2
      Dimensions[l] = k
    }
    CovM <- list()
    CovM[[1]] <- matrix(CovM_vector[1:Dimensions[1]], ncol = n)
    for (l in 2:n) {
      CovM[[l]] <- matrix(CovM_vector[(Dimensions[l-1]+1):Dimensions[l]], ncol = n-l+1)
    }
  } else {
    CovM <- CreateM_fast(J_all, N, K, g, epsilon, m, s_squared, v)
  }

  h <- list()
  for (l in 1:n) {
    if (det(CovM[[l]]) <= 0) {
      temp <- nrow(CovM[[l]])
      h[[l]] <- rep(1/temp, temp)
    } else {
      h[[l]] <- CalculateOptimalWeights((CovM[[l]] + b[[l]] %*% t(b[[l]])), RequirePositiveWeights = RequirePositiveWeights)
    }
  }

  m_hat <- numeric()
  for (k in 1:n) {
    Temp <- t(h[[k]]) %*% m_hat_vec[[k]]
    m_hat <- c(m_hat, Temp)
  }

  # calculate lower bounds for s_k:
  lower_bound_s_sq <- min(v) * lower_bound_cv_S^2 * m^2
  upper_bound_s_sq <- max(v) * upper_bound_cv_S^2 * m^2

  s_hat_sq <- numeric()
  for (k in 1:(n-1)) {
    Temp <- 0
    e <- diag(nrow = n-k+1)
    for (i in 1:(n-k+1)) {
      Temp2 <- t(e[,i]-h[[k]]) %*% (CovM[[k]] + b[[k]] %*% t(b[[k]])) %*% (e[,i]-h[[k]])
      Temp <- Temp + (m_hat_per_cell[i,k]-m_hat[k])^2 / Temp2
    }
    s_hat_sq[k] <- Temp * s_squared[k] / (n-k+1)
    if (is.nan(s_hat_sq[k])) s_hat_sq[k] <- 0
  }
  # s_hat_sq[n] <- min(1, s_hat_sq[n-1] / s_hat_sq[n-2]) * s_hat_sq[n-1]
  if (s_hat_sq[n-2] > 0) {
    s_hat_sq[n] <- min(s_hat_sq[n-1]^2 / s_hat_sq[n-2], s_hat_sq[n-2])
  } else {
    s_hat_sq[n] <- 0
  }

  if (!is.null(setS)) {
    if (is.null(Relativities_s_sq)) {
      Relativities_s_sq <- s_squared_start
    }

    Temp <- 0
    s_sq_relativities <- 0
    for (k in setS) {
      if (Relativities_s_sq[k] > 0) {
        Temp <- Temp + n - k
        s_sq_relativities <- s_sq_relativities + (n-k) / Relativities_s_sq[k] * s_hat_sq[k]
      }
    }
    s_sq_relativities <- s_sq_relativities / Temp
    for (k in setS) {
      s_hat_sq[k] <- s_sq_relativities * Relativities_s_sq[k]
    }
    if (s_hat_sq[n-2] > 0) {
      s_hat_sq[n] <- min(s_hat_sq[n-1]^2 / s_hat_sq[n-2], s_hat_sq[n-2])
    } else {
      s_hat_sq[n] <- 0
    }
  }

  s_hat_sq <- pmax(s_hat_sq, lower_bound_s_sq)
  s_hat_sq <- pmin(s_hat_sq, upper_bound_s_sq)

  if (!is.null(s_sq_ext)) {
    s_hat_sq = s_sq_ext
  }

  if (UseRcpp) {
    J_s <- c(J_all)
    CovQ_vector <- CreateQ_cpp(J_s, N, K, epsilon, m, s_squared, v, n)
    CovQ <- list()
    Dimensions <- numeric(n-1)
    k = 0
    for (i in 1:(n-1)) {
      k <- k +  N[i] * N[i+1]
      Dimensions[i] = k
    }
    CovQ <- list()
    CovQ[[1]] <- matrix(CovQ_vector[1:Dimensions[1]], nrow = N[1])
    for (i in 2:(n-1)) {
      CovQ[[i]] <- matrix(CovQ_vector[(Dimensions[i-1]+1):Dimensions[i]], nrow = N[i])
    }


  } else{
    CovQ <- CreateQ_fast(J, K, epsilon, m, s_squared, v)
    # CovQ <- CreateQ(J, K, epsilon, m, s_squared, v)
  }


  Nenner <- rep(0,n-1)
  Zaehler <- rep(0,n-1)

  for(i in 1:(n-1)) {
    e_1_len_i <- rep(0,length(g[[i]]))
    e_1_len_i[1] <- 1
    e_1_len_ip1 <- rep(0,length(g[[i+1]]))
    e_1_len_ip1[1] <- 1
    #Nenner[i] <- v_hat[i]^2 * t(e_1_len_i-g[[i]]) %*% CovR[[i]] %*% (e_1_len_i-g[[i]])
    Nenner[i] <- v[i]^2 * t(e_1_len_i-g[[i]]) %*% CovR[[i]] %*% (e_1_len_i-g[[i]])
    if (i < n-1) {
      # Nenner[i] <- Nenner[i] - 2 * v_hat[i] * v_hat[i+1] * t(e_1_len_i-g[[i]]) %*% CovQ[[i]] %*% (e_1_len_ip1-g[[i+1]])
      # Nenner[i] <- Nenner[i] + v_hat[i+1]^2 * t(e_1_len_ip1-g[[i+1]]) %*% CovR[[i+1]] %*% (e_1_len_ip1-g[[i+1]])
      Nenner[i] <- Nenner[i] - 2 * v[i] * v[i+1] * t(e_1_len_i-g[[i]]) %*% CovQ[[i]] %*% (e_1_len_ip1-g[[i+1]])
      Nenner[i] <- Nenner[i] + v[i+1]^2 * t(e_1_len_ip1-g[[i+1]]) %*% CovR[[i+1]] %*% (e_1_len_ip1-g[[i+1]])
    }

    # Zaehler[i] <-  v_hat[i]^2 * (1/v[i]-1/v_hat[i])^2
    # if (i < n-1) {
    #   Zaehler[i] <- Zaehler[i] - 2 * v_hat[i] * v_hat[i+1] * (1/v[i]-1/v_hat[i]) * (1/v[i+1]-1/v_hat[i+1])
    #   Zaehler[i] <- Zaehler[i] + v_hat[i+1]^2 * (1/v[i+1]-1/v_hat[i+1])^2
    # }
    Zaehler[i] <- ( v[i] * (1/v_hat_original[i]-1/v_hat[i]) - v[i+1] * (1/v_hat_original[i+1]-1/v_hat[i+1]) )^2
  }

  eps_hat_sq <- sum(eps^2 * Zaehler / Nenner) / (n-1)

  # apply upper and lower bound for eps:
  if (upper_bound_eps^2 < eps_hat_sq) {
    Ratio_eps_cap <- upper_bound_eps / sqrt(eps_hat_sq)
  } else {
    Ratio_eps_cap <- 1
  }
  eps_hat_sq <- max(eps_hat_sq, lower_bound_eps^2)
  eps_hat_sq <- min(eps_hat_sq, upper_bound_eps^2)
  eps_hat <- sqrt(eps_hat_sq)
  epsilon_hat <- c*eps_hat

  # m_hat and h[[k]] correspond to \mtilde_k^{(\nu+1)} and d_k^{(\nu+1)} in the paper; we apply the bias correction factor to obtain
  # \mhat_k^{(\nu+1)} and h_k^{(\nu+1)}
  if (Bias_Correction_m_hat_infty) {
    for (k in 1:n) {
      Bias_Correction_Factor <- 1 / (1 + sum(h[[k]]*a[[k]]))
      m_hat[k] <- m_hat[k] * Bias_Correction_Factor
      h[[k]] <- h[[k]] * Bias_Correction_Factor
    }

  }

  Rslt <- list(eps_new = eps_hat, g_new = g, v_new = v_hat, h_new = h, m_new = m_hat, s_sq_new = s_hat_sq, Ratio_eps_cap = Ratio_eps_cap)

}





CheckK <- function(DecumulatedTriangle, K) {
  n <- ncol(DecumulatedTriangle)
  UsedColumns <- numeric(0)
  for (k in K) {
    if (min(DecumulatedTriangle[1:(n-k+1),k])>0) {UsedColumns <- c(UsedColumns,k)}
  }
  return(UsedColumns)
}




CreateJ <- function(DecumulatedTriangle, K, UseAllEstimators) {
  # 1st column:         i means estimator for 1/v_i
  # 2nd & 3rd column: (0,0) means estimator 1/\hat{v}_i
  #                   (j,k) means estimator const * S_{j,k} / (S_{i,k} * \hat{v}_j)
  n <- ncol(DecumulatedTriangle)
  UsedColumns <- numeric(0)
  for (k in K) {
    if (min(DecumulatedTriangle[1:(n-k+1),k])>0) {UsedColumns <- c(UsedColumns,k)}
  }
  J <- list()
  for (i in 1:(n-1)) {
    J[[i]] <- numeric(0)
    J[[i]] <- c(J[[i]],i,0,0) # estimator 1/\hat{v}_i
    for (k in UsedColumns[UsedColumns<=n-i+1]) {
      for (j in (1:(n-k+1))[-i]) {
        if (UseAllEstimators || j > i) {
          J[[i]] <- c(J[[i]],i,j,k)
        }
      }
    }
    J[[i]] <- matrix(J[[i]], ncol = 3, byrow = T)
  }
  J[[n]] <- matrix(c(n,0,0), ncol = 3)
  return(J)
}




CreateR <- function(J, K, epsilon, m, s_squared, v) {

  n<- length(v)
  epsilon_squared <- epsilon^2
  gam <- matrix(0, ncol = n, nrow = n)
  for (i in 1:n) {
    for (k in K) {
      gam[i,k] <- 1+ s_squared[k]/(v[i]*m[k]^2)
    }
  }

  prod1plusepssq <- rep(1,n)
  for (i in 1:(n-1)) {
    prod1plusepssq[i] <- prod(1+epsilon_squared[i:(n-1)])
  }

  R <- list()
  for(i in 1:(n-1)) {
    NumberOfEstimators <- nrow(J[[i]])
    CovMatrix <- matrix(0, ncol = NumberOfEstimators, nrow = NumberOfEstimators)

    for (nu in 1:NumberOfEstimators) {
      for (mu in nu:NumberOfEstimators) {
        i_1 <- J[[i]][nu, 1]
        i_2 <- J[[i]][mu, 1]
        j_1 <- J[[i]][nu, 2]
        j_2 <- J[[i]][mu, 2]
        k_1 <- J[[i]][nu, 3]
        k_2 <- J[[i]][mu, 3]

        Temp <- 1

        chi <- matrix(c(j_1,k_1,1,i_1,k_1,-1,j_2,k_2,1,i_2,k_2,-1), ncol = 3, byrow = T)
        chi <- as.matrix(stats::aggregate(x=chi[,3], by = list(chi[,1],chi[,2]), FUN="sum"))
        chi <- chi[chi[,2]!=0 & chi[,3]!=0,]
        chi <- matrix(c(chi),ncol = 3)

        if (nrow(chi)>0) {
          rho <- chi
          rho[,3] <- 0.5*chi[,3]*(chi[,3]-1)
          for (l in 1:nrow(rho)) {
            Temp <- Temp * gam[rho[l,1],rho[l,2]]^rho[l,3]
          }
        }

        xi <- matrix(c(i_1,k_1,-1,i_2,k_2,-1), ncol = 3, byrow = T)
        xi <- as.matrix(stats::aggregate(x=xi[,3], by = list(xi[,1],xi[,2]), FUN="sum"))
        xi <- xi[xi[,2]!=0 & xi[,3]!=0,]
        xi <- matrix(c(xi),ncol = 3)

        if (nrow(xi)>0) {
          for (l in 1:nrow(xi)) {
            Temp <- Temp * gam[xi[l,1],xi[l,2]]^xi[l,3]
          }
        }

        if (j_1 == 0) {
          pi_1 <- i_1
        } else {
          pi_1 <- j_1
        }
        if (j_2 == 0) {
          pi_2 <- i_2
        } else {
          pi_2 <- j_2
        }

        CovMatrix[mu,nu] <- 1/v[i]^2 * (prod1plusepssq[max(pi_1,pi_2)] * Temp - 1)
      }
    }
    D <- diag(diag(CovMatrix))
    CovMatrix <- CovMatrix + t(CovMatrix) - D
    R[[i]] <- CovMatrix
  }
  return(R)
}




CreateR_fast <- function(J, K, epsilon, m, s_squared, v) { #(IndexMatrix, epsilon, m, s_sqared, v)
  n<- length(v)
  epsilon_squared <- epsilon^2
  gam <- matrix(0, ncol = n, nrow = n)
  for (i in 1:n) {
    for (k in K) {
      gam[i,k] <- 1+ s_squared[k]/(v[i]*m[k]^2)
    }
  }

  prod1plusepssq <- rep(1,n)
  for (i in 1:(n-1)) {
    prod1plusepssq[i] <- prod(1+epsilon_squared[i:(n-1)])
  }

  R <- list()

  for(i in 1:(n-1)) {
    NumberOfEstimators <- nrow(J[[i]])
    CovMatrix <- matrix(0, ncol = NumberOfEstimators, nrow = NumberOfEstimators)

    for (nu in 1:NumberOfEstimators) {
      for (mu in nu:NumberOfEstimators) {
        i_1 <- J[[i]][nu, 1]
        i_2 <- J[[i]][mu, 1]
        j_1 <- J[[i]][nu, 2]
        j_2 <- J[[i]][mu, 2]
        k_1 <- J[[i]][nu, 3]
        k_2 <- J[[i]][mu, 3]
        if (j_1 == 0 && j_2 == 0) {
          # case j_1 = j_2 = k_1 = k_2 = 0
          CovMatrix[nu,mu] <- 1/v[i_1]^2 * (prod1plusepssq[i_1]-1)
        } else if (j_1 == 0) {
          # case j_1 = k_1 = 0, j_2 <> 0, k_2 <> 0
          CovMatrix[nu,mu] <- 1/v[i_1]^2 * (prod1plusepssq[max(i_1,j_2)]-1)
        } else if (j_2 == 0) {
          # case j_1 <> 0, k_1 <> 0, j_2 = k_2 = 0
          CovMatrix[nu,mu] <- 1/v[i_1]^2 * (prod1plusepssq[max(j_1,i_2)]-1)
        } else if (j_1 != j_2 && k_1 != k_2) {
          # case j_i <> 0, k_i <> 0, j_1 <> j_2, k_1 <> k_2
          CovMatrix[nu,mu] <- 1/v[i_1]^2 * (prod1plusepssq[max(j_1,j_2)]-1)
        } else if (j_1 == j_2 && k_1 != k_2) {
          # case j_i <> 0, k_i <> 0, j_1 = j_2, k_1 <> k_2
          CovMatrix[nu,mu] <- 1/v[i_1]^2 * (prod1plusepssq[j_1]-1)
        } else if (j_1 != j_2 && k_1 == k_2) {
          # case j_i <> 0, k_i <> 0, j_1 <> j_2, k_1 = k_2
          CovMatrix[nu,mu] <- 1/v[i_1]^2 * (prod1plusepssq[max(j_1,j_2)]*gam[i_1,k_1]-1)
        } else if (j_1 == j_2 && k_1 == k_2) {
          # case j_i <> 0, k_i <> 0, j_1 = j_2, k_1 = k_2
          CovMatrix[nu,mu] <- 1/v[i_1]^2 * (prod1plusepssq[j_1]*gam[j_1,k_1]*gam[i_1,k_1]-1)
        } else{
          stop("Formula missing for i_1=i_2.")
        }
      }
    }
    D <- diag(diag(CovMatrix))
    CovMatrix <- CovMatrix + t(CovMatrix) - D
    R[[i]] <- CovMatrix
  }
  return(R)
}




CreateQ <- function(J, K, epsilon, m, s_squared, v) {

  n<- length(v)
  epsilon_squared <- epsilon^2
  gam <- matrix(0, ncol = n, nrow = n)
  for (i in 1:n) {
    for (k in K) {
      gam[i,k] <- 1+ s_squared[k]/(v[i]*m[k]^2)
    }
  }

  prod1plusepssq <- rep(1,n)
  for (i in 1:(n-1)) {
    prod1plusepssq[i] <- prod(1+epsilon_squared[i:(n-1)])
  }

  Q <- list()
  for(i in 1:(n-2)) {
    NumberOfEstimators_i <- nrow(J[[i]])
    NumberOfEstimators_ip <- nrow(J[[i+1]])
    CovMatrix <- matrix(0, nrow = NumberOfEstimators_i, ncol = NumberOfEstimators_ip)

    for (nu in 1:NumberOfEstimators_i) {
      for (mu in 1:NumberOfEstimators_ip) {
        i_1 <- J[[i]][nu, 1]
        i_2 <- J[[i+1]][mu, 1]
        j_1 <- J[[i]][nu, 2]
        j_2 <- J[[i+1]][mu, 2]
        k_1 <- J[[i]][nu, 3]
        k_2 <- J[[i+1]][mu, 3]

        Temp <- 1

        chi <- matrix(c(j_1,k_1,1,i_1,k_1,-1,j_2,k_2,1,i_2,k_2,-1), ncol = 3, byrow = T)
        chi <- as.matrix(stats::aggregate(x=chi[,3], by = list(chi[,1],chi[,2]), FUN="sum"))
        chi <- chi[chi[,2]!=0 & chi[,3]!=0,]
        chi <- matrix(c(chi),ncol = 3)

        if (nrow(chi)>0) {
          rho <- chi
          rho[,3] <- 0.5*chi[,3]*(chi[,3]-1)
          for (l in 1:nrow(rho)) {
            Temp <- Temp * gam[rho[l,1],rho[l,2]]^rho[l,3]
          }
        }

        xi <- matrix(c(i_1,k_1,-1,i_2,k_2,-1), ncol = 3, byrow = T)
        xi <- as.matrix(stats::aggregate(x=xi[,3], by = list(xi[,1],xi[,2]), FUN="sum"))
        xi <- xi[xi[,2]!=0 & xi[,3]!=0,]
        xi <- matrix(c(xi),ncol = 3)

        if (nrow(xi)>0) {
          for (l in 1:nrow(xi)) {
            Temp <- Temp * gam[xi[l,1],xi[l,2]]^xi[l,3]
          }
        }

        if (j_1 == 0) {
          pi_1 <- i_1
        } else {
          pi_1 <- j_1
        }
        if (j_2 == 0) {
          pi_2 <- i_2
        } else {
          pi_2 <- j_2
        }

        CovMatrix[nu, mu] <- 1/(v[i]*v[i+1]) * (prod1plusepssq[max(pi_1,pi_2)] * Temp - 1)
      }
    }
    Q[[i]] <- CovMatrix
  }
  return(Q)
}





CreateQ_fast <- function(J, K, epsilon, m, s_squared, v) { #(IndexMatrix, epsilon, m, s_sqared, v)
  n<- length(v)
  epsilon_squared <- epsilon^2
  gam <- matrix(0, ncol = n, nrow = n)
  for (i in 1:n) {
    for (k in K) {
      gam[i,k] <- 1+ s_squared[k]/(v[i]*m[k]^2)
    }
  }

  prod1plusepssq <- rep(1,n)
  for (i in 1:(n-1)) {
    prod1plusepssq[i] <- prod(1+epsilon_squared[i:(n-1)])
  }

  Q <- list()

  for(i in 1:(n-2)) {
    NumberOfEstimators1 <- nrow(J[[i]])
    NumberOfEstimators2 <- nrow(J[[i+1]])
    CovMatrix <- matrix(0, nrow = NumberOfEstimators1, ncol = NumberOfEstimators2)

    for (nu in 1:NumberOfEstimators1) {
      for (mu in 1:NumberOfEstimators2) {
        i_1 <- J[[i]][nu, 1]
        i_2 <- J[[i+1]][mu, 1]
        j_1 <- J[[i]][nu, 2]
        j_2 <- J[[i+1]][mu, 2]
        k_1 <- J[[i]][nu, 3]
        k_2 <- J[[i+1]][mu, 3]

        # case i_2 = i_1 + 1
        if (j_1 == 0 && j_2 == 0) {
          # case j_1 = j_2 = k_1 = k_2 = 0
          CovMatrix[nu,mu] <- 1/(v[i_1]*v[i_2]) * (prod1plusepssq[max(i_1,i_2)]-1)
        } else if (j_1 == 0) {
          # case j_1 = k_1 = 0, j_2 <> 0, k_2 <> 0
          CovMatrix[nu,mu] <- 1/(v[i_1]*v[i_2]) * (prod1plusepssq[max(i_1,j_2)]-1)
        } else if (j_2 == 0) {
          # case j_1 <>0, k_1 <> 0, j_2 = k_2 = 0
          CovMatrix[nu,mu] <- 1/(v[i_1]*v[i_2]) * (prod1plusepssq[max(j_1,i_2)]-1)
        } else if (j_1 != j_2 && k_1 != k_2) {
          # case j_i <> 0, k_i <> 0, j_1 <> j_2, k_1 <> k_2
          CovMatrix[nu,mu] <- 1/(v[i_1]*v[i_2]) * (prod1plusepssq[max(j_1,j_2)]-1)
        } else if (j_1 == j_2 && k_1 != k_2) {
          # case j_i <> 0, k_i <> 0, j_1 = j_2, k_1 <> k_2
          CovMatrix[nu,mu] <- 1/(v[i_1]*v[i_2]) * (prod1plusepssq[j_1]-1)
        } else if (j_1 != j_2 && k_1 == k_2) {
          # case j_i <> 0, k_i <> 0, j_1 <> j_2, k_1 = k_2
          if (i_1 != j_2 && i_2 != j_1) {
            CovMatrix[nu,mu] <- 1/(v[i_1]*v[i_2]) * (prod1plusepssq[max(j_1,j_2)]-1)
          } else if (i_1 == j_2 && i_2 != j_1) {
            CovMatrix[nu,mu] <- 1/(v[i_1]*v[i_2]) * (prod1plusepssq[max(j_1,j_2)]/gam[i_1,k_1]-1)
          } else if (i_1 != j_2 && i_2 == j_1) {
            CovMatrix[nu,mu] <- 1/(v[i_1]*v[i_2]) * (prod1plusepssq[max(j_1,j_2)]/gam[i_2,k_1]-1)
          } else {
            CovMatrix[nu,mu] <- 1/(v[i_1]*v[i_2]) * (prod1plusepssq[max(j_1,j_2)]/(gam[i_1,k_1]*gam[i_2,k_1])-1)
          }
        } else if (j_1 == j_2 && k_1 == k_2) {
          # case j_i <> 0, k_i <> 0, j_1 = j_2, k_1 = k_2
          CovMatrix[nu,mu] <- 1/(v[i_1]*v[i_2]) * (prod1plusepssq[max(j_1,j_2)]*gam[j_1,k_1]-1)
        } else{
          stop("Formula missing for i_1!=i_2.")
        }
      }
    }
    Q[[i]] <- CovMatrix
  }
  NumberOfEstimators1 <- nrow(J[[n-1]])
  NumberOfEstimators2 <- 1
  Q[[n-1]] <- matrix(0, ncol = NumberOfEstimators2, nrow = NumberOfEstimators1)

  return(Q)
}





CalculateOptimalWeights <- function(CovMatrix, Bias = NULL, RequirePositiveWeights = F) {
  # CovMatrix = Covariance matrix of the estimators (positive definite)
  # Bias = vector with biases of the estimators

  if (sum(abs(CovMatrix-t(CovMatrix))) > 0.0000001*sum(abs(CovMatrix))) {stop("Covariance Matrix not symmetric.")}
  if (!is.null(Bias)) {
    if (length(Bias) != ncol(CovMatrix)) {stop("Dimensions of bias and covariance matrix not compatible.")}
  }
  if (det(CovMatrix)<=0) {stop("Covariance matrix not positive definite.")}

  NumberOfParameters <- ncol(CovMatrix)

  #CovMatrix_inv <- solve(CovMatrix)
  if (is.null(Bias)) {
    A <- CovMatrix
  } else {
    A <- CovMatrix + matrix(c(Bias), ncol = 1) %*% matrix(c(Bias), nrow = 1)
  }
  g <- solve(A) %*% rep(1,NumberOfParameters)
  g <- g / sum(g)
  if (RequirePositiveWeights == F || min(g) >= 0) {
    return(g)
  } else {
    PosNegativeWeights <- which(g<0)
    #ResultRecursion <- NULL
    #for (i in PosNegativeWeights) {
    #  ResultRecursion <- c(ResultRecursion, RecursionOptWeights(A, i))
    #}
    #for (i in 1:length(PosNegativeWeights)) {
    #  ResultRecursion <- c(ResultRecursion, RecursionOptWeights(A, PosNegativeWeights[-i]))
    #}
    ResultRecursion <- RecursionOptWeights(A, PosNegativeWeights)
    ResultRecursion <- matrix(ResultRecursion, byrow = TRUE, ncol = NumberOfParameters+1)

    #g <- matrix(ResultRecursion[ResultRecursion[,NumberOfParameters+1]==min(ResultRecursion[,NumberOfParameters+1])],nrow = NumberOfParameters+1)
    g <- c(ResultRecursion[ResultRecursion[,NumberOfParameters+1]==min(ResultRecursion[,NumberOfParameters+1])])
    #g <- c(g[1,1:NumberOfParameters])
    g <- matrix(g[1:NumberOfParameters],ncol=1)
    #return(ResultRecursion)
    return(g)
  }
}




RecursionOptWeights <- function(A, ExcludedParameters) {
  NumberOfParameters <- ncol(A)
  B <- A[-ExcludedParameters, -ExcludedParameters]
  UsedParameters <- (1:NumberOfParameters)[-ExcludedParameters]
  NumberOfUsedParameters <- length(UsedParameters)
  g <- solve(B) %*% rep(1,NumberOfUsedParameters)
  g <- g/sum(g)
  if (min(g)>=0) {
    Varianz <- t(g) %*% B %*% g
    ResultRecursion <- rep(0,NumberOfParameters)
    ResultRecursion[UsedParameters] <- g
    ResultRecursion <- c(ResultRecursion, Varianz)
    return(ResultRecursion)
  } else {
    h <- rep(0,NumberOfParameters)
    h[UsedParameters] <- g
    PosNegativeWeights <- which(h<0)
    ResultRecursion <- NULL
    #for (i in PosNegativeWeights) {
    #  ResultRecursion <- c(ResultRecursion, RecursionOptWeights(A, c(ExcludedParameters,i)))
    #}
    ResultRecursion <- c(ResultRecursion, RecursionOptWeights(A, c(ExcludedParameters,PosNegativeWeights)))
    return(ResultRecursion)
  }
}




CreateM <- function(J_all, N, K, g, epsilon, m, s_squared, v) {
  # provides for l=1,..,n the covariance martices of the Estimators S_{il} / \breve{v}_i (i=1,...,n-l+1)
  n <- length(v)
  N_cum <- cumsum(N)
  Position1v <- c(1,cumsum(N)+1)

  g_matrix <- matrix(0, nrow = N_cum[n], ncol = n)
  for (i in 1:n) {
    g_matrix[Position1v[i]:(Position1v[i+1]-1),i] <- g[[i]]
  }


  # calculate squared epsilons
  epsilon_sqared <- epsilon^2

  # calculate gammas
  # set gam[i,l] <- 1 if $l\not\in K$
  gam <- matrix(1, ncol = n, nrow = n)
  for (i in 1:n) {
    for (k in K) {
      gam[i,k] <- 1+ s_squared[k]/(v[i]*m[k]^2)
    }
  }

  # calculate \prod_{\nu=i}^{n-1} (1+\varepsilon_\nu^2)
  # set to 1 if i=n
  prod1plusepssq <- rep(1,n)
  for (i in 1:(n-1)) {
    prod1plusepssq[i] <- prod(1+epsilon_sqared[i:(n-1)])
  }


  # initialize list that contains the covariance matrices for the different l
  M  <- list()
  for (l in 1:n) {
    CovMatrix <- matrix(0, ncol = N_cum[n-l+1], nrow =  N_cum[n-l+1])

    for (nu in 1:N_cum[n-l+1]) {
      for (mu in nu:N_cum[n-l+1]) {

        i_1 <- J_all[nu, 1]
        i_2 <- J_all[mu, 1]
        j_1 <- J_all[nu, 2]
        j_2 <- J_all[mu, 2]
        k_1 <- J_all[nu, 3]
        k_2 <- J_all[mu, 3]


        Index_S1 <- matrix(c(i_1,l),ncol = 2, byrow = T)
        Index_S2 <- matrix(c(i_2,l),ncol = 2, byrow = T)
        Index_S12 <- matrix(c(i_1,l,i_2,l),ncol = 2, byrow = T)
        Index_r1 <- matrix(c(i_1,j_1,k_1),ncol = 3)
        Index_r2 <- matrix(c(i_2,j_2,k_2),ncol = 3)
        Index_r12 <- matrix(c(i_1,j_1,k_1,i_2,j_2,k_2),ncol = 3, byrow = T)

        E_1 <- Calculate_E_Sr(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v)
        E_2 <- Calculate_E_Sr(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v)
        E_12 <- Calculate_E_Sr(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v)

        CovMatrix[mu,nu] <- E_12-E_1*E_2
        #print(c(l, mu,nu))
      }
    }

    #print(CovMatrix[1:10,1:10])
    #stop()

    D <- diag(diag(CovMatrix))
    CovMatrix <- CovMatrix + t(CovMatrix) - D

    g_matrix_l <- matrix(g_matrix[1:N_cum[n-l+1],1:(n-l+1)],ncol = n-l+1)
    M[[l]] <- t(g_matrix_l) %*%  CovMatrix %*% g_matrix_l
  }
  return(M)
}




CreateM_fast <- function(J_all, N, K, g, epsilon, m, s_squared, v) {
  # provides for l=1,..,n the covariance martices of the Estimators S_{il} / \breve{v}_i (i=1,...,n-l+1)

  n <- length(v)
  N_cum <- cumsum(N)
  Position1v <- c(1,cumsum(N)+1)

  # add estimator S_n1/v_n

  # NumberOfEstimators = vector containing the numbers of estimators for 1/v_1, 1/v_2, ..., 1/v_{n-1}
  # add one estimator for 1/v_n (1/v_n itself!)
  #NumberOfEstimators <- c(NumberOfEstimators,1)
  #NumberOfEstimatorsCum <- cumsum(NumberOfEstimators)

  # add entry i=n, j=0, k=0
  #IndexMatrix <- rbind(IndexMatrix, c(n,0,0))

  # add weight g_{n,(0,0)}:=1 in g_matrix
  #g_matrix <- rbind(g_matrix,0)
  #g_matrix <- cbind(g_matrix,0)
  #g_matrix[nrow(g_matrix),ncol(g_matrix)] <- 1

  g_matrix <- matrix(0, nrow = N_cum[n], ncol = n)
  for (i in 1:n) {
    g_matrix[Position1v[i]:(Position1v[i+1]-1),i] <- g[[i]]
  }


  # calculate squared epsilons
  epsilon_sqared <- epsilon^2

  # calculate gammas
  gam <- matrix(NA, ncol = n, nrow = n)
  for (i in 1:n) {
    for (k in K) {
      gam[i,k] <- 1+ s_squared[k]/(v[i]*m[k]^2)
    }
  }
  #gam <- matrix(0, ncol = n, nrow = n)
  #for (i in 1:n) {
  #  for (k in 1:n) {
  #    gam[i,k] <- 1+ s_squared[k]/(v[i]*m[k]^2)
  #  }
  #}

  # calculate \prod_{\nu=i}^{n-1} (1+\varepsilon_\nu^2)
  # set to 1 if i=n
  prod1plusepssq <- rep(1,n)
  for (i in 1:(n-1)) {
    prod1plusepssq[i] <- prod(1+epsilon_sqared[i:(n-1)])
  }


  # initialize list that contains the covariance matrices for the different l
  M  <- list()
  for (l in 1:n) {
    CovMatrix <- matrix(0, ncol = N_cum[n-l+1], nrow =  N_cum[n-l+1])

    for (nu in 1:N_cum[n-l+1]) {
      for (mu in nu:N_cum[n-l+1]) {
        i_1 <- J_all[nu, 1]
        i_2 <- J_all[mu, 1]
        j_1 <- J_all[nu, 2]
        j_2 <- J_all[mu, 2]
        k_1 <- J_all[nu, 3]
        k_2 <- J_all[mu, 3]
        if (i_1 == i_2)  {
          # case i_1 = i_2
          if (j_1 == 0 && j_2 == 0) {
            # case j_1 = j_2 = k_1 = k_2 = 0
            CovMatrix[nu,mu] <- (m[l]^2 + s_squared[l]/v[i_1]) * prod1plusepssq[i_1] - m[l]^2
          } else if (j_1 == 0 && k_2 == l) {
            # case ok
            CovMatrix[nu,mu] <- m[l]^2/gam[i_1,l] * (prod1plusepssq[max(i_1,j_2)]-1)
          } else if (j_1 == 0 && k_2 != l) {
            # case ok
            CovMatrix[nu,mu] <- (m[l]^2 + s_squared[l]/v[i_1]) * prod1plusepssq[max(i_1,j_2)] - m[l]^2
          } else if (j_2 == 0 && k_1 == l) {
            # case ok
            CovMatrix[nu,mu] <- m[l]^2/gam[i_1,l] * (prod1plusepssq[max(i_1,j_1)]-1)
          } else if (j_2 == 0 && k_1 != l) {
            # case ok
            CovMatrix[nu,mu] <- (m[l]^2 + s_squared[l]/v[i_1]) * prod1plusepssq[max(i_1,j_1)] - m[l]^2
          } else if (j_1 == j_2 && k_1 == k_2 && k_1 == l) {
            # case ok
            CovMatrix[nu,mu] <- m[l]^2/gam[i_1,l]^2 * (gam[j_1,l]*prod1plusepssq[j_1]-1)
          } else if (j_1 == j_2 && k_1 == k_2 && k_1 != l) {
            # case ok
            CovMatrix[nu,mu] <- (m[l]^2 + s_squared[l]/v[i_1]) * gam[i_1,k_1]*gam[j_1,k_1]*prod1plusepssq[j_1] - m[l]^2
          } else if (j_1 != j_2 && k_1 == k_2 && k_1 != l) {
            # case ok
            CovMatrix[nu,mu] <- (m[l]^2 + s_squared[l]/v[i_1]) * gam[i_1,k_1]*prod1plusepssq[max(j_1,j_2)] - m[l]^2
          } else if (j_1 != j_2 && k_1 == k_2 && k_1 == l) {
            # case ok
            CovMatrix[nu,mu] <- m[l]^2/gam[i_1,l]^2 * (prod1plusepssq[max(j_1,j_2)]-1)
          } else if (k_1 != k_2 && k_1 != l && k_2 != l) {
            # case ok
            CovMatrix[nu,mu] <- (m[l]^2+s_squared[l]/v[i_1]) * prod1plusepssq[max(j_1,j_2)] - m[l]^2
          } else if (k_1 != k_2 && k_1 == l) {
            # case ok
            CovMatrix[nu,mu] <- m[l]^2/gam[i_1,l] * (prod1plusepssq[max(j_1,j_2)]-1)
          } else if (k_1 != k_2 && k_2 == l) {
            # case ok
            CovMatrix[nu,mu] <- m[l]^2/gam[i_1,l] * (prod1plusepssq[max(j_1,j_2)]-1)
          } else  {
            readline("Fehler1")
            readline(c(l,i_1,j_1,k_1,i_2,j_2,k_2))
          }
        } else {
          # case i_1 <> i_2
          if (j_1 == 0 && j_2 == 0) {
            # case
            CovMatrix[nu,mu] <- m[l]^2 * (prod1plusepssq[max(i_1,i_2)]-1)
          } else if (j_2 == 0 && j_1 == i_2 && k_1 == l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2/gam[i_1,l] * (gam[i_2,l] * prod1plusepssq[i_2]-1)
          } else if (j_1 == 0 && j_2 == i_1 && k_2 == l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2/gam[i_2,l] * (gam[i_1,l] * prod1plusepssq[i_1]-1)
          } else if (j_2 == 0 && j_1 != i_2 && k_1 == l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2/gam[i_1,l] * (prod1plusepssq[max(j_1,i_2)]-1)
          } else if (j_1 == 0 && j_2 != i_1 && k_2 == l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2/gam[i_2,l] * (prod1plusepssq[max(j_2,i_1)]-1)
          } else if (j_2 == 0 && k_1 != l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2 * (prod1plusepssq[max(j_1,i_2)]-1)
          } else if (j_1 == 0 && k_2 != l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2 * (prod1plusepssq[max(j_2,i_1)]-1)
          } else if (j_1 == j_2 && k_1 == k_2 && k_1 == l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2/(gam[i_1,l]*gam[i_2,l]) * (gam[j_1,l] * prod1plusepssq[j_1]-1)
          } else if (j_1 == j_2 && k_1 == k_2 && k_1 != l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2 * (gam[j_1,k_1] * prod1plusepssq[j_1]-1)
          } else if (j_1 != j_2 && j_1 != i_2 && j_2 != i_1 && k_1 == k_2 && k_1 != l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2 * (prod1plusepssq[max(j_1,j_2)]-1)
          } else if (j_1 != j_2 && j_1 == i_2 && j_2 == i_1 && k_1 == k_2 && k_1 != l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2 * (prod1plusepssq[max(j_1,j_2)]/(gam[i_1,k_1]*gam[i_2,k_2])-1)
          } else if (j_1 != j_2 && j_1 != i_2 && j_2 == i_1 && k_1 == k_2 && k_1 != l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2 * (prod1plusepssq[max(j_1,j_2)]/gam[i_1,k_1]-1)
          } else if (j_1 != j_2 && j_1 == i_2 && j_2 != i_1 && k_1 == k_2 && k_1 != l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2 * (prod1plusepssq[max(j_1,j_2)]/gam[i_2,k_2]-1)
          } else if (j_1 != j_2 && k_1 == k_2 && k_1 == l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2/(gam[i_1,l]*gam[i_2,l]) * (prod1plusepssq[max(j_1,j_2)]-1)
          } else if (j_1 != i_2 && k_1 == l && k_2 != l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2/gam[i_1,l] * (prod1plusepssq[max(j_1,j_2)]-1)
          } else if (j_2 != i_1 && k_2 == l && k_1 != l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2/gam[i_2,l] * (prod1plusepssq[max(j_1,j_2)]-1)
          } else if (j_1 == i_2 && k_1 == l && k_2 != l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2/gam[i_1,l] * (gam[i_2,l]*prod1plusepssq[max(j_1,j_2)]-1)
          } else if (j_2 == i_1 && k_2 == l && k_1 != l) {
            # case
            CovMatrix[nu,mu] <- m[l]^2/gam[i_2,l] * (gam[i_1,l]*prod1plusepssq[max(j_1,j_2)]-1)
          } else if (k_1 != l && k_2 != l && k_1 != k_2) {
            # case
            CovMatrix[nu,mu] <- m[l]^2 * (prod1plusepssq[max(j_1,j_2)]-1)
          } else {
            print(c(l,i_1,j_1,k_1,i_2,j_2,k_2))
            readline("Fehler2")

          }
        }
      }
    }
    D <- diag(diag(CovMatrix))
    CovMatrix <- CovMatrix + t(CovMatrix) - D

    g_matrix_l <- matrix(g_matrix[1:N_cum[n-l+1],1:(n-l+1)],ncol = n-l+1)
    M[[l]] <- t(g_matrix_l) %*%  CovMatrix %*% g_matrix_l
  }
  return(M)
}


