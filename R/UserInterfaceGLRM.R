#' This function provides the results of the Loss Ratio and the Gen. Loss Ratio Method in the Extended Additive Model
#' @param triangle Numeric matrix (n x n) containing the claims triangle
#' @param volumes Numeric vector of lenght n containing the volume estimates
#' @param is.incremental  Logical. Indicates whether the claims triangle is cumulative (default) or incremental.
#' @param K Numeric vector containing the columns with lognormally distributed increments. If K is set to NULL, then the maximum possible set is used.
#' @param c Vector of lenght n-1 containing the relativities of the epsilon_i. If c = NULL (defalut) then c = (1,...,1) is used, i.e. all epsilon_i are equal. c is normalized such that max(c) = 1.
#' @param eps_ext Numeric. External estimator for epsilon. If eps_ext = NULL (default) then epsilon is estimated in the recursion.
#' @param eps_min Numeric. Lower bound for the estimator for epsilon.
#' @param eps_max Numeric. Upper bound for the estimator for epsilon.
#' @param eps_start Numeric. Start value for epsilon in the recursion.
#' @param s_sq_ext Numeric vector. Contains external estimators for the \eqn{s_k^2}. If s_sq_ext = NULL (default) then epsilon is estimated in the recursion. Not used if eps_ext = 0.
#' @param Iterations Integer containing the number of iterations for the recursion.
#' @param UseAllEstimators Logical. If TRUE then all estimators r_{i,(j,k)} for 1/v_i are used. If FALSE then only those with (j,k)=(0,0) or j>i.
#' @param Weights_LR Either one of the strings 'optimal' or 'canonical' or a square matrix containing the weights to be used.
#' @param Calculate_MSE_GLR Logical. If FALSE then the prediction error of the generalized loss ratio method is not calculated (better performance)
#' @param Calculate_MSE_LR Logical. If FALSE then the prediction error of the loss ratio method is not calculated
#' @param NumberOfSimulations Numeric. Number of simulations for simulation of standard error.
#' @param RunAllSimulations Logical. If TRUE then also the time consuming simulations are run.
#' @param IterationsInSimMSE Integer. Number of iterations in the recursive parameter estimation in each simulation (if RunAllSimulations = TRUE).
#' @param Use_eps_ext_inSimMSE Logical. If TRUE then eps_ext is used in the simulation of the MSE with unfixed weights. If FALSE then eps_ext = NULL is used in the simulation.
#' @param Use_s_sq_ext_inSimMSE Logical. If TRUE then s_sq_ext is used in the simulation of the MSE with unfixed weights. If FALSE then s_sq_ext = NULL is used in the simulation.
#' @param UseRelativities_s_sq Logical. If TRUE (default) then relativities are used for the estimators of \eqn{s_k^2}
#' @param Relativities_s_sq Numeric vector of lenght n which contains the relativities if UseRelativities_s_sq = TRUE. If NULL then the relativities from the estimators of the loss ratio method are used.
#' @param SetIndicesRelativities. Numeric vector of lenght < n. Contains the columns k for which the relativities are applied. If NULL and UseRelativities_s_sq = TRUE, then SetIndicesRelativities is set to 1:(n-1).
#' @param Seed. If Seed != NULL then the value is used as seed for the simulations.
#' @param adjust_v_hat_infty_for_small_eps_ext. Logical. If TRUE then adjustment of revised volume estimates v_hat_infty is reduced if the adjustment is not plausible for a given eps_ext or eps_max.
#' @param Bias_Correction_m_hat_infty Logical. If TRUE (default) then the estimators m_hat_infty are corrected for bias (otherwise weights d are used instead of h)
#' @param UseRcpp If TRUE then fast Rcpp routines are used.

#' @return The function returns a list with various results:
#'
#' \itemize{
#' \item Summary_GLR: Overview to the results of the generalized loss ratio method in the extended additive model
#' \item Summary_LR: Overview to the results of the loss ratio method in the extended additive model
#' \item Predictions_GLR: Predictions of the generalized loss ratio method (cumulative)
#' \item Reserves_GLR: Reserves of the generalized loss ratio method
#' \item Predictions_LR: Predictions of the generalized loss ratio method (cumulative)
#' \item Reserves_LR: Reserves of the loss ratio method
#' \item Param_eps_hat_infty: Estimator for eps from the recursion.
#' \item Param_m_hat_infty: Estimator for m from the recursion.
#' \item Param_s_hat_sq_infty: Estimators for \eqn{s_k^2} from the recursion.
#' \item Param_v_hat_infty: Estimator for the volumes from the recursion.
#' \item Weights_LR: Matrix containing the used weights for the loss ratio method.
#' \item Weights_LR_Type: String containing the type of the chosen weights.
#' \item Recursion_eps_hat: Estimator for \eqn{\epsilon} for each step \eqn{k} of the recursion.
#' \item Recursion_v_hat: Estimators for the volumes \eqn{v_i} for each step of the recursion.
#' \item Recursion_m_hat: Estimators for \eqn{m_k} for each step of the recursion.
#' \item Recursion_s_hat_sq: Estimators for \eqn{s_k^2} for each step of the recursion.
#' \item se_GLR_fixed_weights_calc: Calculated standard errors of the generalized loss ratio method with fixed weights \eqn{g} and \eqn{h}.
#' \item se_GLR_fixed_weights_sim: Simulated standard errors of the generalized loss ratio method with fixed weights \eqn{g} and \eqn{h}.
#' \item se_GLR_sim: Simulated standard errors of the generalized loss ratio method, where the weights \eqn{g} and \eqn{h} are estimated in each simulation.
#' \item se_LR_fixed_weights_calc: Calculated standard errors of the generalized loss ratio method with fixed weights \eqn{w}.
#' \item se_LR_fixed_weights_sim: Simulated standard errors of the generalized loss ratio method with fixed weights \eqn{w}.
#' \item se_LR_sim: Simulated standard errors of the generalized loss ratio method, where the weights \eqn{w} are calculated in each simulation (in case of canonical or optimal weights).
#' \item se_CL_sim: Simulated standard errors of the chain ladder method (simulation based on the extended additive model).
#' }

#' @examples
#' GLRM(glrm_example1_C, glrm_example1_v_hat)
#' GLRM(glrm_example1_S, glrm_example1_v_hat, is.incremental = TRUE, K = 1:2, UseAllEstimators = FALSE,
#'      eps_ext = 0.03, Calculate_MSE_LR = TRUE, Calculate_MSE_GLR = TRUE, Weights_LR = "canonical")
#' GLRM(glrm_example2_C, glrm_example2_v_hat, Iterations = 100, RunAllSimulations = TRUE,
#'      NumberOfSimulations = 1000)

#' @export
GLRM <- function(triangle, volumes, is.incremental = FALSE, K = NULL, c = NULL, eps_ext = NULL, eps_start = NULL, eps_min = 0.001, eps_max = 0.1, s_sq_ext = NULL, Iterations = 50, UseAllEstimators = TRUE, Weights_LR = "optimal", Calculate_MSE_GLR = FALSE, Calculate_MSE_LR = FALSE, NumberOfSimulations = 10000, RunAllSimulations = FALSE, IterationsInSimMSE = 20, Use_eps_ext_inSimMSE = FALSE, Use_s_sq_ext_inSimMSE = FALSE, Seed = NULL, UseRcpp = TRUE, UseRelativities_s_sq = TRUE, Relativities_s_sq = NULL, SetIndicesRelativities = NULL, adjust_v_hat_infty_for_small_eps_ext = TRUE, Bias_Correction_m_hat_infty = TRUE) {
  if (!is.matrix(triangle)) {
    stop("Triangle must be a matrix.")
  } else if (nrow(triangle) != ncol(triangle)) {
    stop("Triangle must be a square matrix.")
  } else if (!is.numeric(triangle)) {
    stop("Triangle must be a numeric matrix.")
  }

  n <- nrow(triangle)
  if (is.integer(triangle)) {
    triangle <- as.double(triangle)
    triangle <- matrix(triangle, ncol = n)
  }
  if (!is.na(triangle[n,2])) {
    if ((is.incremental & triangle[n,2] != 0) | (!is.incremental & triangle[n,2] != triangle[n,1])) {
      warning("The triangle seems to contain data below the diagonal. This data is ignored!")
    }
  }

  if (is.incremental) {
    S <- triangle
    C <- i2c(S)
  } else {
    C <- triangle
    S <- c2i(C)
  }
  if (min(C[row(C)+col(C) <= n+1] <= 0)) {
    stop("Entries of cumulative triangle must be positive.")
  }




  if (!is.numeric(volumes)) {
    stop("Volumes must be provided in a vector!")
  } else {
    v <- c(volumes)
    if (length(volumes) != n) {
      stop("Vector of volumes has wrong size!")
    }
    if (min(volumes) <= 0) {
      stop("Volumes must be positive.")
    }
  }
  if (is.integer(volumes)) {
    volumes <- as.double(volumes)
  }

  if (is.null(c)) {
    c <- rep(1,n-1)
  } else {
    if (!is.numeric(c)) {
      stop("c must be a numeric vector.")
    }
    if (length(c) != n-1) {Dimensions[l-1]
      stop("Vector c has wrong size.")
    }
    if (min(c) <= 0) {
      stop("Vector c must have positive entries.")
    }
    c <- c / max(c)
  }
  if (is.null(K)) {
    K <- 1:(n-1)
  } else {
    if (!is.numeric(K)) {
      stop("K must be a integer vector.")
    }
    if (sum(abs(K - as.integer(K))) > 0) {
      stop("K must be a integer vector.")
    }
    if (!(1 %in% K)) {
      stop("Vector K must contain 1.")
    }
    if (min(K==sort(K))==0) {
      stop("K has to be sorted.")
    }
    if (min(K) < 1 || max(K) > n-1) {
      stop("Vector K must have entries between 1 and n-1.")
    }
  }
  # Check and reduce K if necessary (if not all incremental claims are positive)
  K <- CheckK(S, K)

  if (!is.character(Weights_LR)) {
    if (!is.matrix(Weights_LR)) {
      stop("Weights_LR has to be a square matrix or one of the strings optimal or canonical.")
    } else if (ncol(Weights_LR) != n || nrow(Weights_LR) != n) {
      stop("Weights_LR has to be a square matrix or one of the strings optimal or canonical.")
    }
  } else {
    if (Weights_LR != "optimal" && Weights_LR != "canonical") {
      stop("Weights_LR has to be a square matrix or one of the strings optimal or canonical.")
    }
  }

  if (!is.null(eps_ext)) {
    if (!is.numeric(eps_ext)) {
      stop("eps_ext has to be NULL or numeric.")
    }
    if (eps_ext < 0) {
      stop("eps_ext must be non-negative.")
    }
  }
  if (!is.null(s_sq_ext)) {
    if (!is.numeric(s_sq_ext)) {
      stop("s_sq_ext has to be NULL or numeric.")
    }
    if (min(s_sq_ext) < 0) {
      stop("s_sq_ext must be positive.")
    }
    if (length(s_sq_ext) != n) {
      stop("s_sq_ext must have length n.")
    }
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
  rm(Temp)



  ReturnList <- list()
  class(ReturnList) <- "rl"

  eps_ext_is_0 <- F
  if (!is.null(eps_ext)) {
    if (eps_ext==0) {
      eps_ext_is_0 <- T
    }
  }


  if (!is.null(Seed)) {
    set.seed(seed = Seed)
  }
  #############################################################
  # Fall eps_ext = 0: GLR = LR und klassisches additives Modell
  #############################################################
  if (eps_ext_is_0) {
     if (is.matrix(Weights_LR)) {
      W <- Weights_LR
      ReturnList$Weights_LR_Type = "user"
    } else {
      W <- NULL
      ReturnList$Weights_LR_Type = "canonical"
    }

    print("Applying LR method")
    Results_LR <- LR_Method(S, v, W)

    W <- Results_LR$W_used
    ReturnList$Weights_LR <- W

    eps_hat_infty <- 0
    s_hat_sq_infty <- Results_LR$s_hat_sq
    if (!is.null(s_sq_ext)) {
      s_hat_sq_infty <- s_sq_ext
      ReturnList$s_sq_ext <- s_sq_ext
    }
    m_hat_infty <- Results_LR$m_hat
    v_hat_infty <- v

    ReturnList$Param_eps_hat_infty <- eps_hat_infty
    ReturnList$Param_v_hat_infty <- v
    ReturnList$Param_s_hat_sq_infty <- s_hat_sq_infty
    ReturnList$Param_m_hat_infty <- m_hat_infty

    MSE_LRM_epshatinfty <- rep(NA, n+1)
    RandError_LRM <- rep(NA, n+1)
    EstError_LRM <- rep(NA, n+1)

    if (Calculate_MSE_LR || Calculate_MSE_GLR) {
      print("Calculating MSEs")
      for (NrOfAlpha in 1:(n+1)) {
        MSE_LR_temp <- MSE_LR(S, W, Alpha[[NrOfAlpha]], 0 * c, Results_LR$m_hat, Results_LR$s_hat_sq, v)
        MSE_LRM_epshatinfty[NrOfAlpha] <- MSE_LR_temp$MSE
        EstError_LRM[NrOfAlpha] <- MSE_LR_temp$EstimationError
        RandError_LRM[NrOfAlpha] <- MSE_LR_temp$RandomError
      }
    }


    MSE_GLRM <- MSE_LRM_epshatinfty
    RandError_GLRM <- RandError_LRM
    EstError_GLRM <- EstError_LRM
    Bias_GLRM <- rep(0, n+1)

    ReturnList$Predictions_GLR <- i2c(Results_LR$Predictions, is.square = TRUE)
    ReturnList$se_GLR_fixed_weights_calc <- sqrt(MSE_GLRM)
    ReturnList$EstimationError_GLR_fixed_weights_calc <- sqrt(EstError_GLRM)
    ReturnList$RandomError_GLR_fixed_weights_calc <- sqrt(RandError_GLRM)
    ReturnList$Bias_GLR_fixed_weights_calc <- Bias_GLRM

    ReturnList$Predictions_LR <- i2c(Results_LR$Predictions, is.square = TRUE)
    ReturnList$se_LR_fixed_weights_calc <- sqrt(MSE_LRM_epshatinfty)
    ReturnList$EstimationError_LR_fixed_weights_calc <- sqrt(EstError_LRM)
    ReturnList$RandomError_LR_fixed_weights_calc <- sqrt(RandError_LRM)

    print("Simulating MSEs")

    SimulatedData <- SimulateClaims(v_hat_infty, m_hat_infty, s_hat_sq_infty, eps_hat_infty * c, NumberOfSimulations, K)
    Temp <- sqrt(Sim_MSE_bullet_LR_cpp(c(SimulatedData$Triangle), c(SimulatedData$volumes), n, NumberOfSimulations, W, v_hat_infty, m_hat_infty))
    ReturnList$se_LR_fixed_weights_sim <- Temp[1:(n+1)]
    ReturnList$RandomError_LR_fixed_weights_sim <- Temp[(n+2):(2*(n+1))]
    ReturnList$EstimationError_LR_fixed_weights_sim <- Temp[(2*(n+1)+1):(3*(n+1))]

    ReturnList$se_GLR_fixed_weights_sim <- ReturnList$se_LR_fixed_weights_sim
    ReturnList$RandomError_GLR_fixed_weights_sim <- ReturnList$RandomError_LR_fixed_weights_sim
    ReturnList$EstimationError_GLR_fixed_weights_sim <- ReturnList$EstimationError_LR_fixed_weights_sim

    ReturnList$se_LR_sim <- ReturnList$se_LR_fixed_weights_sim
    ReturnList$RandomError_LR_sim <- ReturnList$RandomError_LR_fixed_weights_sim
    ReturnList$EstimationError_LR_sim <- ReturnList$EstimationError_LR_fixed_weights_sim

    ReturnList$se_GLR_sim <- ReturnList$se_GLR_fixed_weights_sim
    ReturnList$RandomError_GLR_sim <- ReturnList$RandomError_LR_fixed_weights_sim
    ReturnList$EstimationError_GLR_sim <- ReturnList$EstimationError_LR_fixed_weights_sim

    ScalingFactor = 1
  } else {


    #############################################################
    # Fall eps_ext != 0
    #############################################################

    # Skalierung der Daten für mehr numerische Stabilität:
    ScalingFactor <- mean(v) * 10
    # ScalingFactor <- 1
    v <- v / ScalingFactor
    S <- S / ScalingFactor
    C <- C / ScalingFactor
    if (!is.null(s_sq_ext)) {
      ReturnList$s_sq_ext <- s_sq_ext
      s_sq_ext <- s_sq_ext / ScalingFactor
    }

    # Bestimme eps_start (falls NULL)
    if (is.null(eps_start)) {
      Indicator <- matrix(0, ncol = n-1, nrow = n)
      Indicator[row(Indicator)+col(Indicator) < n+1] <- 1
      CL_Factors <- apply(C[ , 2:n] * Indicator, 2, function (x) sum(x, na.rm = T)) / apply(C[ , 1:(n-1)] * Indicator, 2, function (x) sum(x, na.rm = T))
      Ultimates_CL <- c(1, cumprod(rev(CL_Factors))) * diag(C[,n:1])
      eps_start_sq <- sum(1 / c^2 * (Ultimates_CL[1:(n-1)] / Ultimates_CL[2:n] * v[2:n] / v[1:(n-1)] - 1)^2) / (n-1)
      eps_start <- sqrt(eps_start_sq)
      if (is.na(eps_start)) {
        eps_start <- 0.05
      }
    }

    eps_start <- max(eps_min, eps_start)
    eps_start <- min(eps_max, eps_start)

    if (UseRelativities_s_sq) {
      if (is.null(SetIndicesRelativities)) {
        # SetIndicesRelativities <- K
        SetIndicesRelativities <- 1:(n-1)
      }
      if (!is.null(Relativities_s_sq)) {
        if (is.numeric(Relativities_s_sq)) {
          if (length(Relativities_s_sq) != n) {
            stop("Relativities_s must be a numeric vector of length n or NULL")
          }
        } else {
          stop("Relativities_s must be a numeric vector of length n or NULL")
        }
      }
    } else {
      SetIndicesRelativities <- NULL
    }

    # Parameter estimation and generalized loss ratio method
    #Results_GLR <- Generalized_LR_Method(S, v, K = K, c = c, eps_ext = eps_ext, eps_start = eps_start, s_sq_ext = s_sq_ext, lower_bound_eps = eps_min, upper_bound_eps = eps_max, Iterations = Iterations, UseAllEstimators = UseAllEstimators, UseRcpp =  UseRcpp, Bias_Correction_m_hat_infty = Bias_Correction_m_hat_infty)
    Results_GLR <- Generalized_LR_Method(S, v, K = K, c = c, eps_ext = eps_ext, eps_start = eps_start, s_sq_ext = s_sq_ext, lower_bound_eps = eps_min, upper_bound_eps = eps_max, Iterations = Iterations, setS = SetIndicesRelativities, Relativities_s_sq = Relativities_s_sq, UseAllEstimators = UseAllEstimators, UseRcpp =  UseRcpp, adjust_v_hat_infty_for_small_eps_ext = adjust_v_hat_infty_for_small_eps_ext, Bias_Correction_m_hat_infty = Bias_Correction_m_hat_infty)
    g_hat_infty <- Results_GLR$g_hat_infty
    h_hat_infty <- Results_GLR$h_hat_infty
    eps_hat_infty <- Results_GLR$eps_hat_infty
    s_hat_sq_infty <- Results_GLR$s_hat_sq_infty
    m_hat_infty <- Results_GLR$m_hat_infty
    v_hat_infty <- Results_GLR$v_hat_infty
    ReturnList$Recursion_eps_hat <- matrix(Results_GLR$eps_hat, nrow=1)
    ReturnList$Recursion_v_hat <- Results_GLR$v_hat * ScalingFactor
    ReturnList$Recursion_m_hat <- Results_GLR$m_hat
    ReturnList$Recursion_s_hat_sq <- Results_GLR$s_hat_sq * ScalingFactor

    if (Calculate_MSE_GLR) {
      CovdotMatrix <- MSE_GLR(S, Alpha[[n+1]], K,  g_hat_infty, h_hat_infty, eps_hat_infty * c, m_hat_infty, s_hat_sq_infty, v_hat_infty, UseAllEstimators = UseAllEstimators,blnUseRcpp =  UseRcpp)$CovdotMatrix
    }

    MSE_LRM_epshatinfty <- rep(NA, n+1)
    MSE_GLRM <- rep(NA, n+1)
    RandError_LRM <- rep(NA, n+1)
    RandError_GLRM <- rep(NA, n+1)
    EstError_LRM <- rep(NA, n+1)
    EstError_GLRM <- rep(NA, n+1)
    Bias_GLRM <- rep(NA, n+1)


    if (is.character(Weights_LR)) {
      if (Weights_LR == "optimal") {
        W <- OptimalWeights_LR_Method(v_hat_infty, m_hat_infty, s_hat_sq_infty, eps_hat_infty * c)
        #W <- OptimalWeights_LR_Method(v, m_hat_infty, s_hat_sq_infty, eps_hat_infty * c)
        ReturnList$Weights_LR_Type <- "optimal"
      } else if (Weights_LR == "canonical") {
        W <- matrix(v, ncol = 1) %*% matrix(1, nrow=1, ncol = n)
        W[col(W) + row(W)> n+1] <- 0
        for (k in 1:n) {
          W[, k] <- W[, k] / sum(W[, k])
        }
        ReturnList$Weights_LR_Type <- "canonical"
      }
    } else {
      W <- Weights_LR
      W[col(W) + row(W)> n+1] <- 0
      for (k in 1:n) {
        W[, k] <- W[, k] / sum(W[, k])
      }
      ReturnList$Weights_LR_Type <- "user"
    }

    ReturnList$Weights_LR <- W

    Results_LR <- LR_Method(S, v, W)

    if (Calculate_MSE_LR || Calculate_MSE_GLR) {
      print("Calculating MSEs:")
      pb <- txtProgressBar(min = 0, max = 1, style = 3)
      progress <- 0
      setTxtProgressBar(pb, progress)

      for (NrOfAlpha in 1:(n+1)) {
        if (Calculate_MSE_LR) {
          MSE_LR_temp <- MSE_LR(S, W, Alpha[[NrOfAlpha]], eps_hat_infty * c, m_hat_infty, s_hat_sq_infty, v_hat_infty)
          MSE_LRM_epshatinfty[NrOfAlpha] <- MSE_LR_temp$MSE
          EstError_LRM[NrOfAlpha] <- MSE_LR_temp$EstimationError
          RandError_LRM[NrOfAlpha] <- MSE_LR_temp$RandomError
        }
        if (Calculate_MSE_GLR) {
          MSE_GLR_temp <- MSE_GLR(S, Alpha[[NrOfAlpha]], K,  g_hat_infty, h_hat_infty, eps_hat_infty * c, m_hat_infty, s_hat_sq_infty, v_hat_infty, CovdotMatrix = CovdotMatrix, UseAllEstimators = UseAllEstimators, blnUseRcpp =  UseRcpp)
          MSE_GLRM[NrOfAlpha] <- MSE_GLR_temp$MSE
          EstError_GLRM[NrOfAlpha] <- MSE_GLR_temp$PredictionVariance + MSE_GLR_temp$Bias^2
          RandError_GLRM[NrOfAlpha] <- MSE_GLR_temp$ProcessVariance
          Bias_GLRM[NrOfAlpha] <- MSE_GLR_temp$Bias
        }
        progress <- NrOfAlpha / (n+1)
        setTxtProgressBar(pb, progress)
      }
      close(pb)
    }

    print("Simulating MSEs with fixed weights")

    SimulatedData <- SimulateClaims(v_hat_infty, m_hat_infty, s_hat_sq_infty, eps_hat_infty * c, NumberOfSimulations, K)

    # Matrix J[[i]] contains the set $J_i$
    J <- CreateJ(S, K, UseAllEstimators)
    # J[[i]] = used list of indices for the estimators of v_i^{-1}; for \hatr_{i,(j,k)}
    #         - J[[i]][,1] = i, i.e. estimator for v_i^{-1}
    #         - J[[i]][,2] = j
    #         - J[[i]][,3] = k
    J_all <- J[[1]]
    for (i in 2:n) {
      J_all <- rbind(J_all,J[[i]])
    }
    J_s <- c(J_all)
    Number_of_Est_s <- numeric(n)
    g_s <- numeric()
    h_s <- numeric()

    for (i in 1:n) {
      Number_of_Est_s[i] <- nrow(J[[i]])
      g_s <- c(g_s,g_hat_infty[[i]])
      h_s <- c(h_s,h_hat_infty[[i]])
    }

    # ReturnList$se_GLR_fixed_weights_sim <- Sim_MSE_bullet_GLR_cpp(c(SimulatedData$Triangle), rep(v, NumberOfSimulations), n, NumberOfSimulations, K, J_s, Number_of_Est_s, h_s, g_s, m_hat_infty, s_hat_sq_infty)
    Temp <- ScalingFactor * sqrt(Sim_MSE_bullet_GLR_cpp(c(SimulatedData$Triangle), c(SimulatedData$volumes), n, NumberOfSimulations, K, J_s, Number_of_Est_s, h_s, g_s, v_hat_infty, m_hat_infty, s_hat_sq_infty))
    ReturnList$se_GLR_fixed_weights_sim <- Temp[1:(n+1)]
    ReturnList$RandomError_GLR_fixed_weights_sim <- Temp[(n+2):(2*(n+1))]
    ReturnList$EstimationError_GLR_fixed_weights_sim <- Temp[(2*(n+1)+1):(3*(n+1))]

    Temp <- ScalingFactor * sqrt(Sim_MSE_bullet_LR_cpp(c(SimulatedData$Triangle), c(SimulatedData$volumes), n, NumberOfSimulations, W, v_hat_infty, m_hat_infty))
    ReturnList$se_LR_fixed_weights_sim <- Temp[1:(n+1)]
    ReturnList$RandomError_LR_fixed_weights_sim <- Temp[(n+2):(2*(n+1))]
    ReturnList$EstimationError_LR_fixed_weights_sim <- Temp[(2*(n+1)+1):(3*(n+1))]

    if (is.character(Weights_LR)) {
      if (Weights_LR == "canonical") {
        Temp  <- ScalingFactor * sqrt(Sim_MSE_LR_canonical_cpp(c(SimulatedData$Triangle), c(SimulatedData$volumes), n, NumberOfSimulations, v_hat_infty, m_hat_infty))
        ReturnList$se_LR_sim  <- Temp[1:(n+1)]
        ReturnList$RandomError_LR_sim <- Temp[(n+2):(2*(n+1))]
        ReturnList$EstimationError_LR_sim <- Temp[(2*(n+1)+1):(3*(n+1))]
      }
    }

    Temp <- ScalingFactor * sqrt(Sim_MSE_CL_cpp(c(SimulatedData$Triangle), n, NumberOfSimulations, v_hat_infty, m_hat_infty))
    ReturnList$se_CL_sim <- Temp[1:(n+1)]
    ReturnList$RandomError_CL_sim <- Temp[(n+2):(2*(n+1))]
    ReturnList$EstimationError_CL_sim <- Temp[(2*(n+1)+1):(3*(n+1))]


    if (RunAllSimulations) {
      Predictions_GLR_sim <- array(NA, dim = c(n,n,NumberOfSimulations))
      Predictions_LR_sim <- array(NA, dim = c(n,n,NumberOfSimulations))

      print(paste0("Simulating MSEs with unfixed weights"))
      pb <- txtProgressBar(min = 0, max = 1, style = 3)
      progress <- 0
      setTxtProgressBar(pb, progress)

      for (SimNr in 1:NumberOfSimulations) {
        if (Use_eps_ext_inSimMSE && Use_s_sq_ext_inSimMSE) {
          GLR_sim <- Generalized_LR_Method(SimulatedData$Triangle[,,SimNr],SimulatedData$volumes[,SimNr], eps_start = eps_hat_infty, eps_ext = eps_ext, lower_bound_eps = eps_min, upper_bound_eps = eps_max, s_sq_ext = s_sq_ext, K = K, c = c, UseAllEstimators = UseAllEstimators, Iterations = IterationsInSimMSE, setS = SetIndicesRelativities, Relativities_s_sq = Relativities_s_sq, ShowProgressBar = F, UseRcpp = UseRcpp, adjust_v_hat_infty_for_small_eps_ext = adjust_v_hat_infty_for_small_eps_ext, Bias_Correction_m_hat_infty = Bias_Correction_m_hat_infty)
        } else if (Use_eps_ext_inSimMSE) {
          GLR_sim <- Generalized_LR_Method(SimulatedData$Triangle[,,SimNr],SimulatedData$volumes[,SimNr], eps_start = eps_hat_infty, eps_ext = eps_ext, lower_bound_eps = eps_min, upper_bound_eps = eps_max, K = K, c = c, UseAllEstimators = UseAllEstimators, Iterations = IterationsInSimMSE, setS = SetIndicesRelativities, Relativities_s_sq = Relativities_s_sq, ShowProgressBar = F, UseRcpp = UseRcpp, adjust_v_hat_infty_for_small_eps_ext = adjust_v_hat_infty_for_small_eps_ext, Bias_Correction_m_hat_infty = Bias_Correction_m_hat_infty)
        } else if (Use_s_sq_ext_inSimMSE) {
          GLR_sim <- Generalized_LR_Method(SimulatedData$Triangle[,,SimNr],SimulatedData$volumes[,SimNr], eps_start = eps_hat_infty, lower_bound_eps = eps_min, upper_bound_eps = eps_max, s_sq_ext = s_sq_ext, K = K, c = c, UseAllEstimators = UseAllEstimators, Iterations = IterationsInSimMSE, setS = SetIndicesRelativities, Relativities_s_sq = Relativities_s_sq, ShowProgressBar = F, UseRcpp = UseRcpp, adjust_v_hat_infty_for_small_eps_ext = adjust_v_hat_infty_for_small_eps_ext, Bias_Correction_m_hat_infty = Bias_Correction_m_hat_infty)
        } else {
          GLR_sim <- Generalized_LR_Method(SimulatedData$Triangle[,,SimNr],SimulatedData$volumes[,SimNr], eps_start = eps_hat_infty, lower_bound_eps = eps_min, upper_bound_eps = eps_max, K = K, c = c, UseAllEstimators = UseAllEstimators, Iterations = IterationsInSimMSE, setS = SetIndicesRelativities, Relativities_s_sq = Relativities_s_sq, ShowProgressBar = F, UseRcpp = UseRcpp, adjust_v_hat_infty_for_small_eps_ext = adjust_v_hat_infty_for_small_eps_ext, Bias_Correction_m_hat_infty = Bias_Correction_m_hat_infty)
        }
        Predictions_GLR_sim[,,SimNr] <- GLR_sim$Predictions
        W_opt_sim <- OptimalWeights_LR_Method(GLR_sim$v_hat_infty, GLR_sim$m_hat_infty, GLR_sim$s_hat_sq_infty, GLR_sim$eps_hat_infty * c)
        LR_sim <- LR_Method(SimulatedData$Triangle[,,SimNr], SimulatedData$volumes[,SimNr], W_opt_sim)
        Predictions_LR_sim[,,SimNr] <- LR_sim$Predictions
        progress <- SimNr / NumberOfSimulations
        setTxtProgressBar(pb, progress)

      }

      close(pb)

      FutureIncrements_sim <- SimulatedData$Triangle

      for (i in 1:n) {
        Predictions_GLR_sim[i,1:(n-i+1),] <- 0
        Predictions_LR_sim[i,1:(n-i+1),] <- 0
        FutureIncrements_sim[i,1:(n-i+1),] <- 0
      }

      Reserves_GLR_sim <- apply(Predictions_GLR_sim, c(1,3), 'sum')
      Reserves_LR_sim <- apply(Predictions_LR_sim, c(1,3), 'sum')
      Reserves_Sim <- apply(FutureIncrements_sim, c(1,3), 'sum')
      Reserves_GLR_sim <- rbind(Reserves_GLR_sim, apply(Reserves_GLR_sim, 2, 'sum'))
      Reserves_LR_sim <- rbind(Reserves_LR_sim, apply(Reserves_LR_sim, 2, 'sum'))
      Reserves_Sim <- rbind(Reserves_Sim, apply(Reserves_Sim, 2, 'sum'))

      # Ultimates_GLR_sim <- apply(Predictions_GLR_sim, c(1,3), 'sum')
      # Ultimates_LR_sim <- apply(Predictions_LR_sim, c(1,3), 'sum')
      # Ultimates_Sim <- apply(SimulatedData$Triangle, c(1,3), 'sum')
      # Ultimates_GLR_sim <- rbind(Ultimates_GLR_sim, apply(Ultimates_GLR_sim, 2, 'sum'))
      # Ultimates_LR_sim <- rbind(Ultimates_LR_sim, apply(Ultimates_LR_sim, 2, 'sum'))
      # Ultimates_Sim <- rbind(Ultimates_Sim, apply(Ultimates_Sim, 2, 'sum'))

      ExpectedReserves <- v_hat_infty %*% t(m_hat_infty)
      ExpectedReserves[col(ExpectedReserves)+row(ExpectedReserves) < n+2] <- 0
      ExpectedReserves <- apply(ExpectedReserves, 1, 'sum')
      ExpectedReserves <- c(ExpectedReserves, sum(ExpectedReserves))

      qa_GLR <- (Reserves_Sim - Reserves_GLR_sim)^2
      qa_re_GLR <- (Reserves_Sim - ExpectedReserves)^2
      qa_ee_GLR <- (Reserves_GLR_sim - ExpectedReserves)^2

      qa_LR <- (Reserves_Sim - Reserves_LR_sim)^2
      qa_re_LR <- (Reserves_Sim - ExpectedReserves)^2
      qa_ee_LR <- (Reserves_LR_sim - ExpectedReserves)^2

      ReturnList$se_GLR_sim <- sqrt(apply(qa_GLR,1,'sum')/NumberOfSimulations) * ScalingFactor
      ReturnList$RandomError_GLR_sim <- sqrt(apply(qa_re_GLR,1,'sum')/NumberOfSimulations) * ScalingFactor
      ReturnList$EstimationError_GLR_sim <- sqrt(apply(qa_ee_GLR,1,'sum')/NumberOfSimulations) * ScalingFactor

      if (Weights_LR == 'optimal') {
        ReturnList$se_LR_sim <- sqrt(apply(qa_LR,1,'sum')/NumberOfSimulations) * ScalingFactor
        ReturnList$RandomError_LR_sim <- sqrt(apply(qa_re_LR,1,'sum')/NumberOfSimulations) * ScalingFactor
        ReturnList$EstimationError_LR_sim <- sqrt(apply(qa_ee_LR,1,'sum')/NumberOfSimulations) * ScalingFactor
      }

      # qa_GLR <- (Ultimates_Sim - Ultimates_GLR_sim)^2
      # qa_LR <- (Ultimates_Sim - Ultimates_LR_sim)^2
      # ReturnList$se_GLR_sim <- sqrt(apply(qa_GLR,1,'sum')/NumberOfSimulations) * ScalingFactor
      # if (Weights_LR == 'optimal') {
      #   ReturnList$se_LR_sim <- sqrt(apply(qa_LR,1,'sum')/NumberOfSimulations) * ScalingFactor
      # }
      #


    } else {
      ReturnList$se_GLR_sim <- rep(" not simulated", n+1)
      if (is.character(Weights_LR)) {
        if (Weights_LR == "optimal") {
          ReturnList$se_LR_sim <- rep(" not simulated", n+1)
        }
      }
    }

    if (!is.character(Weights_LR)) {
      ReturnList$se_LR_sim <- ReturnList$se_LR_fixed_weights_sim
    }


    ReturnList$Param_eps_hat_infty <- eps_hat_infty
    ReturnList$Param_v_hat_infty <- v_hat_infty * ScalingFactor
    ReturnList$Param_s_hat_sq_infty <- s_hat_sq_infty * ScalingFactor
    ReturnList$Param_m_hat_infty <- m_hat_infty

    ReturnList$Predictions_GLR <- i2c(Results_GLR$Predictions, is.square = TRUE) * ScalingFactor
    ReturnList$se_GLR_fixed_weights_calc <- sqrt(MSE_GLRM) * ScalingFactor
    ReturnList$EstimationError_GLR_fixed_weights_calc <- sqrt(EstError_GLRM) * ScalingFactor
    ReturnList$RandomError_GLR_fixed_weights_calc <- sqrt(RandError_GLRM) * ScalingFactor
    ReturnList$Bias_GLR_fixed_weights_calc <- Bias_GLRM * ScalingFactor

    ReturnList$Predictions_LR <- i2c(Results_LR$Predictions, is.square = TRUE) * ScalingFactor
    ReturnList$se_LR_fixed_weights_calc <- sqrt(MSE_LRM_epshatinfty) * ScalingFactor
    ReturnList$EstimationError_LR_fixed_weights_calc <- sqrt(EstError_LRM) * ScalingFactor
    ReturnList$RandomError_LR_fixed_weights_calc <- sqrt(RandError_LRM) * ScalingFactor

  } # Ende Fall eps_ext != 0


  ReturnList$UseRelativities_s_sq = UseRelativities_s_sq
  ReturnList$Relativities_s_sq = Relativities_s_sq
  ReturnList$SetIndicesRelativities = SetIndicesRelativities

  Summary_GLRM <- data.frame(v_hat_infty = c(v_hat_infty * ScalingFactor, sum(v_hat_infty) * ScalingFactor))
  Summary_GLRM$epsilon_hat_infty <- c(eps_hat_infty * c, NA, NA)
  temp <- rep(0,n+1)
  temp[1:n] <- ReturnList$Predictions_GLR[,n]
  temp[n+1] <- sum(temp[1:n])
  Summary_GLRM$Ultimates <- temp
  for (i in 1:n) {
    temp[i] <- Summary_GLRM$Ultimates[i] - ReturnList$Predictions_GLR[i,n-i+1]
  }
  temp[n+1] <- sum(temp[1:n])
  Summary_GLRM$Reserves <- temp
  ReturnList$Reserves_GLR <- temp
  if (is.na(ReturnList$se_GLR_fixed_weights_calc[1])) {ReturnList$se_GLR_fixed_weights_calc <- rep(" not calculated", n+1)}
  Summary_GLRM$se_fixed_weights_calc <- ReturnList$se_GLR_fixed_weights_calc
  # Summary_GLRM$RandomError <- ReturnList$RandomError_GLR
  # Summary_GLRM$EstimationError <- ReturnList$EstimationError_GLR
  # Summary_GLRM$Bias <- ReturnList$Bias_GLR
  Summary_GLRM$se_fixed_weights_sim <- ReturnList$se_GLR_fixed_weights_sim
  Summary_GLRM$se_sim <- ReturnList$se_GLR_sim
  rownames(Summary_GLRM)[n+1] <- "Total"

  ReturnList$Summary_GLR <- Summary_GLRM

  v <- v * ScalingFactor
  Summary_LRM <- data.frame(v_hat = c(v,sum(v)) , epsilon_hat_infty = c(eps_hat_infty * c, NA, NA))
  temp <- rep(0,n+1)
  temp[1:n] <- ReturnList$Predictions_LR[,n]
  temp[n+1] <- sum(temp[1:n])
  Summary_LRM$Ultimates <- temp
  for (i in 1:n) {
    temp[i] <- Summary_LRM$Ultimates[i] - ReturnList$Predictions_LR[i,n-i+1]
  }
  temp[n+1] <- sum(temp[1:n])
  Summary_LRM$Reserves <- temp
  ReturnList$Reserves_LR <- temp
  if (is.na(ReturnList$se_LR_fixed_weights_calc[1])) {ReturnList$se_LR_fixed_weights_calc <- rep(" not calculated", n+1)}
  Summary_LRM$se_fixed_weights_calc <- ReturnList$se_LR_fixed_weights_calc
  Summary_LRM$se_fixed_weights_sim <- ReturnList$se_LR_fixed_weights_sim
  Summary_LRM$se_sim <- ReturnList$se_LR_sim

  # Summary_LRM$RandomError <- ReturnList$RandomError_LR
  # Summary_LRM$EstimationError <- ReturnList$EstimationError_LR

  rownames(Summary_LRM)[n+1] <- "Total"

  ReturnList$Summary_LR <- Summary_LRM


  ReturnList$K <- CheckK(S, K)
  ReturnList$c <- c
  ReturnList$eps_ext <- eps_ext
  ReturnList$NumberOfSimulations <- NumberOfSimulations
  ReturnList$Iterations <- Iterations
  ReturnList$eps_start <- eps_start

  return(ReturnList)


}


#' @export
print.rl <- function(ReturnList) {
  cat("\nSummary for Generalized Loss Ratio Method:\n\n")
  print(ReturnList$Summary_GLR)
  cat("\nUsed K: ", ReturnList$K,"\n")
  cat("Used c: ", ReturnList$c,"\n")
  if (!is.null(ReturnList$eps_ext)) {
    cat("External estimator for epsilon: ", ReturnList$eps_ext * 100,"%\n")
  } else {
    cat("Used eps_start: ", ReturnList$eps_start * 100,"%\n")
  }
  if (!is.null(ReturnList$s_sq_ext)) {
    cat("External estimators for s_k^2: ", ReturnList$s_sq_ext,"\n")
  } else {
    if (ReturnList$UseRelativities_s_sq) {
      if (is.null(ReturnList$Relativities_s_sq)) {
        cat("Relativities for s_k^2 used:\n      Relativities: from loss ratio method\n      Set of indices: ", ReturnList$SetIndicesRelativities, "\n")
      } else {
        cat("Relativities for s_k^2 used:\n      Relativities: ", ReturnList$Relativities_s_sq, "\n      Set of indices: ", ReturnList$SetIndicesRelativities, "\n")
      }

    }
  }
  cat("Number of Iterations: ", ReturnList$Iterations,"\n")
  cat("Number of Simulations: ", ReturnList$NumberOfSimulations,"\n")

  if (ReturnList$Weights_LR_Type == "optimal") {
    cat("\n\nSummary for Loss Ratio Method in the extended additive model (approximately optimal weights):\n\n")
  } else if (ReturnList$Weights_LR_Type == "canonical") {
    cat("\n\nSummary for Loss Ratio Method in the extended additive model (canonical weights):\n\n")
  } else {
    cat("\n\nSummary for Loss Ratio Method in the extended additive model (user defined weights):\n\n")
  }
  print(ReturnList$Summary_LR)
  cat("\nNumber of Simulations: ", ReturnList$NumberOfSimulations,"\n")

  # cat("\nTo obtain all information use unclass(...).\n")
}
