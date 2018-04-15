#' This function provides the Chain Ladder Method
#' @param triangle Numeric matrix (n x n) containing the claims triangle
#' @param is.incremental  Logical. Indicates whether the claims triangle is cumulative (default) or incremental.

#' @return The function returns a list with various results:
#'
#' \itemize{
#' \item Predictions: Predictions of the chain ladder method (cumulative)
#' \item Reserves: Reserves of the chain ladder method
#' \item MSE: MSE[i] = mean squared error of prediction for accident year i; MSE[n+1] = mean squared error of prediction for the total reserve
#' }

#' @examples
#' CL(glrm_example1_C)
#' CL(glrm_example1_S, is.incremental = FALSE)

#' @export
CL <- function(triangle, is.incremental = FALSE) {
  ###########################################################################################################
  # This function provides the chain ladder method
  ###########################################################################################################
  if (!is.matrix(triangle)) {
    stop("Triangle must be a matrix.")
  } else if (nrow(triangle) != ncol(triangle)) {
    stop("Triangle must be a square matrix.")
  } else if (!is.numeric(triangle)) {
    stop("Triangle must be a numeric matrix.")
  }

  n <- nrow(triangle)

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

  f <- numeric(n-1)
  sigma_sq <- numeric(n-1)
  for (k in 1:(n-1)) {
    f[k] <- sum(C[1:(n-k), k+1]) / sum(C[1:(n-k), k])
  }
  for (k in (1:(n-2))) {
    sigma_sq[k] <- sum(C[1:(n-k), k] * (C[1:(n-k), k+1] / C[1:(n-k), k] - f[k])^2) / (n-k-1)
  }
  sigma_sq[n-1] <- min(sigma_sq[n-2]^2 / sigma_sq[n-3], sigma_sq[n-3])
  Predictions <- C
  for (k in 2:n) {
    for (i in (n-k+2):n) {
      Predictions[i,k] <- Predictions[i, k-1] * f[k-1]
    }
  }

  MSE <- numeric(n+1)

  for (i in 2:n) {
    temp <- 0
    for (k in (n+1-i):(n-1)) {
      temp <- temp + sigma_sq[k] / f[k]^2 * (1 / Predictions[i,k] + 1 / sum(Predictions[1:(n-k),k]))
    }
    MSE[i] <- Predictions[i,n]^2 * temp
  }

  # MSE Total Reserve
  MSE[n+1] <- sum(MSE[1:n])
  for (i in 2:(n-1)) {
    temp <- 0
    for (k in (n+1-i):(n-1)) {
      temp <- temp + 2 * sigma_sq[k] / f[k]^2 / sum(Predictions[1:(n-k),k])
    }
    MSE[n+1] <- MSE[n+1] + Predictions[i,n] * sum(Predictions[(i+1):n,n]) * temp
  }

  Result <- list()
  rownames(Predictions) <- paste0("AY", 1:n)
  colnames(Predictions) <- paste0("DY", 1:n)
  Result$Predictions <- Predictions
  Reserves <- Predictions[, n] - diag(Predictions[, n:1])
  Reserves <- c(Reserves, sum(Reserves))
  names(Reserves) <- paste0(c(rep("AY", n), "total"), c(1:n, ""))
  Result$Reserves <- Reserves
  names(MSE) <- paste0(c(rep("AY", n), "total"), c(1:n, ""))
  Result$MSE <- MSE

  return(Result)
}
