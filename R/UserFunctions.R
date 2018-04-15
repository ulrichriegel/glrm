#' This function calculates the cumulative claims triangle from an incremental claims triangle
#' @param triangle Numeric matrix (n x n). Incremental claims triangle or claims square
#' @param is.square Logical. FALSE if only a triangle (i.e. entries with indices i+j <= n+1) is given. TRUE if a claims square is given
#' @export

i2c <- function (triangle, is.square = FALSE) {
  if (!is.matrix(triangle)) {
    stop("Triangle must be a matrix.")
  }
  if (nrow(triangle) != ncol(triangle)) {
    stop("Triangle must be a square matrix.")
  }
  n <- nrow(triangle)
  cum <- t(apply(triangle, 1, cumsum))
  dimnames(cum) <- dimnames(triangle)
  if (!is.square) {
    cum[row(cum)+col(cum)>n+1] <- NA
  }
  return(cum)
}

#' This function calculates the incremental claims triangle from a cumulative claims triangle
#' @param triangle Numeric matrix (n x n). Cumulative claims triangle or claims square
#' @param is.square Logical. FALSE if only a triangle (i.e. entries with indices i+j <= n+1) is given. TRUE if a claims square is given
#' @export

c2i <- function (triangle, is.square = FALSE) {
  if (!is.matrix(triangle)) {
    stop("Triangle must be a matrix.")
  }
  if (nrow(triangle) != ncol(triangle)) {
    stop("Triangle must be a square matrix.")
  }
  n <- nrow(triangle)
  incr <- cbind(triangle[,1], t(apply(triangle, 1, diff)))
  dimnames(incr) <- dimnames(triangle)
  if (!is.square) {
    incr[row(incr)+col(incr)>n+1] <- NA
  }
  return(incr)
}

