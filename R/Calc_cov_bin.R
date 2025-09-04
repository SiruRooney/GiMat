#' Joint probability between two indicators
#'
#' @param a a vector with an indicator variable.
#' @param b a vector with an indicator variable.
#' @param log logical. If \code{TRUE} returns the log probability. 
#'
#' @return joint probability
#' @export joint_p
#'
#' @examples
#' X <- rbinom(1000, size = 1, p = .25)
#' Y <- rbinom(1000, size = 1, p = .75)
#' joint_p(X, Y)
joint_p <- function(a, b, log = FALSE){
  if(length(a) != length(b)) stop("Vectors not of the same size")
  ans <- log(mean(a*b))
  if(!log) ans <- exp(ans)
  return(ans)
}


#' Covariance between binary variables.
#'
#' @param pij joint probability.
#' @param pi marginal probability of the first variable.
#' @param pj marginal probability of the second variable.
#'
#' @return covariance.
#' @export bin_cov
#'
#' @examples
#' bin_cov(pij = 0.5, pi = .2, pj = .2)
#' bin_cov(pij = 0.0, pi = .2, pj = .2)
bin_cov <- function(pij, pi, pj){
  return( pij - exp(log(pi) + log(pj)) )
}


#' Correlation between binary variables#'
#' @param pij joint probability.
#' @param pi marginal probability of the first variable.
#' @param pj marginal probability of the second variable.
#'
#' @return correlation coefficient.
#' @export bin_corr
#'
#' @examples
#' bin_corr(pij = 0.1, pi = 0.1, pj = .2) 
#' bin_corr(pij = 0.0, pi = .2, pj = .2)
#' # bin_corr(pij = 0.5, pi = .2, pj = .2) ## Error
bin_corr <- function(pij, pi, pj){
  num <- pij - exp(log(pi) + log(pj))
  sgn <- sign(num)
  lnum <- log(abs(num))
  ldenom <- 0.5 * ( log(pi) + log1p(-pi) + log(pj) + log1p(-pj) )
  ans <- sgn * exp(lnum - ldenom)
  if(ans > 1)  stop("Computed correlation is larger than 1, check your probabilities")
  return(ans)
}


#' Covariance matrix for binary variables
#'
#' @param x a matrix or data.frame containing only binary variables.
#' @param ncores (optional) number of cores to be used.
#'
#' @return covariance matrix.
#' @importFrom stats var
#' @export binary_cov_matrix
binary_cov_matrix <- function(x, ncores = 1){
  Var1 <- Var2 <- NULL
  Vs <- apply(x, 2, var)
  if(any(Vs ==0)) {
    warning("There are zero-variance variables!")
    return(NaN)
  }
  N <- ncol(x)
  grid <- subset(expand.grid(1:N, 1:N), Var1 < Var2)
  ps <- as.numeric(colMeans(x))
  covs <- unlist(
    parallel::mclapply(1:nrow(grid), function(row) {
      i <- grid[row, 1]
      j <- grid[row, 2]
      ans <- bin_cov(pij = joint_p(x[, i], x[, j]), pi = ps[i], pj = ps[j])
      return(ans)
    }, mc.cores = ncores)
  )
  all.comps <- cbind(grid, covs)
  cov.mat <- matrix(NA, ncol = N, nrow = N)
  diag(cov.mat) <- Vs
  for(k in 1:nrow(all.comps)){
    i <- all.comps[k, 1]
    j <- all.comps[k, 2]
    cov.mat[i, j] <- cov.mat[j, i] <- all.comps[k, 3]
  }
  return(cov.mat)
}









