#' Estimate the test statistic for the paired data
#'
#' @description
#' Internal function
#'
#' @param X A matrix of functions
#' @param Y A matrix of functions
#'
#' @keywords internal
#'
#'
#'
get_paired_statistic <- function(X, Y) {

  # Get dimensions of X and Y
  p <- ncol(X)
  n1 <- nrow(X)
  n2 <- nrow(Y)

  # Calculate sample covariances, scaled by (n-1)/n
  Sx <- stats::cov(X) * (n1 - 1) / n1
  Sy <- stats::cov(Y) * (n2 - 1) / n2

  # Calculate centered data matrices
  xa <- t(t(X) - colMeans(X))
  ya <- t(t(Y) - colMeans(Y))

  # Calculate terms needed for the denominators of the test statistic (scale)
  vx <- t(xa^2) %*% (xa^2) / n1 - 2/n1 * (t(xa) %*% xa) * Sx + Sx^2
  vy <- t(ya^2) %*% (ya^2) / n2 - 2/n2 * (t(ya) %*% ya) * Sy + Sy^2

  # Calculate psi_xy matrix
  psi_xy <- matrix(0, nrow = p, ncol = p)
  for (j in 1:n1) {
    psi_xy <- psi_xy + ((xa[j, ] %*% t(xa[j, ])) * (ya[j, ] %*% t(ya[j, ])))
  }
  psi_xy <- ((psi_xy / n1) - (Sx * Sy))

  # Calculate the denominator for the test statistic
  denom_new <- vx/n1 + vy/n1 - (2 * psi_xy/n1)

  # Calculate the paired test statistic
  CLX_new <- max((Sx - Sy)^2 / (denom_new))

  return(CLX_new)
}
