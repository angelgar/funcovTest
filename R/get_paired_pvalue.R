#' Generate p-value for the paired test
#'
#' @description
#' Internal function
#'
#' @param X A matrix of functions
#' @param Y A matrix of functions
#' @param n_perm Number of permutations used to calculate p-value for paired test
#'
#' @keywords internal
#'
#'
get_paired_pvalue <- function(X, Y, no_perm = 1000) {

  # Get observed statistic
  observed_statistic <- get_paired_statistic(X, Y)

  # Get vector of permuted statistic
  stat_perm_vals <- replicate(no_perm, {
    perm_index <- sample(c(T, F), size = nrow(X), replace = T)
    X_perm <- rbind(X[perm_index, ], Y[!perm_index, ])
    Y_perm <- rbind(Y[perm_index, ], X[!perm_index, ])
    get_paired_statistic(X_perm, Y_perm)
  })

  # calculate p-value
  p_val <- sum(stat_perm_vals > observed_statistic) / no_perm

  return(p_val)

}
