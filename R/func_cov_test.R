#' Two sample test for eigendecompositions of functional data
#'
#' @description
#' Function to compare the eigendecomposition of two sets of functional data.
#'
#' @param X A matrix of functions
#' @param Y A matrix of functions
#' @param paired A logical indidating wether you want a paired test
#' @param n_perm Number of permutations used to calculate p-value for paired test
#' @param pve proportion varaince explained; a value between 0-1
#'
#' @export
#'
#'
func_cov_test <- function(X, Y, paired = FALSE, n_perm = 1000, pve = 0.99) {

  if (!paired) {
    p.val <- get_indep_pval(X, Y)
    return(p.val)

  } else {
    # Combine the two datasets and convert to a matrix
    all_data <- cbind(X, Y)
    all_data <- as.matrix(all_data)

    # Perform functional principal component analysis
    refund_output <- refund::fpca.face(Y = t(all_data), pve = pve)
    scores_all <- refund_output$scores

    N.subj <- ncol(X)

    # Separate the scores into two datasets
    scores_data1 <- scores_all[1:N.subj,]
    scores_data2 <- scores_all[-(1:N.subj),]

    p.val <- get_paired_pvalue(scores_data1, scores_data2, n_perm)
    return(p.val)
  }

}
