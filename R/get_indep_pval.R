#' Generate p-value for the independent test
#'
#' @description
#' Internal function
#'
#' @param X A matrix of functions
#' @param Y A matrix of functions
#' @param pve proportion varaince explained; a value between 0-1
#'
#' @keywords internal
#'
#'
#'
get_indep_pval <- function(X, Y, pve) {

  # Combine the two datasets and convert to a matrix
  all_data <- cbind(X, Y)
  all_data <- as.matrix(all_data)

  # Perform functional principal component analysis
  refund_output <- refund::fpca.face(Y = t(all_data), pve = pve)
  scores_all <- refund_output$scores

  # Extract the number of subjects in the first dataset
  N.subj <- dim(X)[2]

  # Separate the scores into two datasets
  scores_data1 <- scores_all[1:N.subj, ]
  scores_data2 <- scores_all[-(1:N.subj), ]

  # Perform hypothesis test and extract the p-value
  test_output <- HDtest::testCov(scores_data1, scores_data2, method = "CLX")
  p.val.clx <- test_output$p.value

  return(p.val.clx)
}
