#' @title Estimate TF activity score in a sample
#'
#' @description \code{InferTFact} is an auxiliary function for function \code{sciraRegAct}.
#'
#' @param exp A numeric vector of gene expression levels for all TF target genes in a sample.
#' @param regnet A matrix, the network estimated by function \code{sciraInfNet}, with +1 referring to positive regulation, -1 negative regulation, and 0 no regulation.
#'
#' @return A vector storing the activity scores of all TFs in a specific sample.
#'
#' @details \code{InferTFact} regresses the expression profiles of TF target genes against the binding profile of this TF on these genes. The output t-statistics are taken as the TF activity scores.
#'
#' @importFrom stats lm
#'
#' @examples
#' # TF regulatory network with 100 genes and 5 TFs.
#' regnet.m <- matrix(sample(c(-1, 0, 1), 500, replace = TRUE), nrow = 100)
#' # gene expression vector (for one sample)
#' exp.v <- rnorm(100)
#' # TF activity score
#' TFact.v <- SCIRA:::InferTFact(exp.v, regnet.m)

InferTFact <- function(exp, regnet) {
  act.v <- apply(regnet, 2, function(tmp.v) {
    lm.o <- lm(exp ~ tmp.v)
    act <- summary(lm.o)$coeff[2, 3]
    return(act)
  })

  return(act.v)
}

#' @title Estimate TF activity score in different samples
#'
#' @description \code{InferTFactPRL } is an auxiliary function for function \code{sciraRegAct}.
#'
#' @param idx A numeric vector indicating the samples from which you want to estimate the TF activity score.
#' @param tmp A numeric matrix of gene expression levels for all TF target genes, with rows referring to genes, columns to samples.
#' @param regnet The network estimated from function \code{sciraInfNet}, with +1 referring to positive regulation, -1 negative regulation, and 0 no regulation.
#'
#' @return A vector storing the activity scores of all TFs from user specified sample in the data matrix.
#'
#' @details The user needs to input the index of the sample in the data matrix from which \code{InferTFactPRL} estimates TF activity. \code{InferTFactPRL} accomplishes the work by regressing the expression profile of TF target genes against the binding profile of this TF on these genes. The output t-statistics are taken as the TF activity score.
#'
#' @importFrom stats lm
#'
#' @examples
#' # TF regulatory network with 100 genes and 5 TFs.
#' regnet.m <- matrix(sample(c(-1, 0, 1), 500, replace = TRUE), nrow = 100)
#' # gene expression data matrix with 100 genes and 10 samples.
#' tmp.m <- matrix(rnorm(1000), nrow = 100)
#' # TF activity score
#' TFact.v <- SCIRA:::InferTFactPRL(1, tmp.m, regnet.m)

InferTFactPRL <- function(idx, tmp, regnet) {
  exp.v <- tmp[, idx]
  act.v <- apply(regnet, 2, function(tmp.v) {
    lm.o <- lm(exp.v ~ tmp.v)
    act <- summary(lm.o)$coeff[2, 3]
    return(act)
  })

  return(act.v)
}
