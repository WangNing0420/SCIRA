#' @title Compute partial correlation coefficient
#'
#' @description \code{ComputePCOR} computes partial correlation between transcription regulators and their targets.
#'
#' @param idx Numeric, the index of target genes in the row of a binarized network \code{selbinNET.m}.
#' @param mapTG.idx A vector of indexes output from \code{match} function when mapping the target genes in the network to rows of the original expression data matrix.
#' @param mapTF.idx A vector of indexes output from \code{match} function when mapping the transcription factors (TFs) in the network to rows of the original expression data matrix.
#' @param selbinNET.m A matrix, the binarized network with rows referring to TF target genes, and columns to TFs (the regulators). 0s means no regulation between TF-gene, while 1s means significant regulation (either positive or negative).
#' @param exp A matrix, the original gene expression data matrix across different tissue types with rows referring to genes, columns to samples.
#'
#' @return A matrix with partial correlation coefficients between TF targets and their regulators.
#'
#' @details This function is used to calculate partial correlations between a TF target gene and all its regulators.
#'
#' User needs to first provide the index (\code{idx}) of the gene in the binarized matrix, then all its regulators are detected from the network. Then the expression profile of this gene and its regulators are retrieved from the expression data matrix (\code{exp}) with \code{mapTG.idx} and \code{mapTF.idx} indexes which indicates the number of rows they are in \code{exp} (from \code{match} function).
#'
#' The output could be:
#' 1) a matrix storing the partial correlation coefficients;
#' 2) NULL if the selected gene has only one regulator according to the binarized network.
#'
#' @importFrom stats cor

ComputePCOR <- function(idx,mapTG.idx,mapTF.idx,selbinNET.m,exp.m){
  g <- idx
  reg.idx <- which(selbinNET.m[g,]==1)
  if(length(reg.idx)>=2){
    tmp.idx <- c(mapTG.idx[g],mapTF.idx[reg.idx])
    cor.m <- cor(t(exp.m[tmp.idx,]))
    pcor.m <- cor2pcor(cor.m)
  }
  else{
    pcor.m <- NULL
  }
  return(pcor.m)
}