#' @title Infer TF activity from gene expression/ DNA methylation profile
#'
#' @description \code{sciraRegAct} calculates TF activity scores in user input data set. It could be a single cell gene expression dataset.
#'
#' @param data A gene expression data matrix, with rows referring to genes and columns to samples.
#' @param regnet A matrix, the regulatory network inferred from \code{sciraInfNet} function.
#' @param norm A character indicating the method used to normalize your input data set, "c" for "centering"; "z" for "z-score normalization".
#' @param ncores A numeric, the number of cores to use. See \code{mclapply}.
#'
#' @return A matrix of TF activity score with rows referring to TFs, columns to samples.
#'
#' @details \code{sciraRegAct} is one of the two main functions in \code{SCIRA} package. It takes the output regulatory network from \code{sciraInfNet} as input, and computes the activity of all TFs in this network from user provided \code{data}.
#'
#' The \code{data} matrix should be single cell gene expression data, with rows are genes and columns are samples. Duplicated row names are not allowed, so you should average the these rows before running \code{sciraRegAct}.
#'
#' Note that it's very important that you use the same gene identifier through out the whole analysis.
#'
#' @import parallel
#' @importFrom stats cor
#' @importFrom stats lm
#' @importFrom stats model.matrix
#' @importFrom stats pbinom
#' @importFrom stats pnorm
#' @importFrom stats sd
#'
#' @export sciraRegAct
#'
#' @examples
#' # estimate TF activity
#' data("scLiver")
#' data("netLiver")
#' TFact <- sciraRegAct(scLiver,regnet=netLiver,norm="z",ncores=1)

sciraRegAct <- function(data, regnet, norm = c("c", "z"), ncores = 4) {


  if (is.vector(data)) {
    inter <- intersect(names(data), rownames(regnet))
    data <- data[inter]
    regnet <- regnet[inter, ]
    actTF <- InferTFact(data, regnet)
    names(actTF) <- colnames(regnet)
  }
  else if (is.matrix(data)) {
    common.v <- intersect(rownames(data),rownames(regnet));
    match(common.v,rownames(data)) -> map1.idx;
    match(common.v,rownames(regnet)) -> map2.idx;        
    ndata <- data[map1.idx,] - rowMeans(data[map1.idx,]);
    
    if(norm=="z"){
      sd.v <- apply(data[map1.idx,],1,sd);
      nz.idx <- which(sd.v>0);
      z.idx <- which(sd.v==0);
      ndata <- data[map1.idx,];       
      ndata[nz.idx,] <- (data[map1.idx[nz.idx],] - rowMeans(data[map1.idx[nz.idx],]))/sd.v[nz.idx];
      ndata[z.idx,] <- 0;
    }
    idx.l <- as.list(1:ncol(data));
    prl.o <- mclapply(idx.l,InferTFactPRL,ndata,regnet[map2.idx,],mc.cores=ncores);
    actTF <- matrix(unlist(prl.o), nrow = ncol(regnet), ncol=  length(prl.o), byrow = FALSE)
    rownames(actTF) <- colnames(regnet)
    colnames(actTF) <- colnames(data)
  }

  return(actTF)
}
