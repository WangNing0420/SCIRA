#' @title Do differential gene expression
#'
#' @description \code{LimmaFn} uses functions in \code{limma} package to easily compute the moderated t-statistics and p-values from differential gene expression tests comparing between different phenotypes even when sample size is small.
#'
#' @param pheno A vector of sample phenotypes. Sample phenotype in a scientific research could be treatment/control, normal/cancer or smoker/non-smoker. Different phenotypes could each be encoded as 0/1 when inputting to \code{LimmaFn}, for example, Normal-0; Cancer-1.
#' @param data A matrix, the unnormalized gene expression dataset, should be a numeric matrix, with rows referring to genes and columns to samples. In this matrix you could use any type of gene IDs, like Entrez ID, Ensembl ID, HUGO gene symbol... But make sure to use the same gene annotation through out your analysis.
#'
#' @return A table with rows for all genes (ranked by significance) and columns of log2 fold-change, average expression, moderated t-statistic, p-value, adjusted p-value (Benjamini–Hochberg procedure). The table is the output of \code{\link[limma]{topTable}} function.
#'
#' @details This function computes the moderated t-statistic for users using empirical Bayes method, it is especially useful when the sample size is too small to perform parametric tests.
#'
#' Given a gene expression data matrix and a vector indicating sample phenotype, \code{LimmaFn} first fits a linear model using \code{\link[limma]{lmFit}}, then it refits the model and do comparisons between any two different phenotypes with \code{\link[limma]{contrasts.fit}}, finally it estimates moderated t-statistics for each comparison from the fitted model using empirical Bayes method (\code{\link[limma]{eBayes}}) and output the result from the \code{\link[limma]{topTable}} function.
#'
#' Note that doing the \code{contrasts.fit} step will not make a difference if you do comparison between two different sample status (treatment/control). However, When there are more than two sample status in your data set, this step will do comparison between every two status. And resulted summary tables will be stored in a list.
#'
#' @seealso \code{\link[limma]{lmFit}} for fitting a linear model, \code{\link[limma]{contrasts.fit}} for refitting, \code{\link[limma]{eBayes}} for Bayes method, \code{\link[limma]{topTable}} for the output table.
#'
#' @importFrom stats cor
#' @importFrom stats lm
#' @importFrom stats model.matrix
#' @importFrom stats pbinom
#' @importFrom stats pnorm
#' @importFrom stats sd
#'
#' @export LimmaFn
#'
#' @examples
#' # prepare the phenotype info ("C"-control; "T"-treatment)
#' pheno.v <- c("C","C","C","T","T","T")
#' # prepare your normalized data matrix.
#' data.m <- matrix(rnorm(120),nrow=20,ncol=6)
#' # run function
#' lim.o <- LimmaFn(pheno.v, data.m)

LimmaFn <- function(pheno.v,data.m){
  
  
  
  ### construct model matrix
  
  sampletype.f <- as.factor(pheno.v);
  
  design.sample <- model.matrix(~0 + sampletype.f);
  
  colnames(design.sample) <- levels(sampletype.f);
  
  sampletypes.v <- levels(sampletype.f);
  
  
  
  ### do linear model fit
  
  lmf.o <- lmFit(data.m,design.sample);
  
  
  
  ### construct contrast matrix
  
  ntypes <- length(levels(sampletype.f));
  
  ncomp <- 0.5*ntypes*(ntypes-1);
  
  cont.m <- matrix(0,nrow=ncol(design.sample),ncol=ncomp);
  
  tmp.v <- vector();
  
  c <- 1;
  
  for(i1 in 1:(ntypes-1)){
    
    for(i2 in (i1+1):ntypes){
      
      cont.m[i1,c] <- -1;
      
      cont.m[i2,c] <- 1;
      
      tmp.v[c] <- paste(sampletypes.v[i2],"--",sampletypes.v[i1],sep="");
      
      c <- c+1;
      
    }
    
  }
  
  rownames(cont.m) <- sampletypes.v; # sampletype.v determined separately
  
  colnames(cont.m) <- tmp.v;
  
  
  
  ### do linear model to contrasts
  
  lmf2.o <- contrasts.fit(lmf.o,cont.m);
  
  
  
  ### empirical Bayesian estimation of differentially expressed genes (DEGs)
  
  bay.o <- eBayes(lmf2.o);
  
  
  
  ### build ranked list of DEGs for each comparison
  
  top.lm <- list();
  
  for(c in 1:ncol(cont.m)){
    
    top.lm[[c]] <- topTable(bay.o,coef=c,adjust="fdr",number=nrow(data.m));
    
  }
  
  
  
  return(list(top=top.lm,cont=cont.m));
  
  
  
} ## end of limma function