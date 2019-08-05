#' RNA-seq data matrix from GTEx 
#' 
#' GeneExp is a gene expression data matrix with rows (467 genes) and 520 columns. This data matrix is a subset of GTEx RNA-seq data which is a very big dataset. In this subset we only randomly chose 20 samples from each tissue type with more than 50 samples in the original GTEx data set.
#' 
#' 
#'
#' 
#' @docType data
#' @keywords datasets
#' @name GeneExp
#' @usage data("GeneExp")
#' @format A matrix with 520 rows and 467 columns
#' @source 
#' 
#' GETX URL:\href{https://www.gtexportal.org/home/}{
#' https://www.gtexportal.org/home/}.
#' 
NULL

#' Transcription factor genes annotated using Entrez gene ID 
#' 
#' TFeid is a vector of transcription factor Entrez IDs. This list of TFs is from MSigDB with 1385 TFs in it. 
#' 
#'
#' 
#' @docType data
#' @keywords datasets
#' @name TFeid
#' @usage data("TFeid")
#' @format A vecto with 1385 observations. chr [1:1385] "9684" "59269" "2965" "6935" "9247" "3661" "6239" "11151" ...
#' @Examples 
#' 
#' data(TFeid)
#'## maybe str(TFeid) ; plot(TFeid) ...
#' 
NULL

#' Liver-specific regulatory network
#' 
#' Liver-specific regulatory network inffered with sciraInfNet function, from GETx dataset. 
#' 
#'
#' 
#' @docType data
#' @keywords datasets
#' @name netLiver
#' @usage data("netLiver")
#' @format A matrix, each column represent a TFs, and rows are target genes.
#' @Examples 
#' 
#' data(netLiver)
#' dim(netLiver)
#' 
NULL

#' scRNA-seq data matrix 
#'  
#' scRNA-seq dataset with 70 samples and 12875 genes, from Yang et al. 2017
#' 
#'
#' 
#' @docType data
#' @keywords datasets
#' @name scLiver
#' @usage data("scLiver")
#' @format A matrix with 12875 rows and 70 coloums.
#' @references 
#' Li Yang  Wei-Hua Wang  Wei-Lin Qiu  Zhen Guo  Erfei Bi  Cheng-Ran Xu. 
#' \emph{A single-cell transcriptomic analysis reveals precise pathways and regulatory mechanisms underlying hepatoblast differentiation}
#' Hepatology 11 (2017): 29353.
#' doi:\href{ https://doi.org/10.1002/hep.29353}{
#' 10.1002/hep.29353}.
#' 
NULL

