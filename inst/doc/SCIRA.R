## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(fig.pos="h", out.extra='', fig.align="center")
library(SCIRA)
data("GeneExp")
data.m <- GeneExp
dataVAL.m <- GeneExp
data("TFeid")

## ----warning=F-----------------------------------------------------------
net.o <- sciraInfNet(data=data.m, tissue=colnames(data.m), toi = "Lung", cft = c("Blood","Spleen"),
         TFs = TFeid, sdth = 0.25, sigth = 0.05, capth = 0.01, pcorth = 0.2, degth = c(0.05, 0.05),
         lfcth = c(log2(1.5), 0), minNtgts = 5, ncores = 1)
## Note: `data.m` must be an intra-sample log-normalized bulk-tissue mRNA expression (RNA-Seq) dataset, like the dataset from GTEX.
## Parameters used here are not recommended. See "?sciraInfNet" for more info.

## ------------------------------------------------------------------------
act <- sciraRegAct(data = data.m, regnet = net.o$netTOI, norm = "z", ncores = 1)

## ----warning=F-----------------------------------------------------------
net.o <- sciraInfNet(data=data.m, tissue=colnames(data.m), toi = "Lung", cft = c("Blood","Spleen"),
         TFs = TFeid, sdth = 0.25, sigth = 0.05, capth= 0.01, pcorth = 0.2, degth = c(0.05, 0.05),
         lfcth = c(log2(1.5), 0), minNtgts = 3, ncores = 1)

## ----out.width = 450, fig.retina = NULL, echo=F--------------------------
knitr::include_graphics("Figures/LiverNet.png")

## ------------------------------------------------------------------------
# estimating transcription factor activity in data.m
TFact <- sciraRegAct(data.m,regnet=net.o$netTOI,norm="z",ncores=1)
## Note: `data.m` here should be a normalized scRNA-seq gene expression data set.

## ----out.width = 200, fig.retina = NULL, echo=F--------------------------
knitr::include_graphics("Figures/TFAct.png")

## ----out.width = 500, fig.retina = NULL, echo=F--------------------------
knitr::include_graphics("Figures/Hnf1a.png")

## ----sessionInfo---------------------------------------------------------
sessionInfo()

