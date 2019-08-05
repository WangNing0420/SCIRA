# SCIRA-package
Single Cell Inference of Regulatory Activity

`SCIRA` leverages the power of large-scale bulk RNA-Seq datasets to infer regulatory networks encoded by transcription factor regulon pairs, and subsequently uses these regulons to estimate regulatory activity of the transcription factors in single cells. It encompasses two main steps:

1. Construction of a tissue-specific transcription factor regulatory network, consisting of transcription factors that are more highly expressed in the user-specified tissue type (the 'tissue type of interest') compared to other tissue types, plus an associated set of high-confidence downstream targets (the "regulon").
2. Estimation of transcription factor activity at single-cell resolution using a corresponding scRNA-seq dataset representative of the same tissue-type.


## Usage
#### Inferring tissue-specific network

```{r eval=FALSE}
net.o <- sciraInfNet(data=data.m, tissue=colnames(data.m), toi = "Lung", cft = c("Blood","Spleen"),
         TFs = TFeid, sdth = 0.25, sigth = 0.05, capth=0.01, pcorth = 0.2, degth = c(0.05, 0.05),
         lfcth = c(log2(1.5), 0), minNtgts = 5, ncores = 1)
```
** Note: `data.m` must be an intra-sample log-normalized bulk-tissue mRNA expression (RNA-Seq) dataset, for instance the dataset from GTEX.

#### Estimating transcription factor activity
```{r eval=FALSE}
TFact <- sciraRegAct(data = data.m, regnet = net.o$netTOI, norm = "z", ncores = 1)
```
#### Installation
An easy way to install SCIRA is by facilitating the devtools R package.
```{r eval=FALSE}
#install.packages("devtools")
library(devtools)
install_github("WangNing0420/SCIRA", build_opts = c("--no-resave-data", "--no-manual"))
```

## Getting started
The SCIRA package contains a tutorial showing people how to implement SCIRA in their work. The tutorial can be found in the package-vignette:

```{r eval=FALSE}
library(SCIRA)

vignette("SCIRA")
```

## References

Chen Y, Widschwendter M, and Teschendorff AE. 2017. “Systems-Epigenomics Inference of Transcription Factor Activity Implicates Aryl-Hydrocarbon-Receptor Inactivation as a Key Event in Lung Cancer Development.” Genome Biol 18:236.

Wang, N., & Teschendorff, AE. 2019.  “Leveraging high-powered RNA-Seq datasets to improve inference of regulatory activity in single-cell RNA-Seq data.“ BioRxiv, 553040. 
