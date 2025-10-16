# BIOS60132_Module2
A common task when working with transcriptomic data is the identification of differentially expressed (DE) genes or tags between groups. In this lesson students will learn how to perform biostatistical analysis in the R programming language with the DESeq2 or edgeR packages.

## Learning Goals
- Become comfortable working with quantified transcriptomic data (counts)
- Learn what contrasts are possible based on the experimental design
- Be able to perform single or multi factor analysis with DESeq2
- Be able to summarize and visualize the DE analysis results

## Prerequisites
- This lesson is designed for students who are unfamiliar with statistical analysis.
- Participants are expected to be comfortable working in the R programming language. 

## Primary Tools
- [R software environment](https://cran.rstudio.com/)
- [RStudio desktop software](https://libcal.library.nd.edu/event/9797081)

## Necessary R Packages
```
packageList <- c("ggplot2", "ggplotify", "pheatmap", "RColorBrewer", "ghibli", "BiocManager")
biocList <- c("DESeq2", "apeglm")
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
newBioc <- biocList[!(biocList %in% installed.packages()[,"Package"])]
if(length(newPackages)){
  install.packages(newPackages)
}
if(length(newBioc)){
  BiocManager::install(newBioc)
}
```


