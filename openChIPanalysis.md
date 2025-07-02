---
title: "Open Source ChIP Analaysis: BRCA1 and RNAPII phosphorylation"
author: "Aiola Stoja"
date: "`r Sys.Date()`"
output:
  html_document:
    theme: lumen
    code_folding: hide
    toc: yes
    toc_float:
      collapsed: yes
      smooth_scroll: yes
  pdf_document:
    toc: yes
---

```{r Packages, message=FALSE, warning=FALSE, include=FALSE}

library(BiocManager)
#BiocManager::install("ChIPseeker")
#BiocManager::install("GO.db")
library(ChIPseeker)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#BiocManager::install("clusterProfiler")
library(clusterProfiler)

```



```{r K562 ChIP peaks coverage plot, echo=FALSE, message=FALSE, warning=FALSE}

k562peaks <-lapply(c("pS5_ENCFF060MMW.bed", "pS2_ENCFF951KHS.bed", "BRCA1_ENCFF620FIH.bed"), readPeakFile, header=FALSE)

promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)

K562tagMatrix <- lapply(k562peaks, getTagMatrix, windows=promoter)

K562compareTSS <- plotAvgProf(K562tagMatrix, xlim=c(-1000, 1000))  +
  scale_color_manual(values = c("blue", "red", "green"))

K562compareTSS                    

```

```{r HepG2 ChIP peaks coverage plot, echo=FALSE, message=FALSE, warning=FALSE}

HepG2peaks <-lapply(c("pS5_HepG2.bed", "pS2_HepG2.bed", "BRCA1_HepG2.bed"), readPeakFile, header=FALSE)

HepG2tagMatrix <- lapply(HepG2peaks, getTagMatrix, windows=promoter)

HepG2compareTSS <- plotAvgProf(HepG2tagMatrix, xlim=c(-1000, 1000))  +
  scale_color_manual(values = c("blue", "red", "green"))

HepG2compareTSS                    

```