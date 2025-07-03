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
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(tidyverse)
#BiocManager::install("ReactomePA")
library(ReactomePA)
library(ggpubr)


```



```{r K562 ChIP peaks coverage plot, echo=FALSE, message=FALSE, warning=FALSE}

k562peaks <-lapply(c("pS5_ENCFF060MMW.bed", "pS2_ENCFF951KHS.bed", "BRCA1_ENCFF620FIH.bed"), readPeakFile, header=FALSE)

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

K562tagMatrix <- lapply(k562peaks, getTagMatrix, windows=promoter)

K562compareTSS <- plotAvgProf(K562tagMatrix, xlim=c(-3000, 3000))  +
  scale_color_manual(values = c("blue", "red", "green"))

K562compareTSS                    

```

```{r HepG2 ChIP peaks coverage plot, echo=FALSE, message=FALSE, warning=FALSE}

HepG2peaks <-lapply(c("pS5_HepG2_Myers.bed", "pS2_HepG2_Snyder.bed", "BRCA1_HepG2_Myers.bed"), readPeakFile, header=FALSE)

HepG2tagMatrix <- lapply(HepG2peaks, getTagMatrix, windows=promoter)

HepG2compareTSS <- plotAvgProf(HepG2tagMatrix, xlim=c(-3000, 3000))  +
  scale_color_manual(values = c("blue", "red", "green"))

HepG2compareTSS                    

```


```{r HepG2 ChIP peak annotation, echo=FALSE, message=FALSE, warning=FALSE}


HepG2peakAnnoList <- lapply(HepG2peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), annoDb="org.Hs.eg.db")


plotAnnoBar(HepG2peakAnnoList)

```


```{r K562 ChIP peak annotation, echo=FALSE, message=FALSE, warning=FALSE}


K562peakAnnoList <- lapply(k562peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), annoDb="org.Hs.eg.db")


plotAnnoBar(K562peakAnnoList)

```

```{r Distribution of transcription factor-binding loci relative to TSS, echo=FALSE, message=FALSE, warning=FALSE}


plotDistToTSS(HepG2peakAnnoList)

plotDistToTSS(K562peakAnnoList)

```

```{r Functional enrichment analysis, echo=FALSE, message=FALSE, warning=FALSE}


##HepG2 Global

HepG2peakannolistDF<- lapply(HepG2peakAnnoList, as.data.frame)

HepG2peakAnnoListgeneId<- lapply(HepG2peakannolistDF, pull, geneId)

HepG2enrichlist<- lapply(HepG2peakAnnoListgeneId, enrichPathway)
names(HepG2enrichlist)<-c(names(HepG2peakAnnoListgeneId))

##HepG2 Promoter 

HepG2peakAnnoListPromoter <- lapply(names(HepG2peakAnnoList), function(x) filter(HepG2peakAnnoList[[x]], 
                                                                         grepl("annotation", "Promoter")))
                               
HepG2peakAnnoListPromotergeneId<- lapply(HepG2peakAnnoListPromoter, pull, geneId)

HepG2enrichlistPromoter<- lapply(HepG2peakAnnoListPromotergeneId, enrichPathway)

names(HepG2enrichlistPromoter)<-c(names(HepG2peakannolistDF))

##Tables of associated pathways 

#Global

lapply(enrichlist, as.data.frame) %>%
  kbl(caption = htmltools::tags$caption("Global Reactome Pathway Enrichment", style = "color:black")) %>%
  kable_styling(bootstrap_options = "striped", full_width = T, html_font = "Cambria")

#Promoter

lapply(enrichlistPromoter, as.data.frame) %>%
  kbl(caption = htmltools::tags$caption("Promoter Reactome Pathway Enrichment", style = "color:black")) %>%
  kable_styling(bootstrap_options = "striped", full_width = T, html_font = "Cambria")

##Plotting


#BRCA1 global vs promoter only associated genes dmso

ggarrange(dotplot(HepG2enrichlist[[3]],  showCategory=10))

#BRCA1 global vs promoter only associated genes etoposide

dotplot(HepG2enrichlistPromoter[[3]])


#U2OSRNAPII global vs promoter only associated genes etoposide
ggarrange(dotplot(enrichlist[["U2OSRNAPIIetop_callpeaknoabC"]], showCategory=11), dotplot(enrichlistPromoter[["U2OSRNAPIIetop_callpeaknoabC"]], showCategory=11), nrow = 2, ncol = 1, align = 'v')

#U2OSRNAPII global vs promoter only associated genes dmso

ggarrange(dotplot(enrichlist[["U2OSRNAPIIdmso1_callpeaknoabC"]], showCategory=11), dotplot(enrichlistPromoter[["U2OSRNAPIIdmso1_callpeaknoabC"]], showCategory=11), nrow = 2, ncol = 1, align = 'v')

```




