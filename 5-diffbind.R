## Set working directory
setwd("C:/path/to/your/working/directory")

## Load libraries
library(DiffBind)
library(tidyverse)
library(BiocParallel)
library(ggplot2)

## Reading in peaksets and associated metadata
# Note: control first and sample at the end
samples <- read.csv("meta/samplesheet_stmreps.txt")
dbObj <- dba(sampleSheet=samples)  #build a consensus set
dbObj #have a look at consensus set

samples <- read.csv("meta/samplesheet_onlywt.txt")
dbObj2 <- dba(sampleSheet = samples)
dbObj2

samples <- read.csv("meta/samplesheet_noabs.txt")
dbObj3 <- dba(sampleSheet = samples)
dbObj3

samples4 <- read.csv("meta/samplesheet_stmvwt.txt")
dbObj4 <- dba(sampleSheet = samples4)
dbObj4

## Compute count information for each of the peaks in the consensus set
#Note: increase memory to 500 000 with memory.limit(size = 500000) before executing
dbObj <- dba.count(dbObj, bParallel = FALSE, bUseSummarizeOverlaps = TRUE)
dbObj

dbObj2 <- dba.count(dbObj2, bParallel = FALSE, bUseSummarizeOverlaps = TRUE)
dbObj2

dbObj3 <- dba.count(dbObj3, bParallel = FALSE, bUseSummarizeOverlaps = TRUE)
dbObj3

dbObj4 <- dba.count(dbObj4, bParallel = FALSE, bUseSummarizeOverlaps = TRUE)
dbObj4

dba.plotPCA(dbObj4, attributes=DBA_FACTOR, label = DBA_ID)
plot(dbObj4)

## Establishing a contrast (which samples we want to compare to one another)
dbObj.contrast <- dba.contrast(dbObj, categories = DBA_CONDITION, minMembers = 2)
dbObj.contrast <- dba.analyze(dbObj.contrast, method = DBA_ALL_METHODS)

dbObj.wt <- dba.contrast(dbObj2, categories = DBA_CONDITION, minMembers = 2)
dbObj.wt <- dba.analyze(dbObj.wt, method = DBA_ALL_METHODS)

dbObj.stmvwt <- dba.contrast(dbObj4, categories = DBA_FACTOR, minMembers = 2)
dbObj.stmvwt <- dba.analyze(dbObj.stmvwt, method = DBA_ALL_METHODS)

## Examining which method is better and what peak distribution is between samples
dba.show(dbObj.contrast, bContrasts = T)
dba.plotPCA(dbObj.contrast, contrast = 1, method = DBA_DESEQ2, attributes = DBA_CONDITION, 
            label = DBA_ID)
dba.plotPCA(dbObj.contrast, contrast = 1, method = DBA_EDGER, attributes = DBA_CONDITION, 
            label = DBA_ID)
dba.plotVenn(dbObj.contrast, contrast=1, method = DBA_ALL_METHODS)

dba.plotMA(dbObj.contrast, method = DBA_DESEQ2)
dba.plotMA(dbObj.contrast, bXY = TRUE)

pvals <- dba.plotBox(dbObj.contrast)


## Extracting results
res_deseq <- dba.report(dbObj.contrast, method = DBA_DESEQ2, contrast = 1, th = 1)
# note: result is a GRanges object

res_deseq_wt <- dba.report(dbObj.wt, method = DBA_DESEQ2, contrast = 1, th = 1)

## Write out file
out <- as.data.frame(res_deseq) #convert GRanges object to dataframe
write.table(out, file = "results/STM_IPvNOAB_DeSeq2.txt", sep = "\t", quote = F, row.names = F)

out.wt <- as.data.frame(res_deseq_wt)
write.table(out.wt, file = "results/STMipvWTip_DeSeq2.txt", sep = "\t", quote = F, row.names = F)

## Create bed files for ip-enriched keeping only significant peaks (p < 0.05)
ip_enrich <- out %>%
  filter(FDR < 0.05 & Fold > 0) %>%
  select(seqnames, start, end)

ip_enrich_v2 <- out %>%
  filter(FDR < 0.05 & Fold < 0)

ip_enrich2 <- out %>%
  filter(FDR < 0.05 & Fold < -1) %>%
  select(seqnames, start, end)
ip_enrich2_v2 <- out %>%
  filter(FDR < 0.05 & Fold < -1)


noab_enrich <- out %>%
  filter(FDR < 0.05 & Fold >0) %>%
  select(seqnames, start, end)

stmvwt_enriched <- out.wt %>%
  filter(FDR < 0.05 & Fold < 0) %>%
  select(seqnames, start, end)
stmvwt_enriched_stringent <- out.wt %>%
  filter(FDR < 0.05 & Fold < -1) %>%
  select(seqnames, start, end)
wtvstm_enriched <- out.wt %>%
  filter(FDR < 0.05 & Fold > 0) %>%
  select(seqnames, start, end)

## Write out bed file
write.table(ip_enrich_v2, file = "results/stm_ip_enriched.bed", sep = "\t", quote = F, 
            row.names = F, col.names = F)

write.table(ip_enrich2_v2, file = "results/stm_ip_enriched_stringent.bed", sep = "\t", quote = F, 
            row.names = F, col.names = F)

write.table(noab_enrich, file = "results/stm_noab_enriched.bed", sep = "\t", quote = F, 
            row.names = F, col.names = F)

write.table(stmvwt_enriched, file = "results/stmvwt_ip_enriched.bed", sep = "\t", quote = F,
            row.names = F, col.names = F)

write.table(stmvwt_enriched_stringent, file = "results/stmvwt_ip_enriched_stringent.bed",
            sep = "\t", quote = F, row.names = F, col.names = F)

write.table(wtvstm_enriched, file = "results/wtvstm_ip_enriched.bed",
            sep = "\t", quote = F, row.names = F, col.names = F)
