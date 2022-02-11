## Load libraries
library(ChIPpeakAnno)
library(ChIPseeker)
library(GenomicFeatures)
library(biomaRt)
library(gprofiler2)
library(tidyverse)

## Create tx database for Araport annotation (Note: example for A. thaliana. It will be a different mart if it's not a plant)
txdb <- makeTxDbFromBiomart(biomart = "plants_mart", host = "plants.ensembl.org",
                            dataset = "athaliana_eg_gene")

## Load data
samplefiles <- list.files("C:/path/to/your/directory", pattern = ".bed", full.names = T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("sample", "control")

## Annotation
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb = txdb, tssRegion = c(-5000,5000), 
                       verbose = FALSE)


## Coverage of genome for each sample
peak_sample <- readPeakFile(samplefiles[[1]])
coverage_sample <- covplot(peak_sample, weightCol = peak_sample$V5)
coverage_sample
ggsave("coverage_sample.pdf", plot = coverage_sample, path = "C:/path/to/your/directory", dpi = 800)
ggsave("coverage_sample.png", plot = coverage_sample, path = "C:/path/to/your/directory", dpi = 800)


peak_control <- readPeakFile(samplefiles[[2]])
coverage_control <- covplot(peak_control, weightCol = peak_control$V5)
coverage_control
ggsave("coverage_control.pdf", plot = coverage_control, path = "C:/path/to/your/directory", dpi = 800)
ggsave("coverage_control.png", plot = coverage_control, path = "C:/path/to/your/directory", dpi = 800)

promoter <- getPromoters(TxDb = txdb, upstream = 5000, downstream = 5000)

tagMatrix_sample <- getTagMatrix(peak_sample, windows = promoter)
avgPlot_sample <- plotAvgProf(tagMatrix_sample, xlim = c(-5000, 5000), conf = 0.95, resample = 1000, xlab = "Genomic Region (5' -> 3')", ylab = "Read Count Frequency")
ggsave("avgPlot_sample.pdf", plot = avgPlot_sample, path = "C:/path/to/your/directory", dpi = 800)
ggsave("avgPlot_sample.png", plot = avgPlot_sample, path = "C:/path/to/your/directory", dpi = 800)

tagMatrix_control <- getTagMatrix(peak_control, windows = promoter)
avgPlot_control <- plotAvgProf(tagMatrix_control, xlim = c(-5000, 5000), conf = 0.95, resample = 1000, xlab = "Genomic Region (5' -> 3')", ylab = "Read Count Frequency")
ggsave("avgPlot_control.pdf", plot = avgPlot_control, path = "C:/path/to/your/directory", dpi = 800)
ggsave("avgPlot_control.png", plot = avgPlot_control, path = "C:/path/to/your/directory", dpi = 800)


## Visualization of genomic feature representation
comparison_bar <- plotAnnoBar(peakAnnoList)

comparison_to_TSS <- plotDistToTSS(peakAnnoList, title = "Distribution of TF-binding loci relative to TSS")

bar2_1 <- plotDistToTSS(peakAnnoList[[1]], title = "Distribution of TF-binding loci relative to TSS")
ggsave("bar2_1.pdf", plot = bar2_1, path = "C:/Users/sbitl2/Documents", dpi = 800)
ggsave("bar2_1.png", plot = bar2_1, path = "C:/Users/sbitl2/Documents", dpi = 800)

bar2_2 <- plotDistToTSS(peakAnnoList[[2]], title = "Distribution of TF-binding loci relative to TSS")
ggsave("bar2_2.pdf", plot = bar2_2, path = "C:/Users/sbitl2/Documents", dpi = 800)
ggsave("bar2_2.png", plot = bar2_2, path = "C:/Users/sbitl2/Documents", dpi = 800)

peak1 <- peakAnnoList[[1]]
plotAnnoPie(peak1)

peak2 <- peakAnnoList[[2]]
plotAnnoPie(peak2)


## Retrieve annotation
sample_annot <- data.frame(peakAnnoList[["sample"]]@anno)
control_annot <- data.frame(peakAnnoList[["control"]]@anno)

## Write down gene names in a table
sample.df <- data.frame(peakAnnoList[["sample"]])
write.table(sample.df, "results/beds/sample_peaks.txt", sep = "\t")
control.df <- data.frame(peakAnnoList[["control"]])
write.table(control.df, "results/beds/control_peaks.txt", sep = "\t")

            
## Functional enrichment (GO)
gconv <- gconvert(query = sample.df$geneId, organism = "athaliana", #Note: change to your organism of interest 
                  target = "ENSG", mthreshold = Inf, filter_na = TRUE)
colnames(gconv)[2] <- "geneId"
sample.annot <- merge(sample.df, gconv, by='geneId')
write.table(sample.annot, "results/beds/annotated_sample_peaks.txt", sep = "\t")


gconv <- gconvert(query = control.df$geneId, organism = "athaliana", target = "ENSG",
                  mthreshold = Inf, filter_na = TRUE)
colnames(gconv)[2] <- "geneId"
control.annot <- merge(control.df, gconv, by='geneId')
write.table(control.annot, "results/beds/annotated_control_peaks.txt", sep = "\t")

