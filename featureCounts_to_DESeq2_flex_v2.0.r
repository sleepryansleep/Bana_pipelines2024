#!/usr/bin/env Rscript

#USAGE:
#Rscript ./featureCounts_to_DESeq2_flex_v2.0.r \
#<1. PATH/TO/WORKING/OUTDIRECTORY> \
#<2. PATH/TO/LIST_FILE.list> \
#<3. PATH/TO/COUNTS_FILE.txt>

args = commandArgs(trailingOnly = TRUE)
options(RCurlOptions=list(followlocation=TRUE, postredir=2L))

OUTDIR <- args[1]
OUTDIR
list <- args[2]
list
cfile <- args[3]
cfile
groupX <- args[4]
groupX
groupY <- args[5]
groupY
PREFIX <- args[6]
PREFIX
COUNTTHRESH <- args[7]
COUNTTHRESH
setwd(OUTDIR)
#exp_design <- eval(parse(text = args[XXXX]))

listfile <- file.path(list)
listfile
coldata <- as.data.frame(read.table(listfile, header=FALSE, row.names=1))
colnames(coldata) <- c("contrast")

countsfile <- file.path(cfile)
countsdata <- read.table(countsfile, header=TRUE, row.names=1)
countsdata <- countsdata[ ,6:ncol(countsdata)]
colnames(countsdata) <- c(row.names(coldata))
countsdata <- as.data.frame(countsdata)
countsdata <- round(countsdata)

#export counts file
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countsdata, colData = coldata, design = ~ contrast)
keep <- rowSums(counts(dds)) >= COUNTTHRESH
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized=TRUE)
out_name <- paste(PREFIX, "DESeq2_normalized_counts.txt", sep = "_")
out_file <- file.path(OUTDIR, out_name)
write.csv(norm_counts, file = out_file)

#differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = countsdata, colData = coldata, design = ~ contrast)
keep <- rowSums(counts(dds)) >= COUNTTHRESH
dds <- dds[keep,]

#dds$group <- factor(paste(dds$treatment, dds$cell, sep = "_"))
#design(dds) = ~ group

dds <- DESeq(dds)
resultsNames(dds)

comparison <- paste(groupY, groupX, sep = "_vs_")

res <- results(dds, contrast=c("contrast",groupY,groupX))
summary(res)
resDF <- as.data.frame(res)
out_name <- paste(PREFIX, comparison, "DESeq2.csv", sep = "_")
out_file <- file.path(OUTDIR, out_name)
write.csv(resDF, file = out_file)

library("ggplot2")
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("contrast"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
a <- max(abs(max(pcaData$PC1,pcaData$PC2)),abs(min(pcaData$PC1,pcaData$PC2)))
pca <- ggplot(pcaData, aes(PC1, PC2, color=contrast)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed(ratio = 1, xlim = c(-a,a), ylim = c(-a,a), expand = TRUE, clip = "on")
out_pca <- paste(PREFIX, "PCA.pdf", sep = "_")
ggsave(filename=file.path(OUTDIR, out_pca), plot=pca)

q()
