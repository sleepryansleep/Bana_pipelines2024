#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
options(RCurlOptions=list(followlocation=TRUE, postredir=2L))

#load required packages
require("groHMM")
library("IRanges")
library("GenomicRanges")
library("rtracklayer")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
library("reshape2")
library("ggplot2")

#USAGE:
#Rscript /endosome/work/GCRB/s161282/scripts/box_plot.r \
#<1. project directory> \
#<2. bam1, WT> \
#<3. bam1 name> \
#<4. bam2, Test> \
#<5. bam2 name> \
#<6. peaks1> \
#<7. peaks1 name> \
#<8. file out name>

dir <- args[1]
setwd(dir)

#read bams
bam1.name <- args[3]
bam1.bam <- readGAlignments(file.path(args[2]))
bam1.gr <- granges(bam1.bam)
bam1 <- bam1.gr
#find library size
library_bam1 <- NROW(bam1)

bam2.name <- args[5]
bam2.bam <- readGAlignments(file.path(args[4]))
bam2.gr <- granges(bam2.bam)
bam2 <- bam2.gr
#find library size
library_bam2 <- NROW(bam2)

#read peaks
peaks1.name <- args[7]
peaks1 <- import(file.path(args[6]), format = "BED")

#Calculate RPKM, = numReads / ( geneLength/1000 * totalNumReads/1,000,000 )

rpkm_bam1_peaks1 <- countOverlaps(peaks1, bam1) / (width(peaks1)/1000 * library_bam1/1000000)
sts1 <- boxplot.stats(rpkm_bam1_peaks1)$stats
bam1_peaks1 <- paste(bam1.name, peaks1.name, sep = "_")
rpkm_bam1_peaks1_df <- data.frame(rpkm_bam1_peaks1)
colnames(rpkm_bam1_peaks1_df) <- bam1_peaks1
rpkm_bam1_peaks1_m <- melt(data = rpkm_bam1_peaks1_df, measure.vars = bam1_peaks1)

rpkm_bam2_peaks1 <- countOverlaps(peaks1, bam2) / (width(peaks1)/1000 * library_bam2/1000000)
sts2 <- boxplot.stats(rpkm_bam2_peaks1)$stats
bam2_peaks1 <- paste(bam2.name, peaks1.name, sep = "_")
rpkm_bam2_peaks1_df <- data.frame(rpkm_bam2_peaks1)
colnames(rpkm_bam2_peaks1_df) <- bam2_peaks1
rpkm_bam2_peaks1_m <- melt(data = rpkm_bam2_peaks1_df, measure.vars = bam2_peaks1)

outprefix <- args[8]

rpkm_all <- rbind(rpkm_bam1_peaks1_m,rpkm_bam2_peaks1_m)
colnames(rpkm_all) <- c("Peaks","RPKM")
rpkm_plot <- ggplot(rpkm_all, aes(x = Peaks, y = RPKM))

#write pdf
out_box <- paste(outprefix, "box_plot.pdf", sep = "_")
bp <- rpkm_plot + geom_boxplot(aes(fill = Peaks), notch=TRUE, outlier.colour = NA) + labs(title=outprefix) + theme(legend.position="none", axis.text.x=element_text(size=7)) + scale_fill_brewer(palette="Blues")
# scale y limits based on ylim1
sts_df <- cbind(sts1,sts2)
bp_scaled <- bp + coord_cartesian(ylim = c(0,max(sts_df)*1.05))
ggsave(filename=file.path(dir, out_box), plot=bp_scaled)

out_violin <- paste(outprefix, "violin_plot.pdf", sep = "_")
vp <- rpkm_plot + geom_violin(aes(fill = Peaks),trim = FALSE) + labs(title=outprefix) + scale_fill_brewer(palette="Blues") + theme(legend.position="none") + geom_boxplot(width = 0.2) + scale_y_continuous(trans = 'log10')
ggsave(filename=file.path(dir, out_violin), plot=vp)

out_cox <- paste(outprefix, "wilcox_test.txt", sep = "_")
cox_peaks1 <- wilcox.test(rpkm_bam1_peaks1,rpkm_bam2_peaks1,conf.int=TRUE)

coxs <- rbind(cox_peaks1$p.value,median(rpkm_bam1_peaks1),median(rpkm_bam2_peaks1))
colnames(coxs) <- c(peaks1.name)
rownames(coxs) <- c("p.value",bam1.name,bam2.name)
write.table(coxs, file=file.path(dir, out_cox), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

q()
