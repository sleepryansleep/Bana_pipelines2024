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
#<6. PEAKS> \
#<7. PEAKS name> \
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
PEAKS.NAME <- args[7]
PEAKS <- import(file.path(args[6]), format = "BED")

#Calculate RPK, = numReads / ( geneLength/1000 )

RPK_bam1_PEAKS <- countOverlaps(PEAKS, bam1) / (width(PEAKS)/1000)

RPK_bam2_PEAKS <- countOverlaps(PEAKS, bam2) / (width(PEAKS)/1000)

OUTPREFIX <- args[8]

COX.NAME <- paste(OUTPREFIX, "wilcox_test.txt", sep = "_")
COX.VAL <- wilcox.test(RPK_bam1_PEAKS,RPK_bam2_PEAKS,conf.int=TRUE)

COX.OUT <- cbind(PEAKS.NAME,COX.VAL$p.value,median(RPK_bam1_PEAKS),median(RPK_bam2_PEAKS))
COX.OUT
write.table(COX.OUT, file=file.path(dir, COX.NAME), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

q()
