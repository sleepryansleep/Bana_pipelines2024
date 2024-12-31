#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
options(RCurlOptions=list(followlocation=TRUE, postredir=2L))

library(pheatmap)
library(ggplot2)

OUTDIR <- dirname(args[1])
INNAME <- tools::file_path_sans_ext(basename(args[1]))
OUTNAME <- paste(INNAME,"heatmap.pdf", sep = "_")
ZscoreNormalize <- args[2]
print(paste0("ZscoreNormalize: ", ZscoreNormalize))
ClustCol <- as.logical(args[3])
print(paste0("ClustCol: ", ClustCol))
ClustRow <- as.logical(args[4])
print(paste0("ClustRow: ", ClustRow))
ColorPalette <- args[5]
print(paste0("ColorPalette: ", ColorPalette))

OUTDIR
INNAME
OUTNAME

setwd(OUTDIR)
DATA <- as.matrix(read.table(args[1], header = TRUE, row.names = 1))

paletteLength <- 100

if (ColorPalette == "") {
  myColor <- colorRampPalette(c("#1245BA", "white", "#AB2143"))(paletteLength)
} else if (ColorPalette == "BlueRed") {
  myColor <- colorRampPalette(c("#1245BA", "white", "#AB2143"))(paletteLength)
} else if (ColorPalette == "PurpleGreen") {
  myColor <- colorRampPalette(c("#533A7B", "white", "#00B85C"))(paletteLength)
} else if (ColorPalette == "Plasma") {
  myColor <- colorRampPalette(c("#0A2463", "#D8315B", "#FFC857"))(paletteLength)
} else {
  print("This color palette doesn't exist!")
}

min(DATA)
max(DATA)

if (ZscoreNormalize == "TRUE") {
  Calc_Zscore <- function(x){
    (x - mean(x)) / sd(x)
  }
  DATA_Z <- t(apply(DATA, 1, Calc_Zscore))
  ColorBreaks <- c(seq(min(DATA_Z), 0, length.out=ceiling(paletteLength/2) + 1),
  seq(max(DATA_Z)/paletteLength, max(DATA_Z), length.out=floor(paletteLength/2)))
  ColorBreaks
  H <- pheatmap(DATA_Z,
    scale = "none",
    cluster_cols = ClustCol,
    cluster_rows = ClustRow,
    fontsize = 20,
    fontsize_row = 4,
    color = myColor,
    breaks = ColorBreaks,
    border_color = "NA"
  )
} else if (ZscoreNormalize == "FALSE") {
  ColorBreaks <- c(seq(min(DATA), 0, length.out=ceiling(paletteLength/2) + 1),
  seq(max(DATA)/paletteLength, max(DATA), length.out=floor(paletteLength/2)))
  ColorBreaks
  H <- pheatmap(DATA,
    scale = "none",
    cluster_cols = ClustCol,
    cluster_rows = ClustRow,
    fontsize = 20,
    fontsize_row = 4,
    color = myColor,
    breaks = ColorBreaks,
    border_color = "NA"
  )
}

ggsave(filename=file.path(OUTDIR, OUTNAME), plot = H)
