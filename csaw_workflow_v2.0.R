#!/usr/bin/env Rscript
# csaw_multiTreat_vs_control_cleaned.R
#
# 1) Merge → blacklist filter → recount
# 2) Single outdir for all output
# 3) Tab files sorted by chrom/pos
# 4) Collect significant windows (FDR<0.05) across contrasts
# 5) MA‐plots with unified axes (＾◇＾)
#

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(csaw)
  library(edgeR)
  library(ggplot2)
  library(ggrastr)
  library(Rsamtools)
  library(BiocParallel)
  library(rtracklayer)
  library(tools)
  library(limma)
})

# 1) Command‐line options -----------------------------------------------------
opt_list <- list(
  make_option(c("-o","--outdir"),  type="character", default="csaw_out"),
  make_option(c("-g","--genome"),  type="character", default="mouse"),
  make_option("--control",         type="character",               help="Control BAMs, comma‐sep"),
  make_option("--treat1",          type="character",               help="Treat1 BAMs, comma‐sep"),
  make_option("--treat2",          type="character", default=NULL),
  make_option("--treat3",          type="character", default=NULL),
  make_option("--treat4",          type="character", default=NULL),
  make_option("--treat5",          type="character", default=NULL),
  make_option("--treat6",          type="character", default=NULL),
  make_option("--spike",           type="character", default=NULL, help="Spike‐in genome: ecoli or drosophila"),
  make_option("--workers",         type="integer",   default=1),
  make_option("--batches",         type="character", default=NULL, help="Comma-sep batch labels, one per sample in the SAME order as BAMs")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# 2) Gather & validate BAM groups ---------------------------------------------
t_args <- Filter(Negate(is.null),
                 list(opt$treat1,opt$treat2,opt$treat3,
                      opt$treat4,opt$treat5,opt$treat6))
nTreat <- length(t_args)
stopifnot(nTreat >= 1 && nTreat <= 6)
treats <- lapply(t_args, function(x) strsplit(x, ",")[[1]])
stopifnot(!any(lengths(treats) < 2))
ctrl   <- strsplit(opt$control, ",")[[1]]
stopifnot(length(ctrl) >= 2)

all_bams  <- unlist(c(treats, list(ctrl)))
group_lvl <- c(paste0("treat",1:nTreat), "control")

# --- parse optional batches ---------------------------------
if (!is.null(opt$batches)) {
  batch_vec <- strsplit(opt$batches, ",")[[1]]
  if (length(batch_vec) != length(all_bams)) {
    stop("Must supply one batch label per sample (", length(all_bams), " expected)")
  }
  batch <- factor(batch_vec)
  message(" → Batches defined: ", paste(levels(batch), collapse=", "))
} else {
  batch <- NULL
}

message("⋌(∘⁼̴⃙̀˘⌂˘⁼̴⃙́∘)⋋  Groups:")
for(i in seq_len(nTreat))
  message("   ", group_lvl[i], ": ", length(treats[[i]]), " BAMs")
message("   control: ", length(ctrl), " BAMs")

# 3) Chromosomes & parallel ---------------------------------------------------
std_chrs <- switch(opt$genome,
  mouse = paste0("chr", c(1:19,"X","Y")),
  human = paste0("chr", c(1:22,"X","Y")),
  stop("Unknown genome: ", opt$genome)
)
bpparam  <- MulticoreParam(workers = opt$workers)

# 4) Prefix helper ------------------------------------------------------------
prefix_of <- function(bam) {
  sub("_rep[0-9]+$","",
      file_path_sans_ext(basename(bam)))
}
prefixes <- vapply(treats, function(x) prefix_of(x[1]), "")

# 5) Read & union per‐group merged peaks ------------------------------------
message("(づ ・ ・)づ reading per‐rep peaks…")
read_peaks <- function(bam) {
  base <- file_path_sans_ext(basename(bam))
  bedf <- file.path(dirname(dirname(dirname(bam))),
                    opt$genome, "MACS2_out", "individual_reps",
                    paste0(base, "_peaks.narrowPeak"))
  df   <- read.table(bedf, sep="\t", header=FALSE)[,1:3]
  gr   <- GRanges(seqnames = df[[1]],
                  ranges   = IRanges(df[[2]], df[[3]]))
  keepSeqlevels(gr, std_chrs, pruning.mode="coarse")
}
plist     <- c(lapply(unlist(treats), read_peaks),
               lapply(ctrl,           read_peaks))
all_peaks <- Reduce(union, plist)
message(" → Union windows: ", length(all_peaks))

# 6) Initial count & filter --------------------------------------------------
param <- readParam(max.frag=1000, pe="both", restrict=std_chrs)
rc0   <- regionCounts(all_bams, all_peaks,
                      param = param, BPPARAM = bpparam)
message(" → regionCounts: ", nrow(rc0), "×", ncol(rc0))

y0   <- asDGEList(rc0)
grp0 <- rep(group_lvl, lengths(c(treats, list(ctrl))))
keep <- filterByExpr(y0, group = grp0)
message(" → filterByExpr keeps ", sum(keep), "/", length(keep))
rcf  <- rc0[keep,]

# 7) Merge & blacklist filter -----------------------------------------------
message("乁(  ՞ ▽ ՞)ㄏ merging windows…")
merged <- mergeWindows(rowRanges(rcf),
                       tol=50L, max.width=5000L)
regions <- merged$region
message(" → merged into ", length(regions), " regions")

# blacklist
blist_file <- switch(opt$genome,
  mouse = "/project/GCRB/shared/bowtie2_indexes/mm10/genome_files/mm10_Boyle_blacklist.bed",
  human = "/project/GCRB/shared/bowtie2_indexes/hg38/genome_files/hg38_Boyle_blacklist.bed",
  stop("No blacklist for genome ",opt$genome)
)
bl   <- import(blist_file)
orig <- length(regions)
regions <- regions[!countOverlaps(regions, bl)]
message(" → removed ", orig - length(regions),
        " blacklisted ( ☆ω☆); left ", length(regions))

# 8) Re‐count on merged regions ----------------------------------------------
message("(۶ ＾▽＾)۶ recounting on merged, filtered regions…")
rc2 <- regionCounts(all_bams, regions,
                    param   = param,
                    BPPARAM = bpparam)
message(" → recount: ", nrow(rc2), "×", ncol(rc2))

# 9) Normalize & fit ---------------------------------------------------------
if (!is.null(opt$spike)) {
  spike_bams <- gsub(
    paste0("/",opt$genome,"/bam/"),
    paste0("/",opt$spike,"_spikein/bam/"),
    all_bams
  )
  sc <- unlist(bplapply(spike_bams,
         function(b) countBam(BamFile(b))$records,
         BPPARAM=bpparam))
  sf <- sc / median(sc)
  y2 <- asDGEList(rc2)
  y2$samples$norm.factors <- sf
  message(" → Spike‐in norm factors: ", paste(round(sf,3), collapse=","))
} else {
  y2 <- normOffsets(rc2, se.out=TRUE) |> asDGEList()
  message(" → TMM normalization applied")
}

y2$samples$group <- factor(grp0, levels = c("control", paste0("treat",1:nTreat)))

# Build design with batch if requested
if (!is.null(batch)) {
  y2$samples$batch <- batch
  design <- model.matrix(~ batch + group, data = y2$samples)
  message(" → Design: ~batch + group")
} else {
  design <- model.matrix(~0 + group, data = y2$samples)
  colnames(design) <- levels(y2$samples$group)
  message(" → Design: ~0 + group")
}

colnames(design)

y2              <- estimateDisp(y2, design)
fit             <- glmQLFit(y2, design, robust=TRUE)

# 10) Define outdir & set up containers -------------------------------------
outdir <- file.path(opt$outdir, opt$genome, "csaw_out")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
results <- vector("list", nTreat)
names(results) <- paste0("treat",1:nTreat)

# 11) Loop contrasts & collect results ---------------------------------------
message(" ⋌( ∘⁼̴⃙̀˘ ⌂ ˘⁼̴⃙́∘ )⋋ running contrasts…")
for (i in seq_len(nTreat)) {
  nm  <- paste0("treat",i)
  # figure out which coefficient corresponds to treat_i
  if ("(Intercept)" %in% colnames(design)) {
    # intercept-based design: treatment effect is "group" + nm
    coef_name <- paste0("group", nm)  # e.g. "grouptreat1"
    if (!coef_name %in% colnames(design)) {
      stop("Cannot find coefficient '", coef_name, "' in design columns: ",
           paste(colnames(design), collapse=", "))
    }
    ct <- makeContrasts(contrasts = coef_name, levels = design)

  } else {
    # 0+group design: explicit control column
    contrast_string <- paste0(nm, " - control")
    ct <- makeContrasts(contrasts = contrast_string, levels = design)
  }

  # run the test
  tbl <- glmQLFTest(fit, contrast = ct)$table
  if (!"FDR" %in% colnames(tbl)) {
    tbl$FDR <- p.adjust(tbl$PValue, method="BH")
  }
  stopifnot(nrow(tbl) == length(regions))
  gr <- regions
  mcols(gr) <- DataFrame(tbl)
  results[[i]] <- gr
}
message(" ✓ all contrasts done!")

# 12) Determine MA‐plot scales & significant windows -------------------------
# Combine all logCPM/logFC to get global axis limits
all_lc <- unlist(lapply(results, function(gr) mcols(gr)$logCPM))
all_lf <- unlist(lapply(results, function(gr) mcols(gr)$logFC))
xlim  <- range(all_lc, na.rm=TRUE)
ylim  <- range(all_lf, na.rm=TRUE)

# Global volcano limits (same for all plots)
all_fdr  <- unlist(lapply(results, function(gr) mcols(gr)$FDR))
vol_xmax <- ceiling(max(abs(all_lf), na.rm=TRUE))
vol_ymax <- ceiling(max(-log10(all_fdr), na.rm=TRUE))
vol_ymax <- pmin(pmax(vol_ymax, 0), 350)
vol_ymax
# Collect *any* region with FDR<0.05
# 1) Pull out significant GRanges from each contrast
sig_list <- lapply(results, function(gr) {
  gr[mcols(gr)$FDR < 0.05]
})

# 2) Unlist into one big GRanges
sig_grs <- unlist(GRangesList(sig_list), use.names = FALSE)
message("Collected ", length(sig_grs), " significant windows (pre‐merge)")

# 3) Merge any overlapping windows
sig_grs <- reduce(sig_grs)
message("After reduce(): ", length(sig_grs), " merged significant windows")

# 4) Clean up seqlevels and sort by genomic coordinates
sig_grs <- sortSeqlevels(sig_grs)
sig_grs <- sort(sig_grs)
message("After sortSeqlevels & sort(): first windows:")
print(head(sig_grs))

# Write bed
bed_file <- file.path(outdir, "csaw_significant_windows.bed")
message(" (づ ＾▽＾)づ writing significant windows bed…")
dfb <- as.data.frame(sig_grs)[,c("seqnames","start","end")]
write.table(dfb, bed_file,
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

#12.5 === Export normalized counts per region =======================================
message("(づ  ^ ヮ ^)づ✎ exporting normalized counts...")

# Get normalized rawcpm (you can change log=TRUE if you want log2 CPMs)
rawcpm <- edgeR::cpm(y2, log=FALSE, normalized.lib.sizes=TRUE)

# Rename samples using BAM basenames
bam_names <- basename(tools::file_path_sans_ext(all_bams))
stopifnot(length(bam_names) == ncol(rawcpm))
colnames(rawcpm) <- bam_names

# Extract groupings from the design matrix (assumes sample names are unique)
sample_groups <- y2$samples$group
stopifnot(length(sample_groups) == ncol(rawcpm))

# Compute per-group mean rawcpm
group_levels <- unique(sample_groups)
mean_cpm <- sapply(group_levels, function(g) {
  rows <- which(sample_groups == g)
  rowMeans(rawcpm[, rows, drop=FALSE])
})
colnames(mean_cpm) <- paste0("mean_", group_levels)

# Combine into one table
tab_out <- data.frame(
  as.data.frame(regions)[,1:3],  # seqnames/start/end
  rawcpm,                        # all replicates
  mean_cpm                       # mean per group
)

# Clean up column names
colnames(tab_out)[1:3] <- c("chr", "start", "end")

# Write to file
write.table(tab_out, file=file.path(outdir,"csaw_normalized_counts.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

message(" ✓ wrote normalized counts to csaw_normalized_counts.tsv")

# 13) Write per‐contrast .tab & MA plots --------------------------------------
message("(  ^ ヮ ^)و writing outputs…")
for (i in seq_len(nTreat)) {
  prefix <- prefixes[i]
  gr     <- results[[i]]

  # sort by chrom then start
  ord  <- order(seqnames(gr), start(gr))
  sgr  <- gr[ord]
  tab  <- as.data.frame(sgr)
  tab_file <- file.path(outdir, paste0(prefix, "_vs_control.tab"))
  write.table(tab, tab_file,
              sep="\t", quote=FALSE, row.names=FALSE)

  # Filter significant regions
  sig_gr <- sgr[sgr$FDR < 0.05]

  # If any found, write BED-like file
  if (length(sig_gr) > 0) {
    sig_df <- data.frame(
      seqnames = as.character(seqnames(sig_gr)),
      start = start(sig_gr) - 1,  # BED format: 0-based start
      end = end(sig_gr),          # BED format: 1-based end
      log2FC = signif(sig_gr$logFC, 3),
      FDR = signif(sig_gr$FDR, 3)
    )
    bed_file <- file.path(outdir, paste0(prefix, "_vs_control_significant.bed"))
    write.table(
      sig_df,
      bed_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    message("  ✓ Wrote BED of significant regions: ", bed_file)
  } else {
    message("  ⚠ No significant regions (FDR < 0.05) to write BED file.")
  }

  # MA plot
    df <- as.data.frame(sgr[order(sgr$FDR),])
    df$negLogFDR <- -log10(pmin(pmax(df$FDR, 1e-350), 1))  # clamp to avoid -Inf

    p <- ggplot(df, aes(x = logCPM, y = logFC)) +
      ggrastr::geom_point_rast(
        data = df[df$FDR >= 0.05,],  # nonsig grey
        color = "grey70", size = 2, alpha = 0.6, raster.dpi = 300
      ) +
      ggrastr::geom_point_rast(
        data = df[df$FDR < 0.05,],
        aes(color = negLogFDR),
        size = 2, alpha = 0.6, raster.dpi = 300
      ) +
      scale_color_gradientn(
        colors = c("#3D0C3E", "#BA324F", "#FABF4B"),
        name = expression(-log[10]~FDR)
      ) +
      geom_smooth(method = "loess", se = FALSE, color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(
        title = paste0(prefix, " vs control"),
        x = expression(log[2]~CPM), y = expression(log[2]~FoldChange)
      ) +
      xlim(xlim) + ylim(ylim) +
      theme_bw(base_size = 14) +
      theme(
        legend.position = "right",
        aspect.ratio = 1,
        axis.text = element_text(color = "black")
      )

    ggsave(file.path(outdir, paste0(prefix, "_vs_control_MA.pdf")),
           plot = p, width = 8, height = 8)
  message("  ✓ wrote ", prefix)

  # Volcano plot ---------------------------------------------
    volc_df <- as.data.frame(sgr)
    volc_df$negLogFDR <- -log10(pmin(pmax(volc_df$FDR, 1e-350), 1))

    volc_df$sig <- volc_df$FDR < 0.05
    xmax <- ceiling(max(abs(volc_df$logFC), na.rm = TRUE))

    p_volc <- ggplot(volc_df, aes(x = logFC, y = negLogFDR)) +
      ggrastr::geom_point_rast(
        data = volc_df[!volc_df$sig,],
        color = "grey70", size = 2, alpha = 0.6, raster.dpi = 300
      ) +
      ggrastr::geom_point_rast(
        data = volc_df[volc_df$sig,],
        aes(color = logCPM),
        size = 2, alpha = 0.6, raster.dpi = 300
      ) +
      scale_color_gradientn(
        colors = c("#3D0C3E", "#BA324F", "#FABF4B"),
        name = expression(log[2]~CPM)
      ) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
      coord_cartesian(xlim = c(-vol_xmax, vol_xmax), ylim = c(0, vol_ymax)) +
      labs(
        title = paste0(prefix, " volcano"),
        x = expression(log[2]~FoldChange), y = expression(-log[10]~FDR)
      ) +
      theme_bw(base_size = 14) +
      theme(
        legend.position = "right",
        aspect.ratio = 1,
        axis.text = element_text(color = "black")
      )

    ggsave(file.path(outdir, paste0(prefix, "_vs_control_volcano.pdf")),
           plot = p_volc, width = 8, height = 8)
   message("  ✓ volcano plotted for ", prefix)
}

message("All done! ‎(و̵ ˃̵ ヮ ˂̵)۶ Clean, black‐list‐filtered, consolidated csaw results.")
