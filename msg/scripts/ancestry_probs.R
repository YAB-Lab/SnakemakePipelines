#!/usr/bin/env Rscript
#
# Aggregate per-sample HMM fits into genome-wide ancestry probability tables
# and summary plots.
#
# Replaces msg/summaryPlots.R + scripts/run_summary_plots.sh wrapper.
# Differences from the legacy script:
#   - optparse CLI; no shell wrapper.
#   - Sources hmmlib.R + ded.R from this script's own directory.
#   - All outputs land under explicit --ancestry-dir / --images-dir paths
#     (legacy script wrote to CWD and the rule mv'd files post-hoc).
#   - Cache .rda files (interpolate.probs memoization) live in
#     <images-dir>/cache/, not in CWD.
#   - The dead `plot.correlation.matrix` and `breakpoint.widths` branches
#     are now config-gated via CLI flags (default off, matching legacy).
#
# Output (under --ancestry-dir):
#   ancestry-probs-par1.tsv, ancestry-probs-par2.tsv,
#   ancestry-probs-par1par2.tsv
# Output (under --images-dir):
#   offdiagonal_data.tsv, rhat-offdiag.rda, rlod-offdiag.rda,
#   offdiagonal-lod.pdf, segregation.pdf, missing.pdf

options(error = quote(q("yes")))

suppressPackageStartupMessages(library(optparse))

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

opt_list <- list(
    make_option(c("-d", "--hmm-fit-dir"), type = "character",
                help = "Directory containing per-sample <indiv>/<indiv>-hmmprob.RData"),
    make_option(c("--ancestry-dir"), type = "character",
                help = "Directory to write the 3 ancestry-probs TSVs into"),
    make_option(c("--images-dir"), type = "character",
                help = "Directory to write PDFs + offdiag rda/tsv into"),
    make_option(c("-l", "--chr-lengths"), type = "character",
                help = "CSV file with chr,length header"),
    make_option(c("-c", "--chroms"), type = "character", default = "all"),
    make_option(c("-p", "--chroms2plot"), type = "character", default = "all"),
    make_option(c("-t", "--thinfac"), type = "double", default = 1),
    make_option(c("-f", "--diffac"), type = "double", default = 0.01),
    make_option(c("-n", "--pna-thresh"), type = "double", default = 0.03),
    make_option(c("--plot-correlation-matrix"), type = "integer", default = 0,
                help = "0/1; gates heatmap + thinning branch (heavy)"),
    make_option(c("--breakpoint-widths"), type = "integer", default = 0,
                help = "0/1; sources scripts/breakpoint-widths.R when 1")
)
opts <- parse_args(OptionParser(option_list = opt_list))

stopifnot(
    !is.null(opts$`hmm-fit-dir`),
    !is.null(opts$`ancestry-dir`),
    !is.null(opts$`images-dir`),
    !is.null(opts$`chr-lengths`)
)

dir          <- opts$`hmm-fit-dir`
ancestry_dir <- opts$`ancestry-dir`
images_dir   <- opts$`images-dir`
chrLen       <- opts$`chr-lengths`
thinfac      <- opts$thinfac
difffac      <- opts$diffac
pna.thresh   <- opts$`pna-thresh`
plot.correlation.matrix <- as.logical(opts$`plot-correlation-matrix`)
breakpoint.widths       <- as.logical(opts$`breakpoint-widths`)
write.ancestry.probs    <- TRUE

# ---------------------------------------------------------------------------
# Source HMM math libs (next to this script)
# ---------------------------------------------------------------------------

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "",
                   script_args[grep("^--file=", script_args)])
script_dir  <- dirname(normalizePath(script_path))

source(file.path(script_dir, "ded.R"))
source(file.path(script_dir, "hmmlib.R"))

# ---------------------------------------------------------------------------
# Output dirs + interpolate.probs cache
# ---------------------------------------------------------------------------

dir.create(ancestry_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(images_dir,   recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(images_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Contig setup
# ---------------------------------------------------------------------------

contigs      <- as.character(strsplit(opts$chroms, ",")[[1]])
contigs2plot <- as.character(strsplit(opts$chroms2plot, ",")[[1]])

chrLengths <- read.csv(chrLen, row.names = 1)
chrLengths <- structure(chrLengths$length,
                        names = as.character(rownames(chrLengths)))

if (opts$chroms == "all") {
    contigs <- as.character(names(chrLengths))
} else {
    contigs <- names(chrLengths)[match(contigs, names(chrLengths),
                                       nomatch = "0")]
    contigs <- contigs[chrLengths[contigs] > 0]
}

if (opts$chroms2plot == "all") {
    contigs2plot <- as.character(names(chrLengths))
}

# ---------------------------------------------------------------------------
# Cached interpolate.probs wrapper
# ---------------------------------------------------------------------------

get.ancestry.probs <- function(ancestry, thinfac, difffac, contigs2use,
                               pna.thresh = 1, type = "all") {
    fname <- file.path(cache_dir,
        sprintf("ancestry-probs-%s-%f-%f-%f.%s.rda",
                ancestry, thinfac, difffac, pna.thresh, type))
    if (!file.exists(fname)) {
        pp.all <- lapply(contigs2use, interpolate.probs, dir = dir,
                         thinfac = thinfac, difffac = difffac,
                         ancestry = ancestry)
        names(pp.all) <- contigs2use
        save(pp.all, file = fname)
    } else {
        cat(type, " ", ancestry, ": ", fname, " file exists\n", sep = "")
        pp.all <- read.object(fname)

        # If cached contigs don't match the request, rebuild
        if (sum(contigs2use %in% names(pp.all)) != length(contigs2use)) {
            pp.all <- lapply(contigs2use, interpolate.probs, dir = dir,
                             thinfac = thinfac, difffac = difffac,
                             ancestry = ancestry)
            names(pp.all) <- contigs2use
            save(pp.all, file = fname)
        }
    }
    pp.pna <- lapply(lapply(pp.all, is.na), colMeans)

    cat(type, " ", ancestry,
        ": Numbers of markers after removal for high missing data proportion (%NA<",
        pna.thresh, "):\n", sep = "")
    print(sapply(pp.pna, function(pna)
        paste(sum(pna < pna.thresh), "/", length(pna), sep = "")))
    pp <- mapply(function(p, pp.pna) p[, pp.pna < pna.thresh, drop = FALSE],
                 pp.all, pp.pna, SIMPLIFY = FALSE)
    names(pp) <- contigs2use
    pp
}

# ---------------------------------------------------------------------------
# Full marker set (no thinning) -> 3 ancestry TSVs
# ---------------------------------------------------------------------------

cat("Extracting markers for the following contigs:\n", opts$chroms)
pp1 <- get.ancestry.probs("par1", thinfac = 1, difffac = 0,
                          contigs2use = contigs, type = "all")
pp  <- pp2 <- get.ancestry.probs("par2", thinfac = 1, difffac = 0,
                                 contigs2use = contigs, type = "all")
p1 <- do.call("cbind", pp1); rm(pp1)
p  <- p2 <- do.call("cbind", pp2); rm(pp2)

contig.lengths <- sapply(pp, ncol); rm(pp)
contig.fac <- factor(rep(contigs, contig.lengths))
current_contigs <- levels(contig.fac)[match(contigs, levels(contig.fac),
                                            nomatch = "0")]

pos <- as.integer(colnames(p)); rm(p)
info <- contig.info(pos, contig.fac, chrLengths[current_contigs])

if (write.ancestry.probs) {
    p1.table <- p1; rm(p1)
    p2.table <- p2; rm(p2)
    colnames(p1.table) <- colnames(p2.table) <-
        paste(contig.fac, pos, sep = ":")
    p12.table <- 1 - (p1.table + p2.table)
    msg.write.table(round(p1.table, 6),
                    file = file.path(ancestry_dir, "ancestry-probs-par1.tsv"))
    msg.write.table(round(p2.table, 6),
                    file = file.path(ancestry_dir, "ancestry-probs-par2.tsv"))
    msg.write.table(round(p12.table, 6),
                    file = file.path(ancestry_dir,
                                     "ancestry-probs-par1par2.tsv"))
    rm(p1.table, p2.table, p12.table)
} else {
    rm(p1, p2)
}

# ---------------------------------------------------------------------------
# Off-diagonal LOD profile + summary plots
# ---------------------------------------------------------------------------

genomeplot <- function(y, ...) {
    plot(x = info$genomepos, y = y, type = "l", bty = "n",
         xaxt = "n", xlab = "", ...)
    abline(v = info$boundaries, lty = 2, col = "blue")
    axis(side = 1, at = info$boundaries, labels = FALSE, tick = TRUE)
    axis(side = 1, at = info$midpoints,
         labels = names(info$midpoints), tick = FALSE, las = 2)
}

cat("\nOff-diagonal LOD profile...\n")
pp <- pp2 <- get.ancestry.probs("par2", thinfac = 1, difffac = 0,
                                contigs2use = contigs2plot, type = "plot")
p <- p2 <- do.call("cbind", pp2); rm(p2, pp2)

contig.lengths <- sapply(pp, ncol); rm(pp)
contig.fac <- factor(rep(contigs2plot, contig.lengths))
current_contigs <- levels(contig.fac)[match(contigs2plot, levels(contig.fac),
                                            nomatch = "0")]

pos  <- as.integer(colnames(p))
info <- contig.info(pos, contig.fac, chrLengths[current_contigs])

rhat.offdiag <- est.rf.p.profile(p = p, offdiag = TRUE, lod = FALSE)
save(rhat.offdiag, file = file.path(images_dir, "rhat-offdiag.rda"))

rlod.offdiag <- est.rf.p.profile(p = p, offdiag = TRUE, lod = TRUE)
save(rlod.offdiag, file = file.path(images_dir, "rlod-offdiag.rda"))

write.table(cbind(as.vector(contig.fac), pos, rhat.offdiag, rlod.offdiag),
            file = file.path(images_dir, "offdiagonal_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = c("chrom", "pos", "rhat", "rlod"))

pdf(file.path(images_dir, "offdiagonal-lod.pdf"), width = 10, height = 5)
genomeplot(rlod.offdiag, ylab = "LOD")
abline(h = 25, col = "red")
dev.off()

pdf(file.path(images_dir, "segregation.pdf"), width = 10, height = 5)
genomeplot(colMeans(p, na.rm = TRUE), ylab = "Mean probability",
           main = "Average probability of par2 homozygosity (hemizygosity)")
abline(h = 1 / 2, col = "red", lty = 2)
legend("bottomleft", legend = "Expectation", col = "red",
       lty = 2, bty = "n")
dev.off()

pdf(file.path(images_dir, "missing.pdf"), width = 10, height = 5)
genomeplot(colMeans(is.na(p)), ylab = "Missing proportion")
dev.off()

rm(p, rhat.offdiag, rlod.offdiag)

# ---------------------------------------------------------------------------
# Optional: thinned ancestry probs + correlation heatmap
# ---------------------------------------------------------------------------

if (plot.correlation.matrix) {
    cat("\nThinning ancestry probabilities...\n")

    pp.thin <- pp2.thin <- get.ancestry.probs(
        "par2", thinfac = thinfac, difffac = difffac,
        contigs2use = contigs2plot, pna.thresh = pna.thresh,
        type = "thinned_plot")
    p.thin <- p2.thin <- do.call("cbind", pp2.thin); rm(pp2.thin)
    contig.lengths.thin <- sapply(pp.thin, ncol)
    contig.fac.thin <- factor(rep(contigs2plot, contig.lengths.thin))
    current_contigs <- levels(contig.fac.thin)[
        match(contigs2plot, levels(contig.fac.thin), nomatch = "0")]

    if (length(current_contigs) > 0) {
        pos.thin  <- as.integer(colnames(p.thin))
        info.thin <- contig.info(pos.thin, contig.fac.thin,
                                 chrLengths[current_contigs])

        if (write.ancestry.probs) {
            p2.thin.table <- p2.thin
            colnames(p2.thin.table) <- paste(contig.fac.thin, pos.thin,
                                             sep = ":")
            msg.write.table(round(p2.thin.table, 6),
                file = file.path(ancestry_dir,
                                 "ancestry-probs-thinned-par2.tsv"))
            rm(p2.thin.table)
        }

        n <- nrow(pp.thin[[1]]); rm(pp.thin)
        lod.max <- n * log10(2)
        lod.thin.fname <- file.path(cache_dir,
            sprintf("lod-thin-%f-%f-%f.rda", thinfac, difffac, pna.thresh))

        if (!file.exists(lod.thin.fname)) {
            lod.thin <- est.rf.p(p.thin, lod = TRUE, na.rm = TRUE)
            save(lod.thin, file = lod.thin.fname)
        } else {
            lod.thin <- read.object(lod.thin.fname)
        }
        rm(p.thin)

        cat("Plotting heatmap...\n")
        bitmap(file = file.path(images_dir, "lod-matrix.bmp"),
               width = 100, height = 100, bg = "transparent")
        plot.rf.ded(pos = info.thin$genomepos, lod.thin,
                    zmax = lod.max, info = info.thin)
        dev.off()
    } else {
        cat("Unable to create thinned lod-matrix.bmp plot: ",
            "Try increasing pna.thresh (", pna.thresh,
            ") or decreasing difffac (", difffac, ")\n", sep = "")
    }
}

if (breakpoint.widths) {
    figdir <- images_dir   # breakpoint-widths.R writes pdfs to `figdir`
    source(file.path(script_dir, "breakpoint-widths.R"))
}
