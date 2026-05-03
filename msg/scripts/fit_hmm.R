#!/usr/bin/env Rscript
#
# Fit per-contig HMMs and emit ancestry-probability plots + breakpoints.
#
# Replaces the legacy msg/fit-hmm.R + scripts/run_fit_hmm.sh wrapper.
# Differences from the legacy script:
#   - optparse CLI (one flag per fit_hmm: config key); no shell wrapper.
#   - Sources hmmlib.R + ded.R from this script's own directory; the
#     hmmprobs C binary is also expected alongside (no msg/ dependency).
#   - Drops the for(indiv in indivs) outer loop. Snakemake parallelizes
#     by sample, so each invocation handles one individual.
#   - Drops library(R.methodsS3) / library(R.oo); charToInt() replaced
#     by base utf8ToInt().
#
# Output (per sample, into --outdir/<indiv>/):
#   <indiv>-hmmprob.pdf, <indiv>-hmmprob.RData,
#   <indiv>-breakpoints.csv, <indiv>-matchMismatch.csv,
#   plus per-contig <indiv>-<contig>-breakpoints.gff.

options(error = quote(q("yes")))

suppressPackageStartupMessages(library(optparse))

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

opt_list <- list(
    make_option(c("-i", "--indiv"), type = "character",
                help = "Sample / individual ID"),
    make_option(c("-d", "--hmm-data-dir"), type = "character",
                help = "Directory holding <indiv>/<indiv>-<contig>.hmmdata"),
    make_option(c("-o", "--outdir"), type = "character",
                help = "Directory to write <indiv>/ HMM fit outputs into"),
    make_option(c("-l", "--chr-lengths"), type = "character",
                help = "CSV file with chr,length header"),
    make_option(c("-s", "--sex"), type = "character",
                help = "'male' or 'female'"),
    make_option(c("--deltapar1"), type = "double", default = 0.01),
    make_option(c("--deltapar2"), type = "double", default = NA,
                help = "Defaults to deltapar1 if unset"),
    make_option(c("--rec-rate"), type = "double", default = 0,
                help = "0 -> per-contig 1/length; >0 -> rate/total_length"),
    make_option(c("--rfac"), type = "double", default = 1e-6),
    make_option(c("--chroms"), type = "character", default = "all",
                help = "Comma-separated contigs to fit, or 'all'"),
    make_option(c("--sex-chroms"), type = "character", default = "X",
                help = "Comma-separated sex contigs, or 'all'"),
    make_option(c("--priors"), type = "character", default = "0.25,0.5,0.25",
                help = "Comma-separated prior on (par1/par1, par1/par2, par2/par2)"),
    make_option(c("--theta"), type = "double", default = 1,
                help = "PHRED quality value correction"),
    make_option(c("--gff-thresh-conf"), type = "double", default = 0.95),
    make_option(c("--one-site-per-contig"), type = "integer", default = 1,
                help = "0/1; subsamples one site per read when 1"),
    make_option(c("--chroms2plot"), type = "character", default = "all"),
    make_option(c("--pepthresh"), type = "character", default = "null",
                help = "If not 'null', also write per-contig CSVs"),
    make_option(c("--use-filter-hmmdata"), type = "integer", default = 0,
                help = "0/1; read pre-filtered hmmdata if 1")
)
opts <- parse_args(OptionParser(option_list = opt_list))

stopifnot(
    !is.null(opts$indiv),
    !is.null(opts$`hmm-data-dir`),
    !is.null(opts$outdir),
    !is.null(opts$`chr-lengths`),
    !is.null(opts$sex)
)

indiv  <- opts$indiv
dir    <- opts$`hmm-data-dir`
outdir <- opts$outdir
sex    <- opts$sex
chrLen <- opts$`chr-lengths`

deltapar1 <- opts$deltapar1
deltapar2 <- if (is.na(opts$deltapar2)) deltapar1 else opts$deltapar2
recrate   <- opts$`rec-rate`
rfac      <- opts$rfac
priors    <- strsplit(opts$priors, ",")[[1]]
theta     <- opts$theta
gff_thresh_conf     <- opts$`gff-thresh-conf`
one.site.per.read   <- as.logical(opts$`one-site-per-contig`)
use_filtered_hmmdata <- as.logical(opts$`use-filter-hmmdata`)
pepthresh <- opts$pepthresh

# ---------------------------------------------------------------------------
# Source HMM math libs (next to this script)
# ---------------------------------------------------------------------------

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "",
                   script_args[grep("^--file=", script_args)])
script_dir  <- dirname(normalizePath(script_path))

source(file.path(script_dir, "ded.R"))
source(file.path(script_dir, "hmmlib.R"))

# Pr.y.given.z(C=TRUE) shells out to `<dir>/hmmprobs`; point it at scripts/.
hmmprobs_dir <- script_dir

dir.create(file.path(outdir, indiv), recursive = TRUE, showWarnings = FALSE)

cat("one.site.per.read has been set to", one.site.per.read, "\n")
if (use_filtered_hmmdata) cat("Using pre-filtered hmmdata file\n")

# ---------------------------------------------------------------------------
# Contig setup
# ---------------------------------------------------------------------------

minCoverage <- 0

contigLengths <- read.csv(chrLen, header = TRUE, sep = ",", as.is = TRUE)
contigs <- as.character(sort(as.vector(contigLengths$chr)))
rownames(contigLengths) <- contigLengths$chr

main.contigs <- as.character(strsplit(opts$chroms, ",")[[1]])
plot.contigs <- as.character(strsplit(opts$chroms2plot, ",")[[1]])
sex.chroms   <- as.character(strsplit(opts$`sex-chroms`, ",")[[1]])
if (opts$chroms     == "all") main.contigs <- contigs
if (opts$chroms2plot == "all") plot.contigs <- contigs
if (opts$`sex-chroms` == "all") sex.chroms <- contigs  # haplodiploid

aveSpace <- sum(as.numeric(
    contigLengths[contigLengths$chr %in% plot.contigs, ]$length)) /
    length(plot.contigs)
plotPadding <- 10 ^ (ceiling(log10(aveSpace)) - 2)

alleles <- c("A", "C", "G", "T")

# ---------------------------------------------------------------------------
# Geneious GFF writer (ported verbatim from legacy fit-hmm.R)
# ---------------------------------------------------------------------------

output_geneious_file <- function(gff_thresh_conf, x, y, par1homo_col,
                                 par2homo_col, indiv, contig, contigLengths) {
    gff_thresh_inverse <- 1 - gff_thresh_conf
    gff_data <- breakpoint.width(x, y[, par1homo_col], y[, par2homo_col],
                                 indiv = indiv, contig = contig,
                                 conf1 = gff_thresh_inverse,
                                 conf2 = gff_thresh_conf)

    if (!is.null(gff_data[["bps"]])) {
        gff_for_output <- {}
        start_inner <- start_outer <- end_inner <- end_outer <- NULL

        if (y[1, par1homo_col] > .5) {
            start_inner <- 1
            start_outer <- 1
        }

        for (i in seq_len(nrow(gff_data[["bps"]]))) {
            row <- gff_data[["bps"]][i, ]
            if (!is.null(start_inner) && !is.null(start_outer)) {
                stopifnot(is.null(end_inner), is.null(end_outer))
                end_inner <- as.numeric(as.character(row[1, 6]))
                end_outer <- as.numeric(as.character(row[1, 7]))
            } else {
                stopifnot(is.null(start_inner), is.null(start_outer))
                start_inner <- as.numeric(as.character(row[1, 7]))
                start_outer <- as.numeric(as.character(row[1, 6]))
            }
            if (!is.null(start_inner) && !is.null(end_inner)) {
                gff_for_output <- rbind(gff_for_output,
                    c("", "Geneious", "msg_run", start_inner, end_inner,
                      ".", ".", ".", paste("name", indiv, sep = "=")))
                start_inner <- NULL; end_inner <- NULL
            }
            if (!is.null(start_outer) && !is.null(end_outer)) {
                gff_for_output <- rbind(gff_for_output,
                    c("", "Geneious", "msg_run", start_outer, end_outer,
                      ".", ".", ".", paste("name", indiv, sep = "=")))
                start_outer <- NULL; end_outer <- NULL
            }
        }
        if (!is.null(start_inner) && !is.null(start_outer)) {
            stopifnot(is.null(end_inner), is.null(end_outer))
            end_inner <- contigLengths[contigLengths$chr == contig, "length"]
            end_outer <- contigLengths[contigLengths$chr == contig, "length"]
            gff_for_output <- rbind(gff_for_output,
                c("", "Geneious", "msg_run", start_inner, end_inner,
                  ".", ".", ".", paste("name", indiv, sep = "=")))
            gff_for_output <- rbind(gff_for_output,
                c("", "Geneious", "msg_run", start_outer, end_outer,
                  ".", ".", ".", paste("name", indiv, sep = "=")))
        }
        write.table(gff_for_output,
            file = file.path(outdir, indiv,
                paste(indiv, contig, "breakpoints.gff", sep = "-")),
            append = FALSE, quote = FALSE, na = "NA",
            row.names = FALSE, col.names = FALSE, sep = "\t")
    }
}

# ---------------------------------------------------------------------------
# Fit / load per-contig HMMs
# ---------------------------------------------------------------------------

cat(indiv, "\n")

dataa <- list()
hmmdata.file <- file.path(outdir, indiv, paste(indiv, "hmmprob.RData", sep = "-"))

if (file.exists(hmmdata.file)) {
    cat("HMM fit for indiv", indiv, "already exists\n")
    dataa <- read.object(hmmdata.file)
    if (is.null(names(dataa)) && length(dataa) == length(contigs)) {
        names(dataa) <- as.character(contigs)
    }
} else {
    for (contig in main.contigs) {

        if (sex == "male" && contig %in% sex.chroms) {
            ploidy <- 1
            ancestries <- c("par1", "par2")
            phi <- rep(1 / length(ancestries), length(ancestries))
        } else {
            ploidy <- 2
            ancestries <- c("par1/par1", "par1/par2", "par2/par2")
            phi <- priors
        }
        cat("\t", contig, sex, ploidy, "\n")

        hmmdata_path <- sprintf("%s/%s/%s-%s.hmmdata", dir, indiv, indiv, contig)
        filtered_path <- sprintf("%s/%s/%s-%s.filtered.hmmdata",
                                 dir, indiv, indiv, contig)

        if (!file.exists(hmmdata_path)) {
            cat("MISSING file for CONTIG ", contig, " INDIV ", indiv, "\n")
            cat(hmmdata_path, "\n")
            next
        }
        if (use_filtered_hmmdata && !file.exists(filtered_path)) {
            cat("MISSING file for CONTIG ", contig, " INDIV ", indiv, "\n")
            cat(filtered_path, "\n")
            next
        }

        if (use_filtered_hmmdata) {
            data <- read.data(dir, indiv, contig, filtered = TRUE)
        } else {
            data <- read.data(dir, indiv, contig)
            data$read <- factor.contiguous(data$pos)
        }

        ok <- !is.na(data$bad) | !is.na(data$par1ref) & !is.na(data$par2ref) &
            !is.na(data$cons)
        cat("\tRound 2: Removing", sum(!ok),
            "sites at which par1/par2/cons allele unknown\n")
        data$bad[!ok] <- "par1/par2/cons unknown"

        ok <- data$A + data$C + data$G + data$T > 0
        cat("\tRound 2: Removing", sum(!ok),
            "sites at which cons allele is known but reads are unknown\n")
        data$bad[!ok] <- "reads unknown"

        ok <- !is.na(data$bad) | data$par1ref %in% alleles &
            data$par2ref %in% alleles
        data$bad[!ok] <- "par1/par2 not in ACGT"
        cat("\tRound 2: Removing", sum(!ok),
            "sites at which par1/par2 ref not %in% {",
            paste(alleles, collapse = ", "), "}\n")

        badpos <- data$pos[!is.na(data$bad)]
        data <- data[is.na(data$bad), , drop = FALSE]

        ok <- data$par1ref != data$par2ref
        cat("\tRemoving", sum(!ok), "sites at which par1 == par2\n")
        data <- data[ok, , drop = FALSE]

        data$count <- data$A + data$C + data$G + data$T
        ok <- data$count >= minCoverage
        cat("\tRemoving", sum(!ok), "sites at where coverage is < ",
            minCoverage, "\n")
        data <- data[ok, , drop = FALSE]

        if (nrow(data) == 0) next

        if (one.site.per.read) {
            data$read <- factor(data$read)
            ok <- !duplicated(data$read)
            cat("\tRemoving", sum(!ok), "sites from same reads\n")
            data <- data[ok, , drop = FALSE]
            cat("\tNumber of informative markers:", nrow(data), "\n")
        }

        cat("\tFinal total of", nrow(data),
            "sites at which par1 != par2\n")
        if (nrow(data) == 0) next

        L <- nrow(data)
        K <- length(ancestries)

        # Transition probabilities
        if (recrate == 0) {
            if (contig %in% main.contigs) {
                r <- 1 / contigLengths[contig, "length"]
            } else {
                cat("\tContig ", contig,
                    " not found in main.contigs - defaulting to contig length of ",
                    contigLengths[1, "chr"], "\n")
                r <- 1 / contigLengths[1, "length"]
            }
        } else {
            r <- recrate / sum(contigLengths[, "length"])
        }

        d <- c(NA, diff(data$pos))
        p <- 1 - exp(-r * d * rfac)
        Pi <- array(dim = c(L, K, K),
                    dimnames = list(NULL, ancestries, ancestries))
        if (ploidy == 2) {
            Pi[, "par1/par1", "par1/par1"] <-
                Pi[, "par1/par2", "par1/par2"] <-
                Pi[, "par2/par2", "par2/par2"] <- 1 - p
            Pi[, "par1/par1", "par1/par2"] <-
                Pi[, "par1/par2", "par1/par1"] <-
                Pi[, "par1/par2", "par2/par2"] <-
                Pi[, "par2/par2", "par1/par2"] <- p
            Pi[, "par1/par1", "par2/par2"] <-
                Pi[, "par2/par2", "par1/par1"] <- 0
        } else {
            Pi[, "par1", "par1"] <- Pi[, "par2", "par2"] <- 1 - p
            Pi[, "par1", "par2"] <- Pi[, "par2", "par1"] <- p
        }
        Pi[1, , ] <- NA

        # Allele frequencies in parental backgrounds
        ppar1 <- ppar2 <- matrix(NA, nrow = 4, ncol = 4,
                                 dimnames = list(alleles, alleles))
        ppar1[] <- deltapar1 / 3
        diag(ppar1) <- 1 - deltapar1
        ppar2[] <- deltapar2 / 3
        diag(ppar2) <- 1 - deltapar2

        p1 <- ppar1[data$par1ref, , drop = FALSE]
        p2 <- ppar2[data$par2ref, , drop = FALSE]
        p12 <- array(c(p1, p2), dim = c(dim(p1), 2))
        dimnames(p12) <- list(NULL, alleles, NULL)

        # Take all (<=50) reads for each site
        N <- min(max(data$A + data$C + data$G + data$T + data$N), 50)
        eps  <- paste("eps",  seq_len(N), sep = "")
        read <- paste("read", seq_len(N), sep = "")
        y <- data[, c(alleles, "reads", "quals", "par1ref"), drop = FALSE]
        y$selected_allele <- NA
        y[, eps]  <- rep(0, N)
        y[, read] <- rep(5, N)

        allele_to_int <- function(x) {
            switch(x, "A" = 0, "C" = 1, "G" = 2, "T" = 3, "N" = 5)
        }

        for (i in seq_len(nrow(y))) {
            total.reads <- unlist(strsplit(
                cleanupReadPileup(y[i, "reads"], y[i, "par1ref"]), ""))
            y[i, "selected_allele"] <-
                total.reads[sample(length(total.reads), 1)]
            qual <- NULL
            qual_corrected <- NULL
            for (s in seq_len(min(length(total.reads), N))) {
                y[i, read[s]] <- allele_to_int(total.reads[s])
                qual <- c(qual, utf8ToInt(
                    unlist(strsplit(y[i, "quals"], ""))[s]) - 33)
            }
            for (g in seq_along(qual)) {
                qual_corrected[g] <- qual[g] *
                    (theta ^ (rank(-qual)[g] - 1))
                y[i, eps[g]] <- 10 ^ (-(qual_corrected[g]) / 10)
            }
        }

        data$read_allele <- as.vector(y[, "selected_allele"])

        # Emission probabilities (calls hmmprobs binary in scripts/)
        prob <- Pr.y.given.z(y = y[, read, drop = FALSE], p = p12, n = N,
                             eps = y[, eps, drop = FALSE],
                             ploidy = ploidy, C = TRUE,
                             dir = hmmprobs_dir, chrom = contig, id = indiv)
        colnames(prob) <- paste("Pr(y|", ancestries, ")")
        data <- cbind(data, prob)
        data$est <- apply(prob, 1, which.max)

        # Posterior probability
        hmm <- forwardback.ded(Pi = Pi, delta = phi, prob = prob)
        Pr.z.given.y <- exp(hmm$logalpha + hmm$logbeta - hmm$LL)
        colnames(Pr.z.given.y) <- paste("Pr(", ancestries, "|y)")
        data <- cbind(data, Pr.z.given.y)
        attr(data, "badpos") <- badpos
        dataa[[contig]] <- data

        if (pepthresh != "null") {
            cat("Saving contig data as CSV ...")
            write.csv(dataa[[contig]],
                      file = paste(hmmdata.file, "chrom", contig, "csv",
                                   sep = "."))
            cat("OK\n")
        }
    }
    cat("Saving data...")
    save(dataa, file = hmmdata.file)
    cat("OK\n")
}

# ---------------------------------------------------------------------------
# Plot + summary outputs
# ---------------------------------------------------------------------------

contigLengths <- contigLengths[plot.contigs, ]

breakpoints <- {}
matchMismatch <- {}

cat("Plotting...")
plotfile <- file.path(outdir, indiv, paste(indiv, "hmmprob.pdf", sep = "-"))
if (file.exists(plotfile)) {
    cat("plot already exists\n")
    quit(save = "no", status = 0)
}
pdf(file = plotfile, width = 7, height = 1.5)
par(mar = c(2, 2.5, 0.5, 0.5), bg = "transparent",
    cex.main = .68, cex.lab = .8, font.lab = 2, cex.axis = .38,
    mgp = c(1, .000001, 0), xaxs = "i")

plot(0, 0, xlab = "", ylab = "", col = "transparent",
     xlim = c(1, sum(as.numeric(contigLengths$length)) +
              plotPadding * (length(plot.contigs) + 1)),
     ylim = c(-1.01, 1.01), axes = FALSE)

axis(side = 2, at = c(-1, 0, 1), labels = c("", "", ""), col = "gray38")
mtext(c("par2", "par1"), side = 2, line = .68, at = c(-1, 1),
      font = 2, cex = .8, col = c("blue", "red"), las = 2)
box(col = "gray68")

current_start <- plotPadding
for (contig in plot.contigs) {

    mtext(side = 1, at = current_start, contig, font = 2, cex = .8,
          line = 1, xpd = TRUE, adj = 0)
    current_end <- current_start +
        contigLengths[contigLengths$chr == contig, "length"] - 1

    if (sex == "male" && contig %in% sex.chroms) {
        ploidy <- 1
        ancestries <- c("par1", "par2")
        par1homo_col <- 1
        par2homo_col <- 2
    } else {
        ploidy <- 2
        ancestries <- c("par1/par1", "par1/par2", "par2/par2")
        par1homo_col <- 1
        par2homo_col <- 3
    }

    if (sum(names(dataa) %in% contig) != 0) {
        contig_data <- dataa[[contig]]
        x <- contig_data$pos
        y <- contig_data[, paste("Pr(", ancestries, "|y)")]

        byBlocks <- breakpoint.width(x, y[, par1homo_col], y[, par2homo_col],
                                     indiv = indiv, contig = contig,
                                     conf1 = .05, conf2 = .95)
        if (!is.null(byBlocks[["bps"]])) {
            breakpoints <- rbind(breakpoints, byBlocks[["bps"]])
        }

        output_geneious_file(gff_thresh_conf, x, y, par1homo_col,
                             par2homo_col, indiv, contig, contigLengths)

        like.par1 <- contig_data[
            contig_data$read_allele == contig_data$par1ref, ]$pos
        like.par2 <- contig_data[
            contig_data$read_allele == contig_data$par2ref, ]$pos
        plot.posterior(x + current_start, y, ancestries,
                       like.par1 + current_start,
                       like.par2 + current_start,
                       bounds = c(1, contigLengths[
                           contigLengths$chr == contig, ]$length) +
                           current_start - 1,
                       subtract = current_start, tickwidth = 5 * 10 ^ 7)

        if (nrow(byBlocks[["blocks"]]) > 0) {
            matchMismatch <- rbind(matchMismatch,
                reportCounts(contig_data,
                    as.vector(byBlocks[["blocks"]][, "V1"]),
                    as.vector(byBlocks[["blocks"]][, "V4"]),
                    as.numeric(as.vector(byBlocks[["blocks"]][, "V7"])),
                    as.numeric(as.vector(byBlocks[["blocks"]][, "V8"]))))
        }
    }

    current_start <- current_start +
        contigLengths[contigLengths$chr == contig, "length"] + plotPadding
    if (contig != plot.contigs[length(plot.contigs)]) {
        abline(v = current_start - (plotPadding / 2),
               col = "gray68", lwd = 1)
    }
}

dev.off()
cat("OK\n")

if (!is.null(breakpoints)) {
    write.table(breakpoints,
        file = file.path(outdir, indiv,
            paste(indiv, "breakpoints.csv", sep = "-")),
        append = FALSE, quote = FALSE, na = "NA",
        row.names = FALSE, col.names = FALSE, sep = ",")
}

if (!is.null(matchMismatch)) {
    write.table(as.data.frame(matchMismatch),
        file = file.path(outdir, indiv,
            paste(indiv, "matchMismatch.csv", sep = "-")),
        append = FALSE, quote = FALSE, na = "NA",
        row.names = FALSE, col.names = TRUE, sep = ",")
}
