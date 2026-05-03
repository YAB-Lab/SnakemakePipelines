#!/usr/bin/env Rscript
#
# Build per-contig HMM input data from BWA-MEM pileups.
#
# Replaces the legacy msg/write-hmm-data.R. Differences:
#   - Reads modern 6-column samtools mpileup output (chrom, pos, ref, depth,
#     bases, quals) instead of the legacy 10-column samtools 0.1.x format.
#   - Native R allele counting (replaces the msg/countalleles C binary).
#   - optparse CLI; no source(ded.R) dependency.
#   - Self-validates: every contig that has *-ref.alleles must produce an
#     .hmmdata file. Missing files cause non-zero exit.
#
# Output (per contig, one file): {indiv}-{contig}.hmmdata
#   Columns: pos ref cons reads quals A C G T N bad par1ref par2ref

suppressPackageStartupMessages(library(optparse))

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

opt_list <- list(
    make_option(c("-i", "--individual"), type = "character",
                help = "Sample / individual ID"),
    make_option(c("-d", "--dir"), type = "character",
                help = "Sample directory containing pileups + refs/"),
    make_option(c("-c", "--chroms"), type = "character", default = "all",
                help = "Comma-separated contigs, or 'all' [default %default]")
)
opts <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opts$individual), !is.null(opts$dir))

indiv <- opts$individual
dir   <- opts$dir
species       <- c("par1", "par2")
pupsp         <- "par1"
allele.states <- c("A", "C", "G", "T", "N")

if (opts$chroms == "all") {
    refs_par1 <- list.files(file.path(dir, "refs", "par1"),
                            pattern = "-ref\\.alleles$")
    contigs <- sub("-ref\\.alleles$", "", refs_par1)
} else {
    contigs <- strsplit(opts$chroms, ",")[[1]]
}

cat("Writing HMM input data for", length(contigs), "contigs.\n")
if (length(contigs) == 0) {
    stop("No contigs found in ", file.path(dir, "refs", "par1"))
}

# ---------------------------------------------------------------------------
# Pileup base decoder (ported verbatim from legacy write-hmm-data.R)
# ---------------------------------------------------------------------------

decode_pileup_bases <- function(x, ref) {
    # x: vector of pileup base codes (e.g. c..c   .+2cc   ^fa^f,^f,)
    # ref: vector of reference alleles (one per element of x)
    # http://samtools.sourceforge.net/pileup.shtml
    x <- gsub("\\^.", "", x)   # begin contiguous read
    x <- gsub("\\$",  "", x)   # end contiguous read

    # Strip indels of varying length
    x <- gsub("^[\\+-][ACGTNXMRWSYKVHDBacgtnxmrwsykvhdb]+", "", x)
    x <- gsub(".*\\*.*", "", x)
    len <- 1
    while (length(grep("[\\+-]", x)) > 0) {
        re <- paste0("[\\+-]", len,
                     paste(rep("[ACGTNXMRWSYKVHDBacgtnxmrwsykvhdb]", len),
                           collapse = ""))
        x <- gsub(re, "", x)
        len <- len + 1
    }

    # Replace . and , with reference base, uppercase
    x <- mapply(gsub, pattern = list("[\\.,]"), replacement = ref, x = x)
    x <- toupper(x)

    # Allele counts via base-R (replaces the msg/countalleles C binary)
    counts <- vapply(allele.states, function(b) {
        nchar(gsub(sprintf("[^%s]", b), "", x))
    }, integer(length(x)))
    if (length(x) == 1) counts <- matrix(counts, nrow = 1)
    colnames(counts) <- allele.states
    counts
}

# ---------------------------------------------------------------------------
# Per-contig processing
# ---------------------------------------------------------------------------

read_pileup_modern <- function(pupfile) {
    # Modern samtools mpileup: chrom pos ref depth bases quals
    # We need pos (col 2), ref (col 3), bases (col 5), quals (col 6).
    # Use cut -f to skip lines with extra fields cleanly (indels can pad).
    pipa <- pipe(paste("cut -f2,3,5,6 <", shQuote(pupfile)))
    raw <- scan(pipa, what = "", sep = "\n", quiet = TRUE)
    close(pipa)
    if (length(raw) == 0) return(NULL)
    parts <- strsplit(raw, "\t")
    stopifnot(all(lengths(parts) == 4))
    m <- matrix(unlist(parts), ncol = 4, byrow = TRUE)
    data.frame(
        pos   = as.integer(m[, 1]),
        ref   = toupper(m[, 2]),
        cons  = toupper(m[, 2]),  # placeholder for fit-hmm.R's is.na(cons) check
        reads = m[, 3],
        quals = m[, 4],
        stringsAsFactors = FALSE
    )
}

written_files <- character(0)

for (contig in contigs) {
    cat("Writing HMM input data for", indiv, contig, "\n")

    outfile <- sprintf("%s/%s-%s.hmmdata", dir, indiv, contig)
    if (file.exists(outfile)) {
        cat("File exists -- skipping:", outfile, "\n")
        written_files <- c(written_files, outfile)
        next
    }

    pupfile <- sprintf("%s/aln_%s_%s-filtered-%s-sorted.pileup",
                       dir, indiv, pupsp, contig)
    if (!file.exists(pupfile)) {
        cat("MISSING ", pupfile, "\n")
        next
    }

    pup <- read_pileup_modern(pupfile)
    if (is.null(pup) || nrow(pup) == 0) {
        cat("\tEmpty pileup -- skipping\n")
        next
    }

    cat("\tStarting with a total of", nrow(pup), "positions.\n")
    ok <- pup$ref %in% allele.states
    cat("\tRemoving", sum(!ok), "positions at which ref is not [ACGT].\n")
    pup <- pup[ok, , drop = FALSE]
    pup <- pup[grep("^\\*+$", pup$reads, invert = TRUE), , drop = FALSE]
    cat("\tKeeping", nrow(pup), "non-indels in par1...")

    cat("\tDecoding read bases...\n")
    tab <- decode_pileup_bases(pup$reads, pup$ref)
    if (nrow(tab) != nrow(pup)) {
        stop("nrows differ after decoding bases: tab=", nrow(tab),
             " pup=", nrow(pup))
    }
    pup <- cbind(pup, tab)
    pup$bad <- rep("", nrow(pup))

    for (sp in species) {
        reffile <- sprintf("%s/refs/%s/%s-ref.alleles", dir, sp, contig)
        ref <- read.delim(reffile, header = FALSE, as.is = TRUE,
                          col.names = c("pos", "allele"))

        # remove positions that are indels in either par1 or par2 (based on ref)
        non_indels_pos <- intersect(ref$pos, pup$pos)
        pup <- pup[pup$pos %in% non_indels_pos, , drop = FALSE]
        ref <- ref[ref$pos %in% non_indels_pos, , drop = FALSE]
        cat("\tKeeping", length(non_indels_pos), "non-indels...\n")

        if (nrow(ref) != nrow(pup)) {
            map <- match(pup$pos, ref$pos)
            ref <- ref[map, ]
        }
        stopifnot(nrow(ref) == nrow(pup), ref$pos == pup$pos)

        spref <- sprintf("%sref", sp)
        pup[[spref]] <- ref$allele

        ok <- !is.na(pup[[spref]])
        cat("\tRemoving", sum(!ok), "positions at which", sp,
            "ref was NA or N.\n")
        pup$bad[!ok] <- paste(sp, "NA/N")

        ok <- pup[[spref]] %in% c("A", "C", "G", "T")
        cat("\tRemoving", sum(!ok), "positions at which", sp,
            "ref is not [ACGT].\n")
        pup$bad[!ok] <- paste(sp, "not ACGT")
    }

    ok <- !is.na(pup$ref)
    cat("\tRemoving", sum(!ok),
        "positions at which par1 ref (pileup) is NA\n")
    pup$bad[!ok] <- "par1 ref NA"

    ok <- pup$ref %in% c("A", "C", "G", "T")
    cat("\tRemoving", sum(!ok),
        "positions at which par1 ref (pileup) is not [ACGT].\n")
    pup$bad[!ok] <- "par1 ref not ACGT"

    ok <- pup$ref == pup$par1ref
    cat("\tRemoving", sum(!ok),
        "positions at which par1 ref (pileup) disagrees with par1 ref.\n")
    pup$bad[!ok] <- "par1 ref disagree"

    cat("\tAfter filtering, keeping", nrow(pup), "positions.\n")
    write.table(pup, file = outfile, sep = "\t",
                col.names = TRUE, row.names = FALSE, quote = FALSE)
    written_files <- c(written_files, outfile)
}

# ---------------------------------------------------------------------------
# Self-validate: every contig with a -ref.alleles file should have an .hmmdata
# ---------------------------------------------------------------------------

expected <- file.path(dir, sprintf("%s-%s.hmmdata", indiv, contigs))
missing  <- expected[!file.exists(expected)]

if (length(missing) > 0) {
    cat("\nMissing .hmmdata files:\n")
    for (m in missing) cat("  ", m, "\n")
    quit(status = 2)
}

cat("\nDone. Wrote", length(written_files), ".hmmdata files.\n")
