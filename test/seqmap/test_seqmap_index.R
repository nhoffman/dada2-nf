#!/usr/bin/env Rscript

# A seqmap's sv column must name a row of the published seqtab.
#
# dada$map indexes dada$denoised, which is not abundance-sorted, while
# save_seqtab writes the seqtab sorted by abundance. Joining them as if the
# indices agreed mislabels reads silently, so this reconciles them: the reads
# counted under each emitted index must equal that seqtab row's abundance.
#
# Run against any output tree:  test_seqmap_index.R output-its [output-minimal]

suppressWarnings(suppressMessages(library(dada2, quietly=TRUE)))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  args <- c("output-its", "output-minimal", "output-noindex",
            "output-cmsearch", "output-vsearch")
}

failures <- 0
checked <- 0

for (root in args) {
  if (!dir.exists(file.path(root, "dada"))) {
    cat("skipping", root, "- no dada/ directory\n")
    next
  }
  for (p in Sys.glob(file.path(root, "dada", "*", "*"))) {
    rds <- file.path(p, "dada.rds")
    if (!file.exists(rds)) next
    d <- readRDS(rds)

    for (dr in c("r1", "r2")) {
      obj <- if (dr == "r1") d$f else d$r
      if (is.null(obj$dada)) next
      abundance <- unname(obj$dada$denoised)
      ord <- order(-abundance)
      row_of <- integer(length(abundance))
      row_of[ord] <- seq_len(length(abundance))
      index <- row_of[obj$dada$map[obj$derep[[1]]$map]]
      counted <- as.integer(table(factor(index, levels=seq_along(abundance))))
      checked <- checked + 1
      if (!identical(counted, as.integer(abundance[ord]))) {
        cat("FAIL", p, dr, "\n  reads per index:", counted,
            "\n  seqtab rows    :", abundance[ord], "\n")
        failures <- failures + 1
      }
    }

    merged <- d$merged
    if (is.null(merged) || nrow(merged) == 0) next
    checked <- checked + 1
    # Matched on sequence, so a chimera dropped by removeBimeraDenovo is NA
    # rather than shifting every row after it.
    merged_row <- match(merged$sequence, colnames(d$seqtab.nochim))
    present <- na.omit(merged_row)
    if (!all(present >= 1 & present <= ncol(d$seqtab.nochim)) ||
        length(unique(present)) != length(present)) {
      cat("FAIL", p, "merged\n  rows:", merged_row, "\n")
      failures <- failures + 1
    }
  }
}

cat(sprintf("%d checks, %d failures\n", checked, failures))
if (checked == 0) {
  cat("nothing checked - no dada.rds in:",
      paste(args, collapse=" "), "\n")
  quit(status=1)
}
if (failures > 0) quit(status=1)
