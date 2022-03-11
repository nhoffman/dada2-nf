#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse, quietly = TRUE))


get_unmerged <- function(obj, orientation){
  stopifnot(orientation %in% c('forward', 'reverse'))
  ori = substring(orientation, 1, 1)  # f or r

  ## denoised sequences
  denoised <- obj[[ori]]$dada$sequence
  abundance <- obj[[ori]]$dada[1]$denoised

  ## indices of merged reads
  merged_idx <- obj$merged[[orientation]]

  ## determing zero-padding/naming for use later on in fasta names
  padchars <- ceiling(log10(length(denoised) + 1))
  seqnames <- gettextf(paste0('%0', padchars, 'i'), seq(length(denoised)))

  if(!is.null(merged_idx)){
    ## return a vector of fasta records
    gettextf('>%s%s:%s\n%s', ori, seqnames[-merged_idx], abundance[-merged_idx], denoised[-merged_idx])
  }
}

main <- function(arguments){
  parser <- ArgumentParser(
      description="Write forward and reverse unmerged denoised, dereplicated reads")
  parser$add_argument('rdata', help='RDS files containing dada2 output')
  parser$add_argument(
             '-f', '--forward-seqs', default='unmerged_F.fasta',
             help='output fasta file with unmerged forward sequence variants')
  parser$add_argument(
             '-r', '--reverse-seqs', default='unmerged_R.fasta',
             help='output fasta file with unmerged reverse sequence variants')

  args <- parser$parse_args(arguments)
  obj <- readRDS(args$rdata)

  if(!is.null(args$forward_seqs)){
    writeLines(get_unmerged(obj, 'forward'), args$forward_seqs)
  }

  if(!is.null(args$reverse_seqs)){
    writeLines(get_unmerged(obj, 'reverse'), args$reverse_seqs)
  }
}

main(commandArgs(trailingOnly=TRUE))
## invisible(warnings())
