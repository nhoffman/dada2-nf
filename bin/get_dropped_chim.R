#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse, quietly = TRUE))
suppressPackageStartupMessages(library(tidyr, quietly = TRUE))
suppressPackageStartupMessages(library(tibble, quietly = TRUE))
suppressPackageStartupMessages(library(readr, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))

main <- function(arguments){
  parser <- ArgumentParser(
      description="Write forward and reverse unmerged denoised, dereplicated reads")
  parser$add_argument('rdata', help='RDS files containing dada2 output')
  parser$add_argument(
             '-o', '--outfile', default='chim_dropped.csv',
             help='output csv file with weight and sequence of svs identified as chimeras')

  args <- parser$parse_args(arguments)

  obj <- readRDS(args$rdata)
  if (!is.null(obj$seqtab) && !is.null(obj$seqtab.nochim)) {
    seqtab <- as.data.frame(as.table(obj$seqtab))
    seqtab.nochim <- as.data.frame(as.table(obj$seqtab.nochim))

    ## use of anti_join returns records in seqtab which are not
    ## present in the filtered seqtab.nochim. In other words,
    ## detect and output the 'dropped' seqs and their weights
    tab <- seqtab %>% anti_join(seqtab.nochim, by=c('Var2')) %>%
      arrange(-Freq) %>% rename(sequence=Var2) %>%
      rename(weight=Freq) %>% select(weight, sequence)
  } else {
    tab <- tibble(weight=numeric(),sequence=character())
  }
  tab %>% write_csv(args$outfile)

main(commandArgs(trailingOnly=TRUE))
