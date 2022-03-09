#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(argparse, quietly = TRUE)))
suppressWarnings(suppressMessages(library(tidyr, quietly = TRUE)))
suppressWarnings(suppressMessages(library(dplyr, quietly = TRUE)))
suppressWarnings(suppressMessages(library(purrr, quietly = TRUE)))

main <- function(arguments){
  parser <- ArgumentParser()
  parser$add_argument('tabs', nargs='+')
  parser$add_argument('-o', '--outfile', default='joined.csv')

  args <- parser$parse_args(arguments)
  joined <- lapply(args$tabs, read.csv, colClasses='character', check.names=FALSE) %>%
    purrr::reduce(dplyr::left_join, by=colnames(.)[1]) %>%
    tidyr::replace_na(as.list(sapply(colnames(.), function(x){'0'})))

  write.csv(joined, file=args$outfile, row.names=FALSE)
}

main(commandArgs(trailingOnly=TRUE))
## invisible(warnings())

