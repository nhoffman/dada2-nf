#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(argparse, quietly = TRUE)))
suppressWarnings(suppressMessages(library(jsonlite, quietly = TRUE)))
suppressWarnings(suppressMessages(library(dada2, quietly = TRUE)))

main <- function(arguments){
  parser <- ArgumentParser()
  parser$add_argument('--infiles', nargs='+', help='input R1,R2 fq.gz')
  parser$add_argument('--params',
                      help=paste(
                          'json file containing optional parameters for fastqPairedFilter',
                          '(see README)'))
  parser$add_argument('--outfiles', nargs='+', help='output R1,R2 fq.gz')
  parser$add_argument('--nthreads', type='integer', default=0,
                      help='number of processes; defaults to number available')

  args <- parser$parse_args(arguments)
  multithread <- if(args$nthreads == 0){ TRUE }else{ args$nthreads }

  if(is.null(args$params)){
    params <- list()
  }else{
    params <- fromJSON(args$params)$fastqPairedFilter
  }

  ## params overrides any values in the first argument
  funcargs <- modifyList(
      list(fn=args$infiles,
           fout=args$outfiles,
           compress=TRUE,
           multithread=multithread,
           verbose=TRUE),
      params)

  do.call(dada2::fastqPairedFilter, funcargs)

  ## There is no error if all reads have been eliminated, but no files
  ## are written in this case. Check for the output files, and if they
  ## don't exist, create empty ones.
  for(fn in args$outfiles){
    if(!file.exists(fn)){
      cat(gettextf('creating empty file %s\n', fn))
      gzf = gzfile(fn)
      cat('', file=gzf, fill=FALSE)
      close(gzf)
    }
  }

}

main(commandArgs(trailingOnly=TRUE))
## invisible(warnings())

