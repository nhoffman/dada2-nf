#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(argparse, quietly = TRUE)))
suppressWarnings(suppressMessages(library(dada2, quietly = TRUE)))

main <- function(arguments){
  parser <- ArgumentParser()
  parser$add_argument('--infiles', nargs='+', help='input R1,R2 fq.gz')
  parser$add_argument('--outfiles', nargs='+', help='output R1,R2 fq.gz')
  parser$add_argument('--f-trunc', type='integer', required=TRUE,
                      help='position at which to truncate forward reads')
  parser$add_argument('--r-trunc', type='integer', required=TRUE,
                      help='position at which to truncate reverse reads')
  parser$add_argument('--trim-left', type='integer', required=TRUE,
                      help='position at which to left-trim both F and R reads')
  parser$add_argument('--max-ee', type='integer',
                      help=paste(
                          'After truncation, reads with higher than maxEE ',
                          '"expected errors" will be discarded (default is ',
                          'no filtering using this parameter)'))
  parser$add_argument('--truncq', type='integer', metavar='N', default=2,
                      help=paste('truncate reads at the first instance ',
                                 'of a quality score <= N'))
  parser$add_argument('--nthreads', type='integer', default=0,
                      help='number of processes; defaults to number available')

  args <- parser$parse_args(arguments)
  trim_left <- rep(args$trim_left, 2)
  maxEE <- if(is.null(args$max_ee)){ Inf }else{ args$max_ee }
  truncQ <- args$truncq
  multithread <- if(args$nthreads == 0){ TRUE }else{ args$nthreads }

  dada2::fastqPairedFilter(fn=args$infiles,
                           fout=args$outfiles,
                           trimLeft=trim_left,
                           truncLen=c(args$f_trunc, args$r_trunc),
                           maxN=0,
                           maxEE=maxEE,
                           truncQ=truncQ,
                           compress=TRUE,
                           multithread=multithread,
                           verbose=TRUE)

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

