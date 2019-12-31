#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(ggplot2, quietly = TRUE)))
suppressWarnings(suppressMessages(library(gridExtra, quietly = TRUE)))
suppressWarnings(suppressMessages(library(argparse, quietly = TRUE)))
suppressWarnings(suppressMessages(library(dada2, quietly = TRUE)))


errplot <- function(err, title=""){
  dada2::plotErrors(err, nominalQ=TRUE) +
    ggtitle(title) +
    theme(aspect.ratio=1)
}


main <- function(arguments){
  parser <- ArgumentParser()
  parser$add_argument('--r1', help='file listing R1 fq files for this batch')
  parser$add_argument('--r2', help='file listing R2 fq files for this batch')
  parser$add_argument('--model', default='error_model.rds',
                      help='output .rds file containing the model')
  parser$add_argument('--plots',
                      help='optional file for output of dada2::plotErrors()')
  parser$add_argument('--nthreads', type='integer', default=0,
                      help='number of processes; defaults to number available')

  args <- parser$parse_args(arguments)
  multithread <- if(args$nthreads == 0){ TRUE }else{ args$nthreads }

  fnFs <- readLines(args$r1)
  fnRs <- readLines(args$r2)

  cat('generating error model for forward reads\n')
  errF <- tryCatch(
      dada2::learnErrors(fnFs, multithread=multithread),
      error=function(err){
        cat('Error:', err$message, '\n')
        cat('saving NULL error model for forward reads\n')
        NULL
      })

  cat('generating error model for reverse reads\n')
  errR <- tryCatch(
      dada2::learnErrors(fnRs, multithread=multithread),
      error=function(err){
        cat('Error:', err$message, '\n')
        cat('saving NULL error model for reverse reads\n')
        NULL
      })

  cat(gettextf('saving error model to %s\n', args$model))
  saveRDS(list(errF=errF, errR=errR), file=args$model)

  if(is.null(args$plots) | is.null(errF)){
    cat('Saving empty file', args$plots, '\n')
    file.create(args$plots)  ## an empty file
  }
  else{
    fig <- gridExtra::grid.arrange(errplot(errF, "forward"),
                                   errplot(errR, "reverse"),
                                   nrow=1)
    ggplot2::ggsave(args$plots, fig)
  }
}

main(commandArgs(trailingOnly=TRUE))
## invisible(warnings())

