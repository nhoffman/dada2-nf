#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(argparse, quietly = TRUE)))
suppressWarnings(suppressMessages(library(jsonlite, quietly = TRUE)))
suppressWarnings(suppressMessages(library(dada2, quietly = TRUE)))

getN <- function(x){
  sum(dada2::getUniques(x))
}

save_seqtab <- function(filename, dada, sampleid, orientation){
  if(is.null(dada)){
    file.create(filename)
  }else{
    if(orientation == 'reverse'){
      seq = dada2::rc(dada$sequence)
    } else{
      seq = dada$sequence
    }
    df <- data.frame(
        sampleid=sampleid,
        abundance=dada[1]$denoised,
        seq=seq
    )
    write.table(
       df[order(-df$abundance),],
       file=filename,
       sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)
  }
}

main <- function(arguments){
  parser <- ArgumentParser()
  parser$add_argument('r1', help='path to fastq for R1')
  parser$add_argument('r2', help='path to fastq for R2')
  parser$add_argument('--errors',
                      help=paste('.rds file containing error model',
                                 '(a list with values "errF" and "errR").',
                                 'Calculates error model for this specimen ',
                                 'if not provided.'))

  ## outputs
  parser$add_argument('--data', default='dada.rds',
                      help="output .rds file containing intermediate data structures")
  parser$add_argument('--seqtab', default='seqtab.csv',
                      help="output file containing merged, chimera-checked SVs")
  parser$add_argument('--seqtab-r1', default='seqtab_r1.csv',
                      help="output file containing denoised reads for R1")
  parser$add_argument('--seqtab-r2', default='seqtab_r2.csv',
                      help="output file containing denoised reads for R2")
  parser$add_argument('--counts', default='counts.csv',
                      help="input and output read counts")
  parser$add_argument('--overlaps', default='overlaps.csv',
                      help="distribution of overlaps among merged reads")

  ## parameters
  parser$add_argument('-s', '--sampleid', default='unknown',
                      help='label for this specimen')
  parser$add_argument('--orientation', default='unknown',
                      help='orientation for this specimen')
  parser$add_argument('--params',
                      help=paste(
                          'json file containing optional parameters for',
                          'dada(), mergePairs(), and removeBimeraDenovo (see README)'))
  parser$add_argument('--nthreads', type='integer', default=0,
                      help='number of processes; defaults to number available')

  args <- parser$parse_args(arguments)
  multithread <- if(args$nthreads == 0){ TRUE }else{ args$nthreads }

  if(is.null(args$params)){
    params <- list()
  }else{
    params <- fromJSON(args$params)
  }

  fnFs <- args$r1
  fnRs <- args$r2

  ## if the inputs are empty, create outputs and exit
  if(length(readLines(gzfile(fnFs, 'rt'), n=1)) == 0){
    cat(gettextf('%s is empty\n', fnFs))

    file.create(args$seqtab)  ## an empty file
    file.create(args$seqtab_r1)  ## an empty file
    file.create(args$seqtab_r2)  ## an empty file

    saveRDS(list(
        sampleid=args$sampleid, f=NULL, r=NULL, merged=NULL,
        seqtab=NULL, seqtab.nochim=NULL
    ),
    file=args$data)

    ## read counts for various stages of the analysis
    write.csv(data.frame(
        sampleid=args$sampleid,
        orientation=args$orientation,
        denoised_r1=0,
        denoised_r2=0,
        merged=0,
        no_chimeras=0
    ),
    file=args$counts, row.names=FALSE)

    ## overlaps
    write.csv(data.frame(
        sampleid=args$sampleid, orientation=args$orientation, nmatch=NA, abundance=NA
    ),
    file=args$overlaps, row.names=FALSE)

    quit(status=0)
  }

  if(is.null(args$errors)){
    errors <- list()
  }else{
    cat(gettextf('using errors in %s\n', args$errors))
    errors <- readRDS(args$errors)
  }

  cat('dereplicating and applying error model for forward reads\n')
  dada_params <- if(is.null(params$dada)){list()}else{params$dada}

  ## list instead of bare derep object for single sample
  derepF <- setNames(list(dada2::derepFastq(fnFs)), args$sampleid)
  paramsF <- modifyList(
      c(list(derep=derepF), dada_params),
      if(is.null(errors$errF)){
        list(selfConsist=TRUE)
      }else{
        list(err=errors$errF)
      })

  dadaF <- tryCatch(
      do.call(dada2::dada, c(paramsF, list(multithread=multithread))),
      error=function(err){
        cat('Error:', err$message, '\n')
        cat('saving NULL object for forward reads\n')
        NULL
      })

  save_seqtab(filename=args$seqtab_r1, dada=dadaF, sampleid=args$sampleid, orientation=args$orientation)

  cat('dereplicating and applying error model for reverse reads\n')
  derepR <- setNames(list(dada2::derepFastq(fnRs)), args$sampleid)
  paramsR <- modifyList(
      c(list(derep=derepR), dada_params),
      if(is.null(errors$errR)){
        list(selfConsist=TRUE)
      }else{
        list(err=errors$errR)
      })
  dadaR <- tryCatch(
      do.call(dada2::dada, c(paramsR, list(multithread=multithread))),
      error=function(err){
        cat('Error:', err$message, '\n')
        cat('saving NULL object for reverse reads\n')
        NULL
      })

  save_seqtab(filename=args$seqtab_r2, dada=dadaR, sampleid=args$sampleid, orientation=args$orientation)

  if(is.null(dadaF) || is.null(dadaR)){
    merged <- NULL
  }else{
    cat('merging reads\n')
    merge_args <- modifyList(
        list(dadaF=dadaF,
             derepF=derepF,
             dadaR=dadaR,
             derepR=derepR),
        if(is.null(params$mergePairs)){list()}else{params$mergePairs})

    merged <- do.call(dada2::mergePairs, merge_args)
  }

  ## success
  if(!is.null(merged) && nrow(merged) > 0){
    cat('making sequence table\n')
    seqtab <- dada2::makeSequenceTable(merged)
    rownames(seqtab) <- args$sampleid

    cat('checking for chimeras\n')
    bimera_args <- modifyList(
        list(unqs=seqtab,
             multithread=multithread),
        if(is.null(params$removeBimeraDenovo)){list()}else{params$removeBimeraDenovo})

    seqtab.nochim <- do.call(dada2::removeBimeraDenovo, bimera_args)
    rownames(seqtab.nochim) <- args$sampleid

    if(args$orientation == 'reverse'){
        seq = dada2::rc(colnames(seqtab.nochim))
    } else{
        seq = colnames(seqtab.nochim)
    }

    ## csv with merged reads
    write.table(
        data.frame(sampleid=args$sampleid,
                   count=as.integer(seqtab.nochim),
                   seq=seq),
        file=args$seqtab,
        sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)

    ## saveRDS(seqtab.nochim, file=args$seqtab)
    saveRDS(list(sampleid=args$sampleid,
                 f=list(derep=derepF, dada=dadaF),
                 r=list(derep=derepR, dada=dadaR),
                 merged=merged,
                 seqtab=seqtab,
                 seqtab.nochim=seqtab.nochim),
            file=args$data)

    ## read counts for various stages of the analysis
    counts <- data.frame(
        sampleid=args$sampleid,
        orientation=args$orientation,
        denoised_r1=getN(dadaF),
        denoised_r2=getN(dadaR),
        merged=getN(merged),
        no_chimeras=rowSums(seqtab.nochim)
    )
    write.csv(counts, file=args$counts, row.names=FALSE)

    ## calculate overlaps among merged reads
    overlaps <- data.frame(aggregate(abundance ~ nmatch, merged, sum))
    overlaps$sampleid <- args$sampleid
    overlaps$orientation <- args$orientation
    write.csv(overlaps[, c('sampleid', 'orientation', 'nmatch', 'abundance')],
              file=args$overlaps, row.names=FALSE)
  }else{
    cat(gettextf('Warning: no merged reads in sample %s\n', args$sampleid))

    file.create(args$seqtab)  ## an empty file

    saveRDS(list(
        sampleid=args$sampleid,
        f=list(derep=derepF, dada=dadaF),
        r=list(derep=derepR, dada=dadaR),
        merged=NULL,
        seqtab=NULL,
        seqtab.nochim=NULL
    ),
    file=args$data)

    ## read counts for various stages of the analysis
    write.csv(data.frame(
        sampleid=args$sampleid,
        orientation=args$orientation,
        denoised_r1=if(is.null(dadaF)){0}else{getN(dadaF)},
        denoised_r2=if(is.null(dadaR)){0}else{getN(dadaR)},
        merged=0,
        no_chimeras=0
    ),
    file=args$counts, row.names=FALSE)

    ## overlaps
    write.csv(data.frame(
        sampleid=args$sampleid, orientation=args$orientation, nmatch=NA, abundance=NA
    ),
    file=args$overlaps, row.names=FALSE)
  }

}

main(commandArgs(trailingOnly=TRUE))
## invisible(warnings())

