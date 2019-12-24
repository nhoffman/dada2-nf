#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(argparse, quietly = TRUE)))
suppressWarnings(suppressMessages(library(dada2, quietly = TRUE)))

do_dada <- function(filtered, sample.names, err=NULL, ...){
  derep <- dada2::derepFastq(filtered)

  if(class(derep) == 'derep'){
    ## list instead of bare derep object for single sample
    derep=list(derep)
  }
  names(derep) <- sample.names
  dada <- dada2::dada(derep, err=err, ...)
  list(derep=derep, dada=dada)
}


getN <- function(x){
  sum(dada2::getUniques(x))
}


main <- function(arguments){
  parser <- ArgumentParser()
  parser$add_argument('r1', help='path to fastq for R1')
  parser$add_argument('r2', help='path to fastq for R2')
  parser$add_argument('--errors',
                      help=paste('.rds file containing error model',
                                 '(a list with values "errF" and "errR").',
                                 'Calculates error model for this specimen if not provided.'))

  ## outputs
  parser$add_argument('--data', default='dada.rds',
                      help="output .rds file containing intermediate data structures")
  parser$add_argument('--seqtab', default='seqtab.csv',
                      help="output file containing chimera-checked SVs")
  parser$add_argument('--counts', default='counts.csv',
                      help="input and output read counts")
  parser$add_argument('--overlaps', default='overlaps.csv',
                      help="distribution of overlaps among merged reads")


  ## parameters
  parser$add_argument('-s', '--sampleid', default='unknown', help='label for this specimen')

  ## params for dada2::dada()
  parser$add_argument('--self-consist', default='TRUE', choices=c('TRUE', 'FALSE'),
                      help='value for dad2::dada(selfConsist) [%(default)s]')
  parser$add_argument('--max-indels', type='integer', default=16,
                      help=paste('sets dada2::dada(BAND_SIZE), ',
                                 'essentially the maximum number of cumulative indels',
                                 'a sequence may contain compared to a more',
                                 'abundant variant to be considered for grouping.',
                                 '[%(default)s]'))
  parser$add_argument('--omega-a', type='double',
                      help='dada2::dada(OMEGA_A)')
  parser$add_argument('--omega-c', type='double',
                      help='dada2::dada(OMEGA_C)')

  ## dada2::mergePairs()
  parser$add_argument('--max-mismatch', type='integer', default=0,
                      help='The maximum mismatches allowed in overlap [%(default)s]')

  ## dada2::removeBimeraDenovo()
  parser$add_argument('--min-fold-over-parent', type='integer', default=2,
                      help='value for removeBimeraDenovo(minFoldParentOverAbundance) [%(default)s]')
  parser$add_argument('--chimera-method', default='consensus',
                      help='value for removeBimeraDenovo(method) [%(default)s]')

  parser$add_argument('--nthreads', type='integer', default=0,
                      help='number of processes; defaults to number available')

  args <- parser$parse_args(arguments)
  multithread <- if(args$nthreads == 0){ TRUE }else{ args$nthreads }

  fnFs <- args$r1
  fnRs <- args$r2

  ## if the inputs are empty, create outputs and exit
  if(length(readLines(gzfile(fnFs, 'rt'), n=1)) == 0){
    cat(gettextf('%s is empty\n', fnFs))

    file.create(args$seqtab)  ## an empty file

    saveRDS(list(
        sampleid=args$sampleid, f=NULL, r=NULL, merged=NULL,
        seqtab=NULL, seqtab.nochim=NULL
    ),
    file=args$data)

    ## read counts for various stages of the analysis
    write.csv(data.frame(
        sampleid=args$sampleid,
        filtered_and_trimmed=0,
        denoised_r1=0,
        denoised_r2=0,
        merged=0,
        no_chimeras=0
    ),
    file=args$counts, row.names=FALSE)

    ## overlaps
    write.csv(data.frame(
        sampleid=args$sampleid, nmatch=NA, abundance=NA
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
  f <- do_dada(fnFs,
               sample.names=args$sampleid,
               err=errors$errF,
               ## additional arguments for dada2::dada
               multithread=multithread,
               selfConsist=as.logical(args$self_consist),
               BAND_SIZE=args$max_indels)

  cat('dereplicating and applying error model for reverse reads\n')
  r <- do_dada(fnRs,
               sample.names=args$sampleid,
               err=errors$errR,
               ## additional arguments for dada2::dada
               multithread=multithread,
               selfConsist=as.logical(args$self_consist),
               BAND_SIZE=args$max_indels)

  cat('merging reads\n')
  merged <- dada2::mergePairs(
                       dadaF=f$dada,
                       derepF=f$derep,
                       dadaR=r$dada,
                       derepR=r$derep,
                       maxMismatch=args$max_mismatch,
                       verbose=TRUE)

  if(nrow(merged) > 0){
    cat('making sequence table\n')
    seqtab <- dada2::makeSequenceTable(merged)
    rownames(seqtab) <- args$sampleid

    cat('checking for chimeras\n')
    seqtab.nochim <- dada2::removeBimeraDenovo(
                                seqtab, multithread=multithread, verbose=TRUE)
    rownames(seqtab.nochim) <- args$sampleid

    write.table(
        data.frame(sampleid=args$sampleid,
                   count=as.integer(seqtab.nochim),
                   seq=colnames(seqtab.nochim)),
        file=args$seqtab,
        sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)

    ## saveRDS(seqtab.nochim, file=args$seqtab)
    saveRDS(list(sampleid=args$sampleid, f=f, r=r, merged=merged,
                 seqtab=seqtab, seqtab.nochim=seqtab.nochim),
            file=args$data)

    ## read counts for various stages of the analysis
    counts <- data.frame(
        sampleid=args$sampleid,
        filtered_and_trimmed=getN(f$derep[[1]]),
        denoised_r1=getN(f$dada),
        denoised_r2=getN(r$dada),
        merged=getN(merged),
        no_chimeras=rowSums(seqtab.nochim)
    )
    write.csv(counts, file=args$counts, row.names=FALSE)

    ## calculate overlaps among merged reads
    overlaps <- data.frame(aggregate(abundance ~ nmatch, merged, sum))
    overlaps$sampleid <- args$sampleid
    write.csv(overlaps[, c('sampleid', 'nmatch', 'abundance')],
              file=args$overlaps, row.names=FALSE)
  }else{
    cat(gettextf('Warning: no merged reads in sample %s\n', args$sampleid))

    file.create(args$seqtab)  ## an empty file

    saveRDS(list(
        sampleid=args$sampleid, f=f, r=r, merged=NULL,
        seqtab=NULL, seqtab.nochim=NULL
    ),
    file=args$data)

    ## read counts for various stages of the analysis
    write.csv(data.frame(
        sampleid=args$sampleid,
        filtered_and_trimmed=getN(f$derep[[1]]),
        denoised_r1=getN(f$dada),
        denoised_r2=getN(r$dada),
        merged=0,
        no_chimeras=0
    ),
    file=args$counts, row.names=FALSE)

    write.csv(data.frame(
        sampleid=args$sampleid, nmatch=NA, abundance=NA
    ),
    file=args$overlaps, row.names=FALSE)
  }

}

main(commandArgs(trailingOnly=TRUE))
## invisible(warnings())

