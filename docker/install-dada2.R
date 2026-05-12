library(parallel)
ncores <- min(c(8, parallel::detectCores()))

dada2_commit <- Sys.getenv('DADA2_REF')
if(nchar(dada2_commit) == 0){
  stop('the environment variable DADA2_REF must be set')
}

cran_packages <- c(
    "ape",
    "argparse",
    "dplyr",
    "ggplot2",
    "gridExtra",
    "jsonlite",
    "lattice",
    "latticeExtra",
    "phyloseq",
    "R.utils",
    "readr",
    "remotes",
    "rmarkdown",
    "qrqc",
    "tibble",
    "tidyr"
)

install.packages(
    cran_packages,
    repos="https://cloud.r-project.org",
    Ncpus=ncores,
    clean=TRUE)

# ref: Desired git reference. Could be a commit, tag, or branch name,
# or a call to github_pull. Defaults to "master".
remotes::install_github(
    "benjjneb/dada2",
    ref=dada2_commit,
    threads=ncores)
