library(parallel)
ncores <- min(c(8, parallel::detectCores()))

dada2_commit <- Sys.getenv('DADA2_REF')
if(nchar(dada2_commit) == 0){
  stop('the environment variable DADA2_REF must be set')
}

cran_packages <- c(
    "argparse"
)

install.packages(
    cran_packages,
    repos="http://cran.us.r-project.org",
    Ncpus=ncores,
    clean=TRUE)

# ref: Desired git reference. Could be a commit, tag, or branch name,
# or a call to github_pull. Defaults to "master".
devtools::install_github(
    "benjjneb/dada2",
    ref=dada2_commit,
    threads=ncores)

bioc_packages <- c(
    "qrqc",
    "phyloseq"
)

BiocManager::install(bioc_packages, Ncpus=ncores)

