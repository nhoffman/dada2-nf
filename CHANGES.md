# Change log for dada2-nf

## 1.18

- Optionally include cutadapt to remove adapters (GH55)
- Docker image is built on Ubuntu 22.04 base (GH56)
- Docker image is hosted by GutHub container registry
- Simplify Docker and Singularity configuration; there is no longer a
  Singularity image associated with this project (both Docker and
  Singularity use the same image hosted on ghcr).
- Update GH action to build image and test sequentially; does not push
  latest tag to container registry except on tagged releases.
- Sort unmerged read fasta files in order of descending abundance (GH59)

## 1.17.1

- Support for ITS using vsearch alignments
- Handle null seqtab in `get_chim_dropped.R`

## 1.16.1

- Actually implements get_unmerged task
- Add script ``get_dropped_chim.R `` identifying SV dropped as chimeras

## 1.15.2

- Implement downsampling (GH44)
- get_unmerged task saves empty file when there are no merged reads

## 1.15.1

- bugfix: fix column names in ljoin.R (GH40)

## 1.15

- New downsample: [int] argument in params (GH38)
- use dada2 1.18
- output includes files containing unmerged reads (GH11)
- remove native Singularity build file

## 1.14

- Handle mixed reads in either direction (GH37)

## 1.13

- Replace cmalign 16s filtering step with cmsearch (GH13)
- Update covariance model to use version from https://rfam.xfam.org/family/RF00177
