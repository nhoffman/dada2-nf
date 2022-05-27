# Change log for dada2-nf

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
