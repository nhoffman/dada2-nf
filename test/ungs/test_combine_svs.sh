#!/bin/sh
python3 ../../bin/combine_svs.py --out corrected_seqtab.csv clusters.uc seqs.fa
diff corrected_seqtab.csv expected_seqtab.csv
rm corrected_seqtab.csv
