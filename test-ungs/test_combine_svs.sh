#!/bin/sh
python3 ../bin/combine_svs.py vsearch_out.txt weights.csv --corrected_weights corrected_weights.csv
diff corrected_weights.csv expected_corrected_weights.csv
rm corrected_weights.csv