# Requirements

Python libraries: dadi, numpy, scipy
R libraries: bbmle

# Input file

The input file should be a file containing allele count configurations for each of the populations in the 3-population tree (see test_input.txt for an example file). Note that here, unlike in input files for 3P-CLR, the input file should also include sites where the derived allele is extinct or fixed in population C.

# Example line

Rscript CalcDriftsDadi3P.R test_input.txt test_output.txt
