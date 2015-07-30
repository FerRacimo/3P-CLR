# 3P-CLR

3P-CLR is a method to detect positive selection in a three-population tree using patterns of linked allele frequency differentiation. It can distinguish selective events that occurred in the common ancestral population of the two most closely-related populations from those events that occurred in each of the populations after the split from each other.

To learn more about the method, you can check out the following preprint:
http://biorxiv.org/content/early/2015/05/12/017566.1

# To compile:

cd src

g++ main.cpp ASFfunc.cpp hcubature.c -std=c++0x -o threepclr

cd ..

# Example line for a test file

src/threepclr "test/test_input.txt" "test/test_output.txt" 20 100 0.0025 0.086,0.102,0.17 NA

# Explanation of example line

Input: test/test_input.txt

Output: test/test_output.txt

Space between analyzed SNPs: 20

Number of SNPs sampled per window: 100

Window length in Morgans: 0.0025

Drift A: 0.086

Drift B: 0.102

Drift C (Inner branch drift + outgroup branch drift): 0.17

Note that all drifts are equal to time (in generations) divided by 2*Ne

Location of file with LD weights (if available, otherwise NA): NA


# Input columns 

(The header can be modified as needed)

column 1: Chromosome

column 2: Physical position of SNP

column 3: Genetic position of SNP

column 4: Number of derived alleles in population A

column 5: Number of sampled alleles in population A

column 6: Number of derived alleles in population B

column 7: Number of sampled alleles in population B

column 8: Number of derived alleles in outgroup population C

column 9: Number of sampled alleles in outgroup population C

# Notes: 

When running 3P-CLR, make sure that all SNPs are segregating in the outgroup population (they can be fixed derived, segregating or fixed ancestral in the two test populations, so long as they are segregating in the outgroup).

Also, make sure to use a single chromosome in each run (i.e. do not combine data from different chromosomes into the same input file).

To obtain positions in genetic distance units, one can run programs like LDhat: https://github.com/auton1/LDhat
It is best if the program is run only in the outgroup population.

# Drift values

3P-CLR requires the values of the drift times of the 3-population tree to be entered as input. Drift times are roughly equal to time in generations divided by twice the haploid effective population size. The program requires 3 drift times:

1) The drift time leading from the common ancestor of A & B to node A (Drift A)

2) The drift time leading from the common ancestor of A & B to node B (Drift B)

3) The drift time leading from the common ancestor of A, B & C to node C, plus the drift time leading from the common ancestor of A, B & C to the common ancestor of A & B (Drift C)

These drift times can be accurately calculated by applying programs like MixMapper (Lipson et al. 2013) to genome-wide data, assuming selection is not prevalent enough in the genome to affect these estimates. Alternatively, we provide an R script that uses dadi (Gutenkunst et al. 2009) and a likelihood optimization algorithm (L-BFGS-B) to obtain approximations to these drift times. The program is located in the folder CalcDrifts3P.

# Output columns

Chr: Chromosome

PhysPos: Physical position of central (candidate beneficial) SNP

GenPos: Genetic position of central SNP

PhysStart: Physical left-most SNP of window

PhysEnd: Physical right-most SNP of window

GenStart: Genetic left-most SNP of window

GenEnd: Genetic right-most SNP of window

3PCLR.Anc: likelihood ratio score for 3P-CLR testing for selection in the ancestral population of A and B

3PCLR.Anc.S: maximum likelihood selection coefficient for 3PCLR.Anc

3PCLR.A: likelihood ratio score for 3P-CLR testing for selection in population A

3PCLR.A.S: maximum likelihood selection coefficient for 3PCLR.A

3PCLR.B: likelihood ratio score for 3P-CLR testing for selection in population B

3PCLR.B.S: maximum likelihood selection coefficient for 3PCLR.B

XPCLRAC: likelihood ratio score for XP-CLR testing for selection in A, using C as outgroup and ignoring B

XPCLRAC.S: maximum likelihood selection coefficient for XPCLRAC

XPCLRBC: likelihood ratio score for XP-CLR testing for selection in B, using C as outgroup and ignoring A

XPCLRBC.S: maximum likelihood selection coefficient for XPCLRBC
