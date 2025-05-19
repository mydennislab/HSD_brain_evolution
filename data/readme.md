# Tajima's D values

Genome-wide Tajima's D values calculated using unrelated 1000 Genomes Project individuals (n=2,504) from five continental ancestries (AFR: African; EUR: European; EAS: East Asian; SAS: South Asian; AMR: American).

Tajima's D values were calculated using `vcftools` in 25-kbp sliding windows. Files include only regions with at least 5 SNPs and >50% accessible base pairs according to the short-read accessibility mask calculated for T2T-CHM13 (v1.0) genome.

Columns legend:
1. Chromosome
2. Start
3. End
4. SNPs
5. Tajima's D
6. SD98 bp
7. Accessible bp
8. Overlapping SD98 genes
