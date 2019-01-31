# BIOS_ASE

# Generate count files for haplotype A and B from geneAE data
sh createCountTables.sh

# Perform binominal test on count data

# Run binomial test on the a/b counts
Rscript ASE_binomial_test/binom_test.R

# Make table with ASE genes (bonf. corrected p-value < 0.05
python ASE_outlier_table/make_outlier_table.py
