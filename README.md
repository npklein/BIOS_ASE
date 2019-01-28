# BIOS_ASE

## Generate count files for haplotype A and B from geneAE data
sh createCountTables.sh

# Run binomial test on the a/b counts
Rscript ASE_binomial_test/binom_test.R

# Make table with ASE genes (bonf. corrected p-value < 0.05
python ASE_outlier_table/make_outlier_table.py

# Select samples from the table that don't have more than 1000 ASE genes and are not CODAM
python ASE_outlier_table/select_samples.py


## Create phenotype table
perl createPhenotypeTable.pl

## Merge phenotypedata with the aggregate #genes and #outlier ASE genes from bonferroni corrected logFoldChange matrix
perl createGenesAndOutliersTable.pl


## Select and remove outliers from data, create list with sample IDs to keep
Rscript removeOutliersAndCODAM.R


## Calculate number of ASE genes and outlier genes per sample (using samplelist as obtained from Rscript)
# Only use samples having at least 30X coverage and at least 5X on both haplotypes
# Calculate meand and SD
# Only assess genes for which at least 100 samples show ASE, when a sample is more than 3SD from the mean mark it as outlier for that specific gene

perl createLogFoldTable.pl


## Count per gene the number of homs and hets, output in long format
python allele_count_tables/combine_genes_and_samples.py
