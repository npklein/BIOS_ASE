# BIOS_ASE

<<<<<<< HEAD
# Generate count files for haplotype A and B from geneAE data
sh createCountTables.sh

# Perform binominal test on count data

=======
## Generate count files for haplotype A and B from geneAE data
sh createCountTables.sh

>>>>>>> 2a80f51c302edfb1917d2925ee56cc7cdaab0734
# Run binomial test on the a/b counts
Rscript ASE_binomial_test/binom_test.R

# Make table with ASE genes (bonf. corrected p-value < 0.05
python ASE_outlier_table/make_outlier_table.py

# Select samples from the table that don't have more than 1000 ASE genes and are not CODAM
python ASE_outlier_table/select_samples.py

<<<<<<< HEAD
## Generate count files for haplotype A and B from geneAE data
sh createCountTables.sh


## Perform binominal test on count data
NPK fill in

=======
>>>>>>> 2a80f51c302edfb1917d2925ee56cc7cdaab0734

## Create phenotype table
perl createPhenotypeTable.pl

<<<<<<< HEAD

=======
>>>>>>> 2a80f51c302edfb1917d2925ee56cc7cdaab0734
## Merge phenotypedata with the aggregate #genes and #outlier ASE genes from bonferroni corrected logFoldChange matrix
perl createGenesAndOutliersTable.pl


## Select and remove outliers from data, create list with sample IDs to keep
Rscript removeOutliersAndCODAM.R


## Calculate number of ASE genes and outlier genes per sample (using samplelist as obtained from Rscript)
<<<<<<< HEAD
Only use samples having at least 30X coverage and at least 5X on both haplotypes
Calculate meand and SD
Only assess genes for which at least 100 samples show ASE, when a sample is more than 3SD from the mean mark it as outlier for that specific gene
=======
# Only use samples having at least 30X coverage and at least 5X on both haplotypes
# Calculate meand and SD
# Only assess genes for which at least 100 samples show ASE, when a sample is more than 3SD from the mean mark it as outlier for that specific gene
>>>>>>> 2a80f51c302edfb1917d2925ee56cc7cdaab0734

perl createLogFoldTable.pl


<<<<<<< HEAD
=======
## Count per gene the number of homs and hets, output in long format
python allele_count_tables/combine_genes_and_samples.py



#Minor allele analysis

## Create major/minor allele tables
perl minor_allele_ratio/createTable.pl

## Filter out all lines with only NA values from major/minor tables
perl minor_allele_ratio/filterMatrices.pl

## Annotate table with CADD information
module load HTSlib/1.3.2-foss-2015b
perl minor_allele_ratio/annotateCountsWithCADD.pl

## Create impact category plot, use generated files counts.matrix.m*rAllelle.chrALL.txt.filtered.txt and counts.chr22.addedCADD.txt as input
Rscript minor_allele_ratio/plot_minor_vs_major_20190129.R


#Gene expression analysis and ASE

## Create gene expression table for all samples
perl geneExpressionTables/selectAllSamplesExcludingCODAMand4outliersForAllGenes.pl

>>>>>>> 2a80f51c302edfb1917d2925ee56cc7cdaab0734

