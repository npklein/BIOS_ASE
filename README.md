# BIOS_ASE

# NOTE: We use hardcoded paths to locations of files on our cluster, if you want to 
#       do these analysis on your own data you wil have to change the paths in the scripts

## Generate count files for haplotype A and B from geneAE data
sh createCountTables.sh

# Run binomial test on the a/b counts
Rscript ASE_binomial_test/binom_test.R

# Make table with ASE genes (bonf. corrected p-value < 0.05
python ASE_outlier_table/make_outlier_table.py

# Select samples from the table that don't have more than 1000 ASE genes and are not CODAM
python ASE_outlier_table/select_samples.py

## Perform binominal test on count data
Rscript binom_sample_ASE_test.all.R

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


#Carriers per disease/inheritance (fig 2 and 3)
## Overlap the heterozygous SNPs with the OMIM data to know in which OMIM gene the SNP is located
python figure_2_and_3/OMIM_enrichment_hetsOnly.py

##Get the 3 star clinvar variants
Rscript figure_2_and_3/get_clinvar_pathogenics.R

##Calcualte enrichment in disease genes
python figure_2_and_3/enrichment_disease_genes_in_outliers_per_category_hetsOnly.py

##Make figure 2 and 3
Rscript figure_2_and_3/plot_carriers_per_clinvar_hetsOnly.R

# Analysis of AD pathogenic variants (fig 5)
## plot AD genes
Rscript figure_5/plot_AD_genes.R
