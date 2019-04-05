# BIOS_ASE

## NOTE: We use hardcoded paths to locations of files on our cluster, if you want to do these analysis on your own data you will have to change the paths in the scripts

### Generate count files for haplotype A and B from geneAE data
`sh createCountTables.sh`

### Run binomial test on the a/b counts
`Rscript ASE_binomial_test/binom_test.R`

### Make table with ASE genes (FDR < 0.05)
`python ASE_outlier_table/make_outlier_table.py`

### Create phenotype table
`perl createPhenotypeTable.pl`

### Merge phenotype table with stats
Merge phenotypedata with the aggregate #genes and #outlier ASE genes from FDR < 0.05 logFoldChange matrix

`perl createGenesAndOutliersTable.pl`

### Remove outliers
Select and remove outliers from data, create list with sample IDs to keep

`Rscript removeOutliersAndCODAM.R`

<br><br>
<br><br>
<br><br>


# ASE outlier (genes showing strong ASE) detection

### Remove outliers
Select samples from the table that don't have more than 1000 ASE genes and are not CODAM

`python ASE_outlier_table/select_samples.py`

### Calculate number of ASE genes and outlier genes per sample (using samplelist as obtained from Rscript)
- Only use samples having at least 30X coverage and at least 5X on both haplotypes
- Calculate meand and SD
- Only assess genes for which at least 100 samples show ASE, when a sample is more than 3SD from the mean mark it as outlier for that specific gene
- uses input from figure_2_and_3/select_samples_from_outlierTable.py
`perl createLogFoldChangeTable.pl`

### Count per gene the number of homs and hets, output in long format
`python allele_count_tables/combine_genes_and_samples.py`

<br><br>
<br><br>
<br><br>


# Minor allele analysis

### Create major/minor allele tables
`perl minor_allele_ratio/createTable.pl`

### Filter out all lines with only NA values from major/minor tables
`perl minor_allele_ratio/filterMatrices.pl`

### Annotate table with CADD information
`module load HTSlib/1.3.2-foss-2015b`

`perl minor_allele_ratio/annotateCountsWithCADD.pl`

### Create impact category plot
Use generated files as input
    counts.matrix.m*rAllelle.chrALL.txt.filtered.txt  (from minor_allele_ratio/minor_allele_ratio/filterMatrices.pl)
    and counts.chr22.addedCADD.txt (from minor_allele_ratio/annotateCountsWithCADD.pl)

`Rscript minor_allele_ratio/plot_minor_vs_major_20190129.R`

<br><br>
<br><br>
<br><br>

# Allele frequency comparison
## Get the allele frequencies from the RNAseq based genotypes
TODO: we compare gonl wgs allele freq with rnaseq genotype allele freq. this gets allele freq for
all samples
`python figure_1/get_allele_frequency.py`

# Carriers per disease/inheritance (fig 2 and 3)

### Overlap the heterozygous SNPs with the OMIM data to know in which OMIM gene the SNP is located
`python figure_2_and_3/OMIM_enrichment_hetsOnly.py`

### Get the 3 star clinvar variants
`Rscript figure_2_and_3/get_clinvar_pathogenics.R`

### Calcualte enrichment in disease genes
`python figure_2_and_3/enrichment_disease_genes_in_outliers_per_category_hetsOnly.py`

### Make figure 2 and 3
`Rscript figure_2_and_3/plot_carriers_per_clinvar_hetsOnly.R`

<br><br>
<br><br>
<br><br>


# Gene expression analysis and ASE

### Create gene expression table for all samples
`perl geneExpressionTables/selectAllSamplesExcludingCODAMand4outliersForAllGenes.pl`

<br><br>
<br><br>
<br><br>


# Combined gene expression and ASE analysis and plots

### Filter count lists
This is a manual step, can otherwise be done using awk for example.

### Create files used as input for figure 5
`perl figure_5/createASEandGeneExpressionTable/createGeneExpressionAndMinorAlleleRatioTables.ListInput.pl`

<br><br>
<br><br>
<br><br>


# CSV files for ASE-browser
There are 3 tables needed to populate the database
- ase_ase
- ase_sampleAse
- ase_genes

### Perform binominal test on individual count data
`Rscript ASE_binomial_test/binom_sample_ASE_test.all.R`

### Create files for tables
`perl ASEbrowserplots/createASEbrowserTablesCsv.pl`

### Create table including all counts
`perl ASEbrowserplots/createSampleAseEntityWithAllCounts.pl`

### Run binomial tests on ase and sampleAse table
`Rscript ASE_binomial_test/binom_snp_aggregate_test.R`

`Rscript ASE_binomial_test/binom_sample_ASE_test.R`

### Split ase_samlpeASE table in smaller chunks, this to produce plots from own laptop (issues with graphial R libraries on cluster)
`perl ASE_binomial_test/splitAse_sampleAseTable.pl`

### Run Rscript to produce plots
`Rscript ASE_binomial_test/manuscript_ASEbrowserPlots.R`

<br><br>
<br><br>
<br><br>


# Plot concordance GoNL DNA vs RNA
`Rscript concordance_GoNL_DNA_vs_RNA.R`

<br><br>
<br><br>
<br><br>


# Pathogenic alleles stratified per disease database

### Create pathogenic allele counts for all variants

`perl pathogenicAlleles/integrateVariantInformationWithAlleleCounts.V2.pl`

`perl pathogenicAlleles/integrateVariantInformationWithAlleleCounts.V3.pl`

### Extract only variants which are present in specific database of interest

`perl pathogenicAlleles/extractOMIMgenesFromAlleleCounts.pl`

`perl pathogenicAlleles/extractOMIMgenesFromAlleleCounts.V3.pl`

`perl pathogenicAlleles/extractCGDgenesFromAlleleCounts.pl`

`perl pathogenicAlleles/extractCGDgenesFromAlleleCounts.V3.pl`

`perl pathogenicAlleles/extractDDG2PgenesFromAlleleCounts.pl`

`perl pathogenicAlleles/extractDDG2PgenesFromAlleleCounts.V3.pl`

### Extract p-values and create plots

`Rscript pathogenicAlleles/enrichment_alt_alleles_per_impact_category.R`

`Rscript pathogenicAlleles/proportion_alt_alleles_per_variant_impact_category.R`

<br><br>
<br><br>
<br><br>


# Analysis of AD pathogenic variants (fig 5)

### plot AD genes

`Rscript figure_5/plot_AD_genes.R`

<br><br>
<br><br>
<br><br>


# Comparison with other ASE and eQTL datasets
## LCL ASE (from https://www.ncbi.nlm.nih.gov/pubmed/25954321, https://molgenis56.target.rug.nl/)
Merge the allelic counts (that are in batches) per sample to create for each chr. one file per sample

`bash ASE_replications/LCL_replication/mergeAllelicCountsPerSample.sh`

### Sum the hapA and hapB counts per SNP and calculate the log fold change over the summed counts

`python ASE_replications/LCL_replication/merge_allelic_counts_per_snp.py`

### Merge our ASE results from merge_allelic_counts_per_snp.py with those of https://www.ncbi.nlm.nih.gov/pubmed/25954321

`python ASE_replications/LCL_replication/merge_with_LCL_ASE.py`

### Plot the single SNP concordance (sup figure ?) and output numbers and p-values used in the manuscript

`Rscript ASE_replications/LCL_replication/plot_single_snp_concordance.R`

<br><br>
<br><br>
<br><br>


## eQTLgen eQTLs (from https://www.biorxiv.org/content/10.1101/447367v1, eqtlgen.org)
### Merge our ASE results from merge_allelic_counts_per_snp.py with the eQTLs from eQTLgen

`python ASE_replications/merge_snpCounts_with_eQTLgen.py`

### plot the concordance between eqtlGen eQTLs and our ASE results

`Rscript ASE_replications/plot_snp_concorance.R`

<br><br>
<br><br>
<br><br>

## Comparison vs GTEx

### Annotate counts with additional information
`perl createAnnotationTableNew.pl`

### Create major/minor table for non-ASE samples
`perl minor_allele_ratio/createTableNonASEsamples.pl`

### Filter table, remove lines containing only NA's
`perl filterMatrices.NonASEsamples.pl`

### Create cumulative count matrices
`perl createCountMatricesCumulativeVariants.pl --postfix nonASEsamples.chrALL.txt.filtered.txt`
`perl createCountMatricesCumulativeVariants.pl --postfix chrALL.txt.filtered.txt`

### Create table including GTEx and our counts/ratios
`perl createTables.AlleleAdded.pl`

### Concordance observed ASE vs GTEx
This step uses the output from createTables.AlleleAdded.pl as input

`Rscript concordance_GTEx_ASE/test_concordance_GTEx_ASE.R`


<br><br>
<br><br>
<br><br>

# Analysis of AD pathogenic variants (fig 5)

### plot AD genes

`Rscript figure_5/plot_AD_genes.R`

<br><br>
<br><br>
<br><br>

# Rare variant enrichment (sup fig 9)

### Calculate enrichment of ASE variants per MAF bin
Take output from allele_count_tables/combine_genes_and_samples.py and calculate for different maf bins
the number of genes that show ASE or not ASE for all the variants in that MAF bin.

`python sup_figure_9/rare_variant_enrichment.py`

### Plot the enrichment

`Rscript sup_figure_9/plot_stratification.R`

<br><br>
<br><br>
<br><br>

# Lipid level correlations
## Correlate lipid level measurements with absolute log fold change
`Rscript logFC_correlations/correlate_logFC.R`

