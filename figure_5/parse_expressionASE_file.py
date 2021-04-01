
import gzip
import numpy as  np
from scipy import stats
# /groups/umcg-bios/tmp04/projects/copy_from_tmp03/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.samplesFILTERED.txt
outliers = set()
gene_has_an_outlier = set()
outlier_per_gene = {}
with open('genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.samplesFILTERED.txt') as input_file:
    samples = input_file.readline().strip().split('\t')
    for line in input_file:
        line = line.strip().split('\t')
        gene = line[0]
        for index, element in enumerate(line):
            if element == 'YES':
                outliers.add(samples[index]+'_'+gene)
                gene_has_an_outlier.add(gene)
                if gene not in outlier_per_gene:
                    outlier_per_gene[gene] = 0
                outlier_per_gene[gene] += 1


x = 0
gene_seen_last = {}
prev_gene = None
expression = []
expression_sample_seen = set()
outlier_expression = []
sample_list = []
outlier_info = []
with gzip.open('geneExpressionAndMajorMinorAlleleCounts.allGenes.sorted.txt.gz','rt') as input_file,  open('expression_ase_vs_non_ase.txt','w') as out:
    header = input_file.readline().strip().split('\t')
    out.write('gene\toutlier\tpercentile\tmajor_count\tminor_count\tratio\n')
    for line in input_file:
        x += 1
        if x % 1000000 == 0:
            print(x,'lines processed')
        line = line.strip().split('\t')

        # because this file has every sample multiple times for each snp, but we only want each sample 1 time per gene
        if line[7] + '_' + line[5] in expression_sample_seen:
            continue
        expression_sample_seen.add(line[7] + '_' + line[5])

        if line[8] == '':
            continue
        gene = line[5]
        if prev_gene and gene != prev_gene:
            for index, expr in enumerate(outlier_expression):
                out.write(gene+'\t'+sample_list[index] +'\t'+str(stats.percentileofscore(expression, expr)))
                out.write('\t'+outlier_info[index][0]+'\t'+outlier_info[index][1]+'\t'+outlier_info[index][2]+'\n')
            expression = []
            outlier_expression = []
            sample_list = []
            outlier_info = []
            if gene in gene_seen_last:
                print(gene,'alread seen, last pos:',gene_seen_last[gene],'. Current pos:',x)
                raise RuntimeError('Genes should be sorted')
            gene_seen_last[gene] = x
        expression.append(float(line[8]))
        if line[12] != 'NA':
            outlier_expression.append(float(line[8]))
            outlier_info.append([line[10], line[11], line[12]])
        sample_list.append(line[7])
        prev_gene = gene

