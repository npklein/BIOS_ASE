import glob
import gzip
annotation_dir = '/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing_annotation/annotatedWith.snpEff.closest.VEP/'
def read_gene_info():
    # from the annotation table extract which SNPs are located within each gene
    gene_info = {}
    for table in glob.glob(annotation_dir+'BIOS_LLDeep_noRNAeditSites_phASER.snpEff.closest.VEP.chr*.annotation.table'):
        with open(table) as input_file:
            for line in input_file:
                line = line.strip().split('\t')
                gene = line[11]
                chr = line[0]
                pos = line[1]
                if gene not in gene_info:
                    gene_info[gene] = []
                gene_info[gene].append(chr+'_'+pos)
    return(gene_info)

def read_gene_outliers(foldChange):
    # from the outlier table, extract for each gene which samples are outliers and which are not outliers
    gene_outlier = {}
    # this file is output of ../allele_count_tables/combine_genes_and_samples.py
    with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.samplesFILTERED.txt') as input_file:
        header = input_file.readline().strip().split('\t')
        header_index = {}
        for index, element in enumerate(header):
            header_index[index] = element
        for line in input_file:
            line = line.strip().split('\t')
            gene = line[0]
            gene_outlier[gene] = {'outlier':[],'not_outlier':[]}
            for index, element in enumerate(line):
                # if the sample does not have a fold change measurement (either cause nothing measured or cause Inf/-Inf
#                if header_index[index] not in foldChange[gene]:
#                    continue
                if element == 'NO':
                    gene_outlier[gene]['not_outlier'].append(header_index[index])
                elif element == 'YES':
                    gene_outlier[gene]['outlier'].append(header_index[index])
    return(gene_outlier)

def read_vcf():
    # from the VCF, extract for each SNP which samples are 0 (hom ref), 1 (het) or 2 (hom alt). Save as this dosage
    dosage_per_sample_per_snp = {}
    snp_info = {}
    #out = open('snp_info.txt','w')
    #out.write('snp\t0\t1\t2\n')
    for f in glob.glob(annotation_dir+'BIOS_LLDeep_noRNAeditSites_phASER.snpEff.closest.VEP.chr*.vcf.gz'):
        print(f)
        with gzip.open(f, 'rb') as f:
            header_index = {}
            for line in f:
                line = line.decode('ascii')
                if line.startswith('#CHR'):
                    line = line.strip().split('\t')
                    for index, element in enumerate(line):
                        header_index[index] = element
                    continue
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                chrom = line[0]
                pos = line[1]
                snp = chrom+'_'+pos
                ref = line[3]
                alt = line[4]
                #out.write(snp)
                if snp not in dosage_per_sample_per_snp:
                    dosage_per_sample_per_snp[snp] = {}
                summed_dosage = 0
                zeros = []
                ones = []
                twos = []
                for index, element in enumerate(line):
                    if index < 9:
                        continue
                    genotype = element.split(':')[0]
                    if genotype == '1|0' or genotype == '0|1':
                        dosage = 1
                        ones.append(header_index[index])
                    elif genotype == '0|0':
                        dosage = 0
                        zeros.append(header_index[index])
                    elif genotype == '1|1':
                        dosage = 2
                        twos.append(header_index[index])

                    summed_dosage += dosage
                    dosage_per_sample_per_snp[snp][header_index[index]] = dosage
                #out.write('\t'+','.join(zeros)+'\t'+','.join(ones)+'\t'+','.join(twos)+'\n')
                AF = summed_dosage / float((len(line)-9) * 2)
                snp_info[snp] = {'ref':ref, 'alt':alt,'AF':AF}
    #out.close()
    return(dosage_per_sample_per_snp, snp_info)

def read_fc():
    # per gene extract which genes have a log fold change (so not NA, Inf or -Inf)
    foldChange = {}
    with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.samplesFILTERED.values.txt') as input_file:
        header = input_file.readline().strip().split('\t')
        header_index = {}
        for index, element in enumerate(header):
            header_index[index] = element
        for line in input_file:
            line = line.strip().split('\t')
            gene = line[0]
            foldChange[gene] = set([])
            for index, element in enumerate(line):
                try:
                    float(element)
                except ValueError:
                    continue
                foldChange[gene].add(header_index[index])                

    return(foldChange)

#print('read log FC')
#foldChange = read_fc()
print('read gene inf')
gene_info = read_gene_info()
print('read vcf')
dosage_per_sample_per_snp, snp_info = read_vcf()
#gene_outlier = read_gene_outliers(foldChange)
print('read outliers')
gene_outlier = read_gene_outliers(None)


gene_outlier_sum = {}

with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/infoTables/alleleCountPerGroupPerGene.outliers.binom.Bonferroni.txt','w') as out:
    out.write('gene\tsnp\tAF\tref\talt\t0\t1\t2\ttotal\n')
    # per gene, sum for the outlier and not_outlier group the allele counts. Divide by number of sampls (that have logFC count) to get normalized number of alleles per gene
    for gene in gene_outlier:
        for g in gene.split(';'):
            if g not in gene_info:
                print('gene missing from gene_info:',g)
                continue
            for snp in gene_info[g]:
                out.write(gene+'\t'+snp+'\t'+str(snp_info[snp]['AF'])+'\t'+snp_info[snp]['ref']+'\t'+snp_info[snp]['alt'])
                zeros = 0
                ones = 0
                twos  = 0
                total = 0
                for sample in gene_outlier[gene]['outlier']:
                    dosage =  int(dosage_per_sample_per_snp[snp][sample])
                    if dosage == 0:
                        zeros += 1
                    if dosage == 1:
                        ones += 1
                    if dosage == 2:
                        twos += 1
                    total += 1
                out.write('\t'+str(zeros)+'\t'+str(ones)+'\t'+str(twos)+'\t'+str(total)+'\n')


with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/infoTables/alleleCountPerGroupPerGene.not_outliers.binom.Bonferroni.txt','w') as out:
    out.write('gene\tsnp\tAF\tref\talt\t0\t1\t2\ttotal\n')
    # per gene, sum for the outlier and not_outlier group the allele counts. Divide by number of sampls (that have logFC count) to get normalized number of alleles per gene
    for gene in gene_outlier:
        for g in gene.split(';'):
            if g not in gene_info:
                print('gene missing from gene_info:',g)
                continue
            for snp in gene_info[g]:
                out.write(gene+'\t'+snp+'\t'+str(snp_info[snp]['AF'])+'\t'+snp_info[snp]['ref']+'\t'+snp_info[snp]['alt'])
                zeros = 0
                ones = 0
                twos = 0
                total = 0
                for sample in gene_outlier[gene]['not_outlier']:
                    dosage =  int(dosage_per_sample_per_snp[snp][sample])
                    if dosage == 0:
                        zeros += 1
                    if dosage == 1:
                        ones += 1
                    if dosage == 2:
                        twos += 1
                    total += 1
                out.write('\t'+str(zeros)+'\t'+str(ones)+'\t'+str(twos)+'\t'+str(total)+'\n')




