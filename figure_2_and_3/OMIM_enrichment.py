
import random
import glob
import gzip
import os
from random import shuffle


round_maf_number = 5
OMIM_genes = set([])
gene_to_symbol = {}
# Genes downloaded from OMIM
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/omim_enrichment/OMIM.20171220.ensembleGenes.txt') as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        OMIM_genes.add(line[0])
        gene_to_symbol[line[0]] = line[1]


annotation_dir = '/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing_annotation/annotatedWith.snpEff.closest.VEP/'

# from the annotation table extract which SNPs are located within each gene
sample_genotypes = {}
snp_info = {}
set_of_samples = set([])
for vcf in glob.glob(annotation_dir+'BIOS_LLDeep_noRNAeditSites_phASER.snpEff.closest.VEP.removedCODAM.4outliersRemoved.chr*.vcf.gz'):
    print(vcf)
    header_index = {}
    with gzip.open(vcf) as input_file:
        for line in input_file:
            line = line.decode('utf-8')
            if line.startswith('#CHR'):
                header = line.strip().split('\t')
                for index, element in enumerate(header):
                    header_index[index] = element
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            info = line[7]
            gene_info = info.split(';')[6]
            gene = info.split('|')[4]

            impact = info.split('|')[2]

            MAF = float(line[7].split('MAF=')[1].split(';')[0])
                
            snp = line[0]+'_'+line[1]

            if gene not in sample_genotypes:
                sample_genotypes[gene] = {}
            if snp not in sample_genotypes[gene]:
                sample_genotypes[gene][snp] = {}

            if snp in snp_info:
                assert(snp_info[snp] == [MAF, impact, gene in OMIM_genes])
            else:
                snp_info[snp] = [MAF, impact, gene in OMIM_genes]
    
            for index, element in enumerate(line):
                if index < 9:
                    continue
                element = element.split(':')[0]
                sample = header_index[index]
                set_of_samples.add(sample)
                if element == '0|0':
                    dosage = ''
                elif element == '1|0' or element == '0|1':
                    dosage = '1'
                elif element == '1|1':
                    dosage = '2'
                else:
                    raise RuntimeError('Genotype not expected: '+element)
                sample_genotypes[gene][snp][sample] = dosage

samples = sorted(set_of_samples)
print('Got',len(samples),'samples')
print('start writing')

outfile = '/groups/umcg-bios/tmp03/projects/outlierGeneASE/omim_enrichment/omim_carriers/OMIM_carriers_allIMPACT.txt'
with open(outfile,'w') as out:
    out.write('gene\tsnp\tMAF\timpact\tisOmimGene')
    for sample in samples:
        out.write('\t'+sample)
    out.write('\n')
    for gene in sorted(sample_genotypes):
        for snp in sorted(sample_genotypes[gene]):
            out.write(gene+'\t'+snp+'\t'+str(snp_info[snp][0])+'\t'+str(snp_info[snp][1])+'\t'+str(snp_info[snp][2]))
            for sample in samples:
                out.write('\t'+sample_genotypes[gene][snp][sample])
            out.write('\n')    

print('written to '+outfile)
