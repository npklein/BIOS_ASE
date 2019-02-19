
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
carriers = {}
info_per_maf = {}
count_per_gene = {}
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
            if gene not in carriers:
                carriers[gene] = {}
            impact = info.split('|')[2]

            MAF = float(line[7].split('MAF=')[1].split(';')[0])
                
            rounded_maf = round(MAF, round_maf_number)
            snp = line[0]+'_'+line[1]

            snp_info = []
            for index, element in enumerate(line):
                if index < 9:
                    continue
                element = element.split(':')[0]
                sample = header_index[index]
                if element == '1|0' or element == '0|1':
                    chr = line[0]
                    pos = line[1]
                    if snp not in carriers[gene]:        
                         carriers[gene][snp] = []
                    carriers[gene][snp].append([header[index], element, MAF, impact, gene in OMIM_genes])

    
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/omim_enrichment/omim_carriers/OMIM_carriers_allIMPACT.hetsOnly.txt','w') as out:
    out.write('gene\tsnp\tsampleName\tgenotype\tgeneSymbol\tENSG\tMAF\timpact\tisOmimGene\n')
    for gene in carriers:
        for snp in carriers[gene]:
            for carrier in carriers[gene][snp]:
                gene_symbol = 'NA'
                if gene in gene_to_symbol:
                    gene_symbol = gene_to_symbol[gene]
                out.write(gene+'\t'+snp+'\t'+carrier[0]+'\t'+carrier[1]+'\t'+gene_symbol+'\t'+gene+'\t'+str(carrier[2]))
                out.write('\t'+str(carrier[3])+'\t'+str(carrier[4])+'\n')
    
