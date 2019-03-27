import glob
from scipy.stats import binom_test
import math
import gzip
from _operator import pos

print('read gtex data')
gtex_info = {}
set_of_snps = set([])
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/concordanceGTEx/counts.matrix.AlleleAdded.txt') as input_file:
    header = input_file.readline().strip().split('\t')
    for line in input_file:
        line = line.strip().split('\t')
        snp = '_'.join(line[0].split('_')[:2])
        ref = line[3]
        alt = line[4]
        minor_allele = line[6]
        major_count = line[16]
        minor_count = line[17]
        tissue = line[14]
        if tissue != 'WHLBLD':
            continue
        if snp in gtex_info:
            print(snp)
            raise RuntimeError('SNP multiple times in file')
        gtex_info[snp] = [snp, ref, alt, minor_allele, major_count, minor_count]
        set_of_snps.add(snp)

print('read top eQTL')
top_eqtl = set([])
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithEqtlGen/data/eqtlRank.txt') as input_file:
    header = input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        top_eqtl.add('_'.join(line))

# read in gtf info to check if snp affects the gene the SNP is located in
gene_loc = {}
with open('/apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf') as input_file:
    for line in input_file:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[2] != 'exon' or line[1] != 'protein_coding':
            continue
        chr = line[0]
        gene = line[8].split('gene_id "')[1].split('"')[0]
        start = int(line[3])
        stop = int(line[4])
        if gene not in gene_loc:
            gene_loc[gene] = []
        gene_loc[gene].append([chr, start, stop])
        

ref_alt_per_snp = {}
f = '/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithBios/data/count_data/counts_summed_per_snp.binom.txt'
print(f)

with open(f) as input_file:
    input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        snp = line[0]
        if snp not in ref_alt_per_snp:
            ref_alt_per_snp[snp] = [line[1],line[2]]
        set_of_snps.add(snp)
        if snp not in different_snpCounts:
            different_snpCounts[snp] = {}

        different_snpCounts[snp] = line[1:]


print('start writing')
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithGTEx/ASE_GTex.binom.txt','w') as out:
    out.write('snp\tmajorAllele\tminorAlle')
    out.write('\tmajorCount'+'\tminorCount'+'\tlogFC'+'\tbinom')
    
    out.write('\tfdr\tsnp_type\tallele_assessed\tzScore\tzScore_alt\tzScore_alt_swapped\n')

    n_snps = len(set_of_snps)
    number_of_snps_total = 0
    for snp in set_of_snps:
        number_of_snps_total += 1
        if number_of_snps_total % 10000 == 0:
            print(str(number_of_snps_total)+'/'+str(n_snps))
        eqtl_snp = '_'.join(snp.split('_')[:2])
        if eqtl_snp not in eqtl_per_snp:
            continue 
        eqtl_data = eqtl_per_snp[eqtl_snp]
        gtex_data = gtex_info[snp]  
        gtex_ref = gtex_data[1]
        gtex_alt = gtex_data[2]
        minor_allele = gtex_data[3]
        eqtl_ref = eqtl_data[1].split('/')[0]
        eqtl_alt = eqtl_data[1].split('/')[1]
        allele_assessed = eqtl_data[2]
        if allele_assessed == eqtl_ref:
            zscore_alt = str(-1*float(eqtl_data[3]))
        elif allele_assessed == eqtl_alt:
            zscore_alt = eqtl_data[3]
        else:
            raise RuntimeError('assessed allele not in SNP type')
        
        zscore_swapped = zscore_alt

        if not (eqtl_alt == gtex_alt and eqtl_ref == gtex_ref) or (eqtl_ref == gtex_alt and eqtl_alt == gtex_ref):
            print(eqtl_ref+'/'+eqtl_alt,gtex_ref+'/'+gtex_alt)
            print('Genotypes not the same')
            continue

        if eqtl_ref == minor_allele:
            zscore_swapped = str(-1*float(zscore_swapped))
        elif eqtl_alt != minor_allele:
            print(eqtl_ref+'/'+eqtl_alt, minor_allele)
            raise RuntimeError('should have een caught at earlier exception')
                    
        out.write(snp+'\t'+ref+'\t'+alt)
        major_count = int(gtex_data[4])
        minor_count = int(gtex_data[5])

        if int(major_count) == 0 or int(minor_count) == 0:
            binom = 'NA'
            logFC = 'NA'
        else:
            binom = binom_test(int(minor_count), int(minor_count)+int(major_count), p=0.5,alternative='two-sided')
            logFC = math.log(float(minor_count)/float(major_count))


        
        out.write('\t'+str(major_count)+'\t'+str(minor_count)+'\t'+str(logFC)+'\t'+str(binom))
        out.write('\t'+eqtl_data[0]+'\t'+eqtl_data[1]+'\t'+eqtl_data[2]+'\t')
        out.write(eqtl_data[3]+'\t'+zscore_alt+'\t'+zscore_swapped+'\n')

