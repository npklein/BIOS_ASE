import glob
from scipy.stats import binom_test
import math


different_snpCounts = {}
set_of_snps = set([])
list_of_names = []

rank_of_eqtl = {}
top_eqtl = set([])
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithEqtlGen/data/eqtlRank.txt') as input_file:
    for line in input_file:
        if len(line.strip().split('\t')) == 1:
            continue
        line = line.strip().split('\t')
        top_eqtl.add(line[0]+'_'+line[1])
        rank_of_eqtl[line[0]+'_'+line[1]] = line[2]
        
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


print('read eqtl data')
eqtl_per_snp = {}
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithEqtlGen/data/2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved.aseSNPs.snpsInProtCodingExon.txt') as input_file:
    head = input_file.readline().strip().split('\t')
    for line in input_file:
        line = line.strip().split('\t')
        gene = line[0]
        fdr = line[2]
        chr = line[3]
        pos = line[4]
        qtl = gene+'_'+chr+'_'+pos
        rs = line[1]
        if qtl not in top_eqtl:
            continue
        if gene not in gene_loc:
            continue
        snp_in_gene = False
        for location in  gene_loc[gene]:
            if location[0] != chr:
                continue
            if int(pos) > int(location[0]) and int(pos) < int(location[1]):
                snp_in_gene = True
        if not snp_in_gene:
            continue
        
        snp_type = line[5]
        allele_assessed = line[6]
        zScore = line[7]
        if chr+'_'+pos not in eqtl_per_snp:
            eqtl_per_snp[chr+'_'+pos] = []
        eqtl_per_snp[chr+'_'+pos].append([fdr,snp_type,allele_assessed,zScore,rs, rank_of_eqtl[qtl]])
    

ref_alt_per_snp = {}
f = '/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithBios/data/count_data/counts_summed_per_snp.binom.txt'
print(f)

with open(f) as input_file:
    #['snp', 'ref', 'alt', 'refCount', 'altCount']
    #['15_102516492_C_G', 'C', 'G', '3124', '1272']
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
list_of_snps = sorted(list(set_of_snps))
print('done reading')



print('start writing')
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithEqtlGen/data/ASE_eQTL.topEqtlOnly.binom.all_eQTL.txt','w') as out:
    out.write('snp\tref\talt')
    out.write('\trefCount'+'\taltCount'+'\tlogFC'+'\tbinom')
    
    out.write('\tfdr_eqtlGen\tsnp_type\tallele_assessed\tzScore\tzScore_alt\tzScore_alt_swapped\teqtl_rank\n')

    number_of_snps_total = 0
    number_of_snps_diff_genotype = 0
    number_of_snps_common = 0
    n_snps = len(set_of_snps)
    for snp in set_of_snps:
        number_of_snps_total += 1
        if number_of_snps_total % 10000 == 0:
            print(str(number_of_snps_total)+'/'+str(n_snps))
        eqtl_snp = '_'.join(snp.split('_')[:2])
        if eqtl_snp not in eqtl_per_snp:
            continue 
        eqtl_dataset = eqtl_per_snp[eqtl_snp]
        number_of_snps_common += 1
        for eqtl_data in eqtl_dataset:
        
            counts = different_snpCounts[snp]
            ref = ref_alt_per_snp[snp][0]
            alt = ref_alt_per_snp[snp][1]
            
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

            if eqtl_ref != ref or eqtl_alt != alt:
                if eqtl_ref == alt and eqtl_alt == ref:
                    zscore_swapped = str(-1*float(zscore_swapped))
                else:
                    number_of_snps_diff_genotype += 1
                    continue
            out.write('_'.join(snp.split('_')[:2])+'\t'+ref+'\t'+alt)
            refCount = counts[2]
            altCount = counts[3]
            if int(refCount) == 0 or int(altCount) == 0:
                binom = 'NA'
                logFC = 'NA'
            else:
                binom = binom_test(int(altCount), int(altCount)+int(refCount), p=0.5,alternative='two-sided')
                logFC = math.log(float(altCount)/float(refCount))

            rs = eqtl_data[4]
            
            out.write('\t'+refCount+'\t'+altCount+'\t'+str(logFC)+'\t'+str(binom))
            out.write('\t'+eqtl_data[0]+'\t'+eqtl_data[1]+'\t'+eqtl_data[2]+'\t'+eqtl_data[3]+'\t')
            out.write(zscore_alt+'\t'+zscore_swapped+'\t'+str(eqtl_data[5])+'\n')

print('total number of snps:',  number_of_snps_total)
print('common SNPs snps:',  number_of_snps_common)
print(str(number_of_snps_diff_genotype)+'/'+str(number_of_snps_common)+' genotypes not the same')
