import glob
import sys
from scipy.stats import binom_test
import math


print('Read cumulative counts')
ratio_per_gene_snp = {}
set_of_gene_snp = set([])
ref_alt_per_gene_snp = {}

for f in glob.glob('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.cumulativeVariants.ALLcounts*.chrALL.txt.filtered.txt'):
    print('parsing',f)
    with open(f) as input_file:
        header = input_file.readline()
        for line in input_file:
            line = line.strip().split('\t')
            gene = line[6]
            snp = '_'.join(line[0].split('_')[0:2])
            ratio_per_gene_snp[gene+'_'+snp] = [line[8],line[9],line[10]]
            set_of_gene_snp.add(gene+'_'+snp)
            ref_alt_per_gene_snp[gene+'_'+snp] = [line[0].split('_')[2],line[0].split('_')[3]]


## Below is not necesarry anymore, because it is tested that eqtlgen snp/gene combi is same as for ASE. keep in for now for reference until testing is finished, remove later
# read in gtf info to check if snp affects the gene the SNP is located in
#gene_loc = {}
#/apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf
#with open('Homo_sapiens.GRCh37.75.gtf') as input_file:
#    for line in input_file:
#        if line.startswith('#'):
#            continue
#        line = line.strip().split('\t')
#        if line[2] != 'exon' or line[1] != 'protein_coding':
#            continue
#        chr = line[0]
#        gene = line[8].split('gene_id "')[1].split('"')[0]
#        start = int(line[3])
#        stop = int(line[4])
#        if gene not in gene_loc:
#            gene_loc[gene] = []
#        gene_loc[gene].append([chr, start, stop])


print('read eqtl data')
eqtl_per_gene_snp = {}
#/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithEqtlGen/data/2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved.aseSNPs.snpsInProtCodingExon.txt
with open('2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved.aseSNPs.snpsInProtCodingExon.eqtlGenRank.txt') as input_file:
    head = input_file.readline().strip().split('\t')
    for line in input_file:
        line = line.strip().split('\t')
        gene = line[0]
        fdr = line[2]
        chr = line[3]
        pos = line[4]
        qtl = gene+'_'+chr+'_'+pos
        rs = line[1]
        eqtl_rank = line[8]
        if gene+'_'+chr+'_'+pos not in set_of_gene_snp:
            continue
        # Since we test that the gene-snp combination is present in the ASE data, we dont need to test if the snp is in the gene
#        if gene not in gene_loc:
#            continue
#        snp_in_gene = False
#        for location in  gene_loc[gene]:
#            if location[0] != chr:
#                continue
#            if int(pos) > int(location[0]) and int(pos) < int(location[1]):
#                snp_in_gene = True
#        if not snp_in_gene:
#            continue
        
        snp_type = line[5]
        allele_assessed = line[6]
        zScore = line[7]
        if gene+'_'+chr+'_'+pos in eqtl_per_gene_snp:
            raise RuntimeError('Combination of gene_snp already seen: '+gene+'_'+chr+'_'+pos)

        eqtl_per_gene_snp[gene+'_'+chr+'_'+pos] = [fdr,snp_type,allele_assessed,zScore,rs,eqtl_rank]
        
    


list_of_gene_snps = sorted(list(set_of_gene_snp))
n_snps = len(set_of_gene_snp)
print('done reading')



print('start writing')

with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithEqtlGen/data/ASE_eQTL.binom.all_eQTL.txt','w') as out:
    out.write('snp\tgene\tref\talt')
    out.write('\trefCount'+'\taltCount'+'\tlogFC'+'\tbinom')
    out.write('\tfdr_eqtlGen\tsnp_type\tallele_assessed\tzScore\tzScore_alt\tzScore_alt_swapped\teqtl_rank\n')

    number_of_snps_total = 0
    number_of_snps_diff_genotype = 0
    number_of_snps_common = 0
    
    for gene_snp in list_of_gene_snps:
        number_of_snps_total += 1
        if number_of_snps_total % 10000 == 0:
            print(str(number_of_snps_total)+'/'+str(n_snps))

        if gene_snp not in eqtl_per_gene_snp:
            continue
        eqtl_data = eqtl_per_gene_snp[gene_snp]
        number_of_snps_common += 1
      
        counts = ratio_per_gene_snp[gene_snp]
        ref = ref_alt_per_gene_snp[gene_snp][0]
        alt = ref_alt_per_gene_snp[gene_snp][1]
        
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
        out.write('_'.join(gene_snp.split('_')[1:])+'\t'+gene_snp.split('_')[0]+'\t'+ref+'\t'+alt)
        refCount = counts[0]
        altCount = counts[1]
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
