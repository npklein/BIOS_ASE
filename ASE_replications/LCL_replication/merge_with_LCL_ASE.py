import scipy.stats
import math
import glob
import os
from scipy.stats import binom_test


ratio_per_snp = {}
set_of_snp = set([])
ref_alt_per_snp = {}
for f in glob.glob('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.cumulativeVariants.ALLcounts*.chrALL.txt.filtered.txt'):
    print('parsing',f)
    with open(f) as input_file:
        header = input_file.readline()
        for line in input_file:
            line = line.strip().split('\t')
            snp = '_'.join(line[0].split('_')[0:2])
            ratio_per_snp[snp] = [line[8],line[9],line[10]]
            set_of_snp.add(snp)
            ref_alt_per_snp[snp] = [line[0].split('_')[2],line[0].split('_')[3]]

LCL_ase = {}
# LCL_ASE.txt downloaded from https://molgenis56.gcc.rug.nl/
with open('LCL_ASE.txt') as input_file:
    header = input_file.readline().strip().split('\t')
    for line in input_file:
        line = line.strip().split('\t')
        if line[4]+'_'+line[5] not in set_of_snp:
            continue
        LCL_ase[line[4]+'_'+line[5]] = [line[0],line[2], line[7], line[8]]

data_dir = '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/merged_count_data/'
if not os.path.exists(data_dir):
    os.mkdir(data_dir)

list_of_snps = sorted(list(set_of_snp))

with open(data_dir+'/ASE.mergedWithLCL.txt','w') as out:
    out.write('chr\tpos\tase_genotype\tase_ref\tase_alt\tase_assessed_allele')
    out.write('\tase_other_alle\tlogFC\tbinom_pval\tLCL_pval\tLCL_ratio\tsummed_logFC\tASE_samples\n')
    y = 0
    for snp in list_of_snps:
        chr = snp.split('_')[0]
        pos = snp.split('_')[1]
        counts = ratio_per_snp[snp]
        ref = ref_alt_per_snp[snp][0]
        alt = ref_alt_per_snp[snp][1]

        refCount = counts[0]
        altCount = counts[1]
        if int(refCount) == 0 or int(altCount) == 0:
            binom = 'NA'
            logFC = 'NA'
        else:
            binom = binom_test(int(altCount), int(altCount)+int(refCount), p=0.5,alternative='two-sided')
            logFC = math.log(float(altCount)/float(refCount))
        if snp not in LCL_ase:
            continue
        LCL_pval = LCL_ase[snp][0]
        if alt == LCL_ase[snp][2] and ref == LCL_ase[snp][3]:
            LCL_ratio = 1-float(LCL_ase[snp][1])
        elif alt == LCL_ase[snp][3] and ref == LCL_ase[snp][2]:
            LCL_ratio = LCL_ase[snp][1]
        else:
            print('ASE genotype: '+ref+'/'+alt)
            print('LCL genotype: '+LCL_ase[snp][2]+'/'+LCL_ase[snp][3] )
            print('-'*20)
            continue
            #raise RuntimeError('genotype not same')


        out.write(chr+'\t'+pos+'\t'+line[1]+'/'+line[2]+'\t')
        out.write(line[3]+'\t'+line[4]+'\t'+line[2]+'\t'+line[1]+'\t'+str(logFC)+'\t'+str(binom))
        out.write('\t'+LCL_pval+'\t'+str(LCL_ratio)+'\t'+line[5]+'\t'+line[6]+'\n')
        y += 1
    print(y,'snps overlapping')
