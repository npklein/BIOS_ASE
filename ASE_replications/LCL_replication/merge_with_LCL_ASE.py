import scipy.stats
import math
import glob
import os
from scipy.stats import binom_test


allele_freqs = {}
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/annotatedWith.snpEff.closest.VEP/chrALL.AFsFromData.txt') as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        allele_freqs[line[0]] = line[1]

ratio_per_snp = {}
set_of_snp = set([])
ref_alt_per_snp = {}
ASE_snp_name = {}
for f in glob.glob('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.cumulativeVariants.ALLcounts*.chrALL.txt.filtered.txt'):
    print('parsing',f)
    with open(f) as input_file:
        header = input_file.readline()
        for line in input_file:
            line = line.strip().split('\t')
            snp = '_'.join(line[0].split('_')[0:2])
            ASE_snp_name[snp] = line[0]
            ratio_per_snp[snp] = [line[8],line[9],line[10]]
            set_of_snp.add(snp)
            ref_alt_per_snp[snp] = [line[0].split('_')[2],line[0].split('_')[3], line[7]]

LCL_ase = {}
# LCL_ASE.txt downloaded from https://molgenis56.gcc.rug.nl/
with open('LCL_ASE.txt') as input_file:
    header = input_file.readline().strip().split('\t')
    for line in input_file:
        line = line.strip().split('\t')
        if line[4]+'_'+line[5] not in set_of_snp:
            continue
        LCL_ase[line[4]+'_'+line[5]] = [line[0],line[2], line[7], line[8], line[6]]

data_dir = '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/merged_count_data/'
if not os.path.exists(data_dir):
    os.mkdir(data_dir)

list_of_snps = sorted(list(set_of_snp))

with open(data_dir+'/ASE.mergedWithLCL.txt','w') as out:
    out.write('chr\tpos\tLCL_pval\tLCL_ratio\tLCL_ratio_swapped\t')
    out.write('bios_minor_allele\tbios_major_allele\tmajorCount\tminorCount\tLCL_ref\tLCL_alt\t')
    out.write('BIOS_N_samples\tLCL_N_samples\n') 
    y = 0
    for snp in list_of_snps:
        chr = snp.split('_')[0]
        pos = snp.split('_')[1]
        counts = ratio_per_snp[snp]
        ref = ref_alt_per_snp[snp][0]
        alt = ref_alt_per_snp[snp][1]
        BIOS_samples = ref_alt_per_snp[snp][2]

        major = ref
        minor = alt
        if float(allele_freqs[ASE_snp_name[snp]]) > 0.5:
            minor = ref
            major = alt

        refCount = counts[0]
        altCount = counts[1]

        if snp not in LCL_ase:
            continue
        LCL_pval = LCL_ase[snp][0]
        if minor == LCL_ase[snp][2] and major == LCL_ase[snp][3]:
            LCL_ratio_swapped = 1-float(LCL_ase[snp][1])
        elif minor == LCL_ase[snp][3] and major == LCL_ase[snp][2]:
            LCL_ratio_swapped = LCL_ase[snp][1]
        else:
            print('ASE genotype: '+ref+'/'+alt)
            print('LCL genotype: '+LCL_ase[snp][2]+'/'+LCL_ase[snp][3] )
            print('-'*20)
            continue
            #raise RuntimeError('genotype not same')


        out.write(chr+'\t'+pos+'\t'+LCL_pval+'\t'+str(LCL_ase[snp][1])+'\t'+str(LCL_ratio_swapped))
        out.write('\t'+minor+'\t'+major+'\t'+refCount+'\t'+altCount+'\t'+LCL_ase[snp][2]+'\t'+LCL_ase[snp][3])
        out.write('\t'+BIOS_samples+'\t'+LCL_ase[snp][4]+'\n')
        y += 1
    print(y,'snps overlapping')
