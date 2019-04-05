import scipy.stats
import math
import glob
import os

LCL_ase = {}
# LCL_ASE.txt downloaded from https://molgenis56.gcc.rug.nl/
with open('LCL_ASE.txt') as input_file:
    header = input_file.readline().strip().split('\t')
    for line in input_file:
        line = line.strip().split('\t')
        LCL_ase[line[4]+'_'+line[5]] = [line[0],line[2], line[7], line[8]]


data_dir = '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/merged_count_data/'
if not os.path.exists(data_dir+'ASE_counts_merged_with_LCL/'):
    os.mkdir(data_dir+'ASE_counts_merged_with_LCL/')
# TODO: add code to make these files
for f in glob.glob(data_dir+'/counts_summed*txt'):
    print(f)
    name = f.split('.txt')[0].split('per_snp.')[1]
    with open(f) as input_file, open(data_dir+'ASE_counts_merged_with_LCL/ASE.'+name+'.mergedWithLCL.txt','w') as out:
        header = input_file.readline()
        out.write('chr\tpos\tase_genotype\tase_ref\tase_alt\tase_assessed_allele')
        out.write('\tase_other_alle\tlogFC\tbinom_pval\tLCL_pval\tLCL_ratio\tsummed_logFC\tASE_samples\n')
        y = 0
        for line in input_file:
            line = line.strip().split('\t')
            alt = line[2]
            ref = line[1]
            chr = line[0].split('_')[0]
            pos = line[0].split('_')[1]
            snp = '_'.join(line[0].split('_')[:2])
                    
            binom = scipy.stats.binom_test(x=int(line[4]), n=int(line[4])+int(line[3]), p=0.5)
    
    
            if snp not in LCL_ase:
                continue
            else:
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
    
            if int(line[4]) == 0:
                continue
            elif int(line[3]) == 0:
                continue
            else:
                logFC = math.log(float(line[4])/float(line[3]))
    
                    
    
            out.write(chr+'\t'+pos+'\t'+line[1]+'/'+line[2]+'\t')
            out.write(line[3]+'\t'+line[4]+'\t'+line[2]+'\t'+line[1]+'\t'+str(logFC)+'\t'+str(binom))
            out.write('\t'+LCL_pval+'\t'+str(LCL_ratio)+'\t'+line[5]+'\t'+line[6]+'\n')
            y += 1
        print(y,'snps overlapping')
