# Take the merged allelic count data from mergeAllelicCountsPerSample.sh, sum ref/alt counts over all samples
# and then calc. binom p-value and log fold change over those 

import glob
import scipy.stats
import math
import sys


allelic_counts_dir = '/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/readbackedPhasing/allelic_counts_mergedPerSample/'
outdir = '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/merged_count_data/'

# These are all the samples after filtering, have to select those also from output of this script
samples_to_keep_file = '/groups/umcg-bios/tmp03/projects/outlierGeneASE/samples_NOUTLIERS1000.depthFiltered.bonferroni.txt'
with open(samples_to_keep_file) as input_file:
    samples_to_keep = input_file.read().split('\n')        


count_per_snp = {}
samples_per_snp = {}

# There are several filter levels that can be applied on the input data, have the following:
# binom: Only sum ref/alt for samples where the  SNP has a (binom.test) p-value < 0.05
# cov30: ref+alt >= 30
# cov20: ref+alt >= 20 and ref > 0 and alt > 0 (used in manuscript)
# noHighImbalance: imbalance between ref/alt not > 95%
# cov30.noHighImbalance: cov30 and noHighImbalance
# noFilter: not setting any of these filters
# In the manuscript the results shown are for cov20, so other filter options are commented out.
filter_levels = ['cov20'] #'noFilter', #'binom','cov30','noHighImbalance','cov30.noHighImbalance', 'noFilter']

for f_level in filter_levels:
    count_per_snp[f_level] = {}
    samples_per_snp[f_level] = {}

def update_counts(dict_to_update, f_level, snp, refCount, altCount):
    logFC = math.log(float(altCount)/float(refCount), 2)
    dict_to_update[f_level][snp]['ref'] += refCount
    dict_to_update[f_level][snp]['alt'] += altCount
    dict_to_update[f_level][snp]['n_samples'] += 1
    dict_to_update[f_level][snp]['logFC'] += logFC
    return dict_to_update

# The files in allelic_counts_dir are out
files = glob.glob(allelic_counts_dir+'/chr*/*txt')
n_files = len(files)
x = 0


# Below is for temporary tst file, remove code later
out = open('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/tmp/sample_snp_pass_filter.txt','w')
out.write('sample\tSNP\tcov20_min1readPerAllele\n')
for f in files:
    x += 1
    if x % 100 == 0:
        print(str(x)+'/'+str(n_files))
        sys.stdout.flush()
    with open(f) as input_file:
        header = input_file.readline().split('\t')
        for line in input_file:
            line = line.strip().split('\t')
            snp,ref,alt,refCount,altCount = line[2], line[3], line[4], int(line[5]), int(line[6])
            sample = f.split('phASER.')[1].split('.chr')[0]
            if sample not in samples_to_keep:
                continue
            for f_level in filter_levels:
                if snp not in count_per_snp[f_level]:
                    count_per_snp[f_level][snp] = {'ref':0, 'alt':0, 'logFC':0,'n_samples':0}

            if snp.split('_')[2] != ref or snp.split('_')[3] != alt:
                raise RuntimeError('genotypes not as in snp name')
            if refCount == 0 or altCount == 0:
                continue

            if 'cov20' in filter_levels:
                if refCount + altCount >= 20 and refCount > 0 and altCount > 0:
                    count_per_snp['cov20'] = update_counts(count_per_snp, 'cov20', snp, refCount, altCount)['cov20']
                    out.write(sample+'\t'+snp+'\tTRUE\n')
                else:
                    out.write(sample+'\t'+snp+'\tFALSE\n')
            if 'noFilter' in filter_levels:
                count_per_snp['noFilter'] = update_counts(count_per_snp, 'noFilter', snp, refCount, altCount)['noFilter']
        
            if 'cov30' in filter_levels:
                if refCount + altCount >= 30:
                    count_per_snp['cov30'] = update_counts(count_per_snp, 'cov30', snp, refCount, altCount)['cov30']
        
            if 'noHighImbalance' in filter_levels:
                if refCount/(refCount+altCount) > 0.95 or refCount/(refCount+altCount) < 0.05:
                    count_per_snp['noHighImbalance'] = update_counts(count_per_snp, 'noHighImbalance', snp, refCount, altCount)['noHighImbalance']
        
            if 'binom' in filter_levels:
                if scipy.stats.binom_test(x=refCount, n=refCount+altCount, p=0.5) > 0.05:
                    count_per_snp['binom'] = update_counts(count_per_snp, 'binom', snp, refCount, altCount)['binom']
        
            if 'cov30.noHighImbalance' in filter_levels:
                if refCount + altCount >= 30 and (refCount/(refCount+altCount) > 0.95 or refCount/(refCount+altCount) < 0.05):
                    count_per_snp['cov30.noHighImbalance'] = update_counts(count_per_snp, 'cov30.noHighImbalance', snp, refCount, altCount)['cov30.noHighImbalance']

for filter in filter_levels:
    with open(outdir+'/counts_summed_per_snp.'+filter+'.txt','w') as out:
        out.write('snp\tref\talt\tsummedRefCount\tsummedAltCount\tsummedLogFC\tn_samples\n')
        counts = count_per_snp[filter]
        for snp in counts:
            ref = snp.split('_')[2]
            alt = snp.split('_')[3]
            out.write(snp+'\t'+ref+'\t'+alt+'\t'+str(counts[snp]['ref'])+'\t'+str(counts[snp]['alt'])+'\t'+str(counts[snp]['logFC'])+'\t'+str(counts[snp]['n_samples'])+'\n')
