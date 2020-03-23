import glob
import gzip
clinvar_info = {}
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/clinvar/clinvar_data/clinvar_20180128.criteriaProvided_multiple_submitters_no_conflicts.txt') as input_file:
    input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        snp = line[1]+'_'+line[2]
        clinstat  = line[15]
        if 'benign' in clinstat.lower():
            clinstat = 'Benign'
        elif 'pathogenic' in clinstat.lower():
            clinstat = 'Pathogenic'
        elif 'uncertain_significance' in clinstat.lower():
            clinstat = 'VUS'
        else:
            raise RuntimeError(clinstat+' not benign, pathogenic, or uncertain significance')
        clinvar_info[snp] = clinstat

def read_isOutlier():
    outliers = {}
    not_outliers = {}
    na = {}
    input_dir = '/groups/umcg-bios/tmp03/projects/outlierGeneASE//logFoldChangeTables/'
    # Input file made by ../createLogFoldChangeTable.pl
    with open(input_dir+'genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.samplesFILTERED.txt') as input_file:
        header = input_file.readline().strip().split('\t')
        for line in input_file:
            line = line.strip().split('\t')
            gene = line[0]
            outliers[gene] = 0
            not_outliers[gene] = 0
            na[gene] = 0
            for index, element in enumerate(line[1:]):
                if element == 'YES':
                    outliers[gene] += 1
                elif element == 'NO':
                    not_outliers[gene] += 1
                elif element == 'NA':
                    na[gene] += 1
                else:
                    raise RuntimeError(element+' should not be one of the options')
    return(outliers, not_outliers, na)
outliers, not_outliers, na = read_isOutlier()



set_of_snps = set()
snp_not_in_clinvar = 0
snp_count = {}
snp_gene = {}
for f in glob.glob('/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing_annotation/annotatedWith.snpEff.closest.VEP/*removedCODAM*vcf.gz'):
    print(f)
    with gzip.open(f,'rt') as input_file:
        for line in input_file:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            snp = line[0]+'_'+line[1]
            gene_info = info.split(';')[6]
            gene = info.split('|')[4]

            if snp in snp_gene_info:
                raise RuntimeError('each SNP should only be in 1 time')
            snp_gene[snp] = gene

            if snp not in clinvar_info:
                snp_not_in_clinvar += 1
            else:
                if snp in snp_count:
                    raise RuntimeError('SNP should only be in once')
                set_of_snps.add(snp)
                snp_count[snp] = {'het':0,'ref':0,'alt':0}
                for element in line[9:]:
                    genotype = element.split(':')[0]
                    if genotype == '1|0' or genotype == '0|1':
                        snp_count[snp]['het'] += 1
                    elif genotype == '1|1':
                        if snp == '13_77574606':
                            print(snp, genotype)
                        snp_count[snp]['alt'] += 1
                    elif genotype == '0|0':
                        snp_count[snp]['ref'] += 1
                    else:
                        raise RuntimeError(genotype+' not expected genotype')

outfile='/groups/umcg-bios/tmp03/projects/outlierGeneASE/clinvar/clinvar_overlapped_with_SNPs.txt'
with open(outfile,'w') as out:
    out.write('snp\tgene\tclinstat\tref\thet\talt\toutlier\tnot_outlier\tna\n')
    for snp in set_of_snps:
        if snp in clinvar_info:
            clinstat = clinvar_info[snp]
        else:
            clinstat = 'not_in_clinvar'
        if snp == '13_77574606':
            print(snp_count[snp])
        out.write(snp+'\t'+gene+'\t'+clinstat+'\t'+str(snp_count[snp]['ref'])+'\t'+str(snp_count[snp]['het'])+'\t'+str(snp_count[snp]['alt']))
        out.wite('\t'+str(outliers[gene])+'\t'+str(not_outliers[gene])+'\t'+str(na[gene])+'\n')

print('writtent to '+outfile)
