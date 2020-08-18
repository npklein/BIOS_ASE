
import gzip
snps = {'10':set(['10_45891347','10_45938895']),
        '15':set(['15_77323533']),
        '17':set(['17_1564671']),'22':set(['22_19951742']),
        '4':set(['4_2833299'])}

variant_sample = {'10_45891347':['AC1JL5ACXX-6-12_AD1GWFACXX-2-12'],
                  '10_45938895':['BD2DCDACXX-2-15'],
                  '15_77323533':['BD2CPRACXX-3-15','BD1NW4ACXX-7-20'],
                  '17_1564671':['BD1NYRACXX-3-8'],
                  '22_19951742':['BC1KAVACXX-2-2'],
                  '4_2833299':['AC1JL5ACXX-7-14_BC52YAACXX-6-14_AD1NRAACXX-3-14']}

#with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/phenotypeTables/allPhenoData/all_sample_ids_table.txt') as input_file:
    

def select_snps(vcf, out_postfix):
    for chr in snps:
        print(chr)
        chr_vcf = vcf.replace('REPLACECHR', chr)
        with gzip.open(chr_vcf,'rt') as input_file,open('vcfs/chr'+chr+'_'+out_postfix+'.txt','w') as out:
            index_to_write = {}
            for line in input_file:
                line = line.strip().split('\t')
                if line[0].startswith('#CHROM'):
                    header = line
                    out.write('\t'.join(header[0:9]))
                    for index, element in enumerate(header):
                        for snp in snps[chr]:
                            for sample in variant_sample[snp]:
                                if element == sample:
                                    out.write('\t'+element)
                                    if snp not in index_to_write:
                                        index_to_write[snp] = []
                                    index_to_write[snp].append(index)
                    out.write('\n')    
                elif line[0].startswith('#'):
                    continue
                else:
                    snp = line[0]+'_'+line[1]
                    if snp in snps[chr]:
                        out.write('\t'.join(line[:9]))
                        print(index_to_write)
                        print(snp)
                        print(index_to_write[snp])
                        print('-'*20)
                        for index in index_to_write[snp]:
                            out.write('\t'+line[index])
                        out.write('\n')


VCF_imputed = '/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing_annotation/annotatedWith.snpEff.closest.VEP/BIOS_LLDeep_noRNAeditSites_phASER.snpEff.closest.VEP.removedCODAM.4outliersRemoved.chrREPLACECHR.vcf.gz'
VCF_not_imputed = '/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chrREPLACECHR.filter.GQ20_callRate50.PASSonly.BiallelicSNVsOnly.noDiagnostics.gg.vcf.gz'

select_snps(VCF_imputed, 'afterImputation_subsetted')
select_snps(VCF_not_imputed, 'beforeImputation_subsetted')
