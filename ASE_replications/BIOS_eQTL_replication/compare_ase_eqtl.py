

import gzip

# ['3.27167E-310', 'rs1131017', '12', '56435929', 'ENSG00000197728', '12', '56435637', 'cis', 'C/G', 'C', '76.5071443', '8', '-', '3746', '-', '-', 'RPS26', '-', '1.1142535 (0.014564)', '-', '-', '0.0']
qtl = {}
dir = "/groups/umcg-bios/tmp03/projects/2020-03-05-ASE-BIOS-eqtl-comparison/"
with gzip.open(dir+'BIOS_eqtls/eQTLsFDR0.05-ProbeLevel.txt.gz','rt') as input_file:
    header = input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        qtl[line[2]+'_'+line[3]] = [line[8], line[9], line[10], line[-1]]

allele_freqs = {}
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/annotatedWith.snpEff.closest.VEP/chrALL.AFsFromData.txt') as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        allele_freqs[line[0]] = line[1]

#10_100176104_A_G    10    100176104    A    G    HPS1    ENSG00000107521    232    42205    56252    98457    0.571335710005383
with open(dir+'data/counts.matrix.cumulativeVariants.ALLcounts.chrALL.txt.filtered.txt') as input_file, open(dir+'ase_and_eqtl.txt','w') as out:
    header = input_file.readline().strip('\n')
    out.write(header+'\tzscore_qtl\tzscoreSwapped_qtl\tqtl_fdr\tqtl_genotype\tqtl_assessed_allele\tmajor-allele\tminor-allele\n')
    print(header.split('\t'))
    for line in input_file:
        #['A/G', 'G', '17.4462413', '0.0']
        qtl_name = line.split('_')[0]+'_'+line.split('_')[1]
        if qtl_name not in qtl:
            continue
        this_qtl = qtl[qtl_name]
        qtl_genotype = this_qtl[0]
        qtl_assessed_allele = this_qtl[1]
        zscore = this_qtl[2]
        ref = line.split('\t')[3]
        alt = line.split('\t')[4]

        major = ref
        minor = alt
        if float(allele_freqs[line.split('\t')[0]]) > 0.5:
            minor = ref
            major = alt


        zscoreSwapped = zscore
#        if ref+'/'+alt != qtl_genotype:
#            zscoreSwapped = str(float(zscoreSwapped)*-1)
            


        if qtl_genotype != ref+'/'+alt and qtl_genotype != alt+'/'+ref:
            print('qtl not same, skip '+line.split('\t')[0])
            continue


        if qtl_assessed_allele != minor:
            zscoreSwapped = str(float(zscoreSwapped)*-1)

        if line.split('\t')[0] == '10_14561500_G_A':
            print(line.strip().split('\t'))
            print(this_qtl)
        out.write(line.strip('\n')+'\t'+zscore+'\t'+zscoreSwapped+'\t'+this_qtl[3]+'\t'+qtl_genotype+'\t'+qtl_assessed_allele+'\t'+major+'\t'+minor+'\n')
