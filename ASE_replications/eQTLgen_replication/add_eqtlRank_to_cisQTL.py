
eqtl_rank = {}
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithEqtlGen/data/eqtlRank.txt') as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        if len(line) == 1:
            continue
        eqtl_rank[line[0]+'_'+line[1]] = line[2]

f_in = '/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithEqtlGen/data/2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved.aseSNPs.snpsInProtCodingExon.txt'
f_out = '/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithEqtlGen/data/2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved.aseSNPs.snpsInProtCodingExon.eqtlGenRank.txt'
with open(f_in) as input_file,  open(f_out,'w') as out:
    header = input_file.readline().strip()+'\teqtlGenRank\n'
    out.write(header)
    for line in input_file:
        line = line.strip()
        chr = line.split('\t')[3]
        pos = line.split('\t')[4]
        snp = chr+'_'+pos
        eqtlGen_rank = eqtl_rank[line.split('\t')[0]+'_'+snp]
        out.write(line+'\t'+eqtlGen_rank+'\n')
