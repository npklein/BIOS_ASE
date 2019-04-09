import gzip
import glob

het_per_sample = {}
input_dir = '/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50_noRnaEditSites/'
for vcf in glob.glob(input_dir+'*genesAnnotated.vcf.gz'):
    print(vcf)
    chrom = vcf.split('chr')[1].split('.')[0]
    het_per_sample[chrom] = {}
    with gzip.open(vcf) as input_file:
        sample_index = {}
        for line in input_file:
            line = line.decode('ascii')
            if line.startswith('##'):
                continue
            data = line.strip().split('\t')
            info = data[7]
            data = data[9:]
            if line.startswith('#CHR'):
                for index, sample in enumerate(data):
                    sample_index[index] = sample
                    if sample not in het_per_sample[chrom]:
                        het_per_sample[chrom][sample] = {}
            else:
                if 'GENE=' in info:
                    gene = line.split('\t')[7].split('GENE=')[1]
                else:
                    gene = 'noGene'
                for index, info in enumerate(data):
                    if gene not in het_per_sample[chrom][sample_index[index]]:
                        het_per_sample[chrom][sample_index[index]][gene] = 0
                    if info.split(':')[0] == '0/1':
                        het_per_sample[chrom][sample_index[index]][gene] += 1

with open('hets_per_sample.txt','w') as out:
    out.write('sample\tchr\tgene\tnumberOfHets\n')
    for chr in het_per_sample:
        print(chr)
        for sample in het_per_sample[chr]:
            for gene in het_per_sample[chr][sample]:
                if het_per_sample[chr][sample][gene] != 0:
                    out.write(sample+'\t'+str(chr)+'\t'+gene+'\t'+str(het_per_sample[chr][sample][gene])+'\n')
