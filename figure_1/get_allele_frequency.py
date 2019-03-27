import gzip
import glob

with open('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig1/data/all_SNPs_AF_shapeitGenotypes.txt','w') as out:
    out.write('CHROM\tPOS\tAF\n')
    for f in glob.iglob('/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/shapeitVCFwithAC/chr*/*vcf.gz'):
        print(f)
        with gzip.open(f) as input_file:
            for line in input_file:
                line = line.decode('utf-8')
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                out.write(line[0]+'\t'+line[1]+'\t'+line[7].split('AF=')[1].split(';')[0]+'\n')
