import glob

pos_rs = {}
datadir = "/groups/umcg-bios/tmp03/projects/2020-03-05-ASE-BIOS-eqtl-comparison/data/"
for f in glob.glob(datadir+'genotypes-hrc-imputed-trityper/*/SNPMappings.txt'):
    with open(f) as input_file:
        for line in input_file:
            line = line.strip().split('\t')
            pos_rs[line[0]+'_'+line[1]] = line[2]

with open('counts.matrix.cumulativeVariants.ALLcounts.chrALL.txt.filtered.variants.txt') as input_file, open(datadir+'counts.matrix.cumulativeVariants.ALLcounts.chrALL.txt.filtered.variants.rsnumbers.txt','w') as out:
    out.write(input_file.readline())
    for line in input_file:
        line = line.split('_')
        if line[0]+'_'+line[1] in pos_rs:
            line[0] = pos_rs[line[0]+'_'+line[1]]
            out.write(line[0]+'\t'+line[-1].split('\t')[1])
