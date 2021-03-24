dir='/groups/umcg-bios/tmp04/projects/copy_from_tmp03/outlierGeneASE/geneAndVariantLists/'
with open(dir+'OMIM_20180211.mim2gene.txt') as input_file, open(dir+'OMIM.20171220.ensembleGenes.txt','w') as out:
    for line in input_file:
        line = line.strip()
        gene = line.split('\t')[-1]
        try:
            symbol = line.split('\t')[-2]
            if len(symbol.strip()) == 0:
                symbol = 'NA'
        except:
            symbol = 'NA'
        if gene.startswith('ENSG'):
            out.write(gene+'\t'+symbol+'\n')
