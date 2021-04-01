


import gzip
x = 1
with gzip.open('geneExpressionAndMajorMinorAlleleCounts.allGenes.txt.gz','rt') as input_file, gzip.open('geneExpressionAndMajorMinorAlleleCounts.allGenes.sorted.txt.gz','wt') as out:
    gene_info = {}
    out.write(input_file.readline())
    for line in input_file:
        x+= 1
        if x % 1000000 == 0:
            print(x,'lines processed')
        gene = line.split('\t')[5]
        if gene not in gene_info:
            gene_info[gene] = []
        gene_info[gene].append(line)
    for gene in gene_info:
        for line in gene_info[gene]:
            out.write(line)