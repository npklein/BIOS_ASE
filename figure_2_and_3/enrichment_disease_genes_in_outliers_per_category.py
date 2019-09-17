

inheritance = {}
manifestation = {}
manifestations = set([])
inheritances = set([])
# TOD: how was this file made
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/geneAndVariantLists/CGD.20171220.txt') as input_file:
    header = input_file.readline().strip().split('\t')
    for line in input_file:
        line = line.strip().split('\t')
        if line[0] not in inheritance:
            inheritance[line[0]] = set([])
        inheritance[line[0]].add(line[4].strip())
        inheritances.add(line[4].strip())
        for manifestationName in line[7].split(';'):
            manifestationName = manifestationName.strip()
            if line[0] not in manifestation:
                manifestation[line[0]] = set([])
            manifestation[line[0]].add(manifestationName)
            manifestations.add(manifestationName)

logFC_per_sample_per_gene = {}
# TODO: how was this file made 
input_dir = '/groups/umcg-bios/tmp03/projects/outlierGeneASE//logFoldChangeTables/'
with open(input_dir+'genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.samplesFILTERED.values.txt') as input_file:
    header = input_file.readline().strip().split('\t')
    header_index = {}
    for index, element in enumerate(header):
        header_index[index] = element
        logFC_per_sample_per_gene[element] = {}
    for line in input_file:
        line = line.strip().split('\t')
        gene = line[0]
        for index, element in enumerate(line):
            sample = header_index[index]
            logFC_per_sample_per_gene[sample][gene] = element

outliers = {}
not_outliers = {}
na = {}
# TODO: how was this file made 
with open(input_dir+'genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.samplesFILTERED.txt') as input_file:
    header = input_file.readline().strip().split('\t')
    for line in input_file:
        line = line.strip().split('\t')
        gene = line[0]
        outliers[gene] = set([])
        not_outliers[gene] = set([])
        na[gene] = set([])
        for index, element in enumerate(line):
            if element == 'YES':
                outliers[gene].add(header[index])

            elif element == 'NO':
                not_outliers[gene].add(header[index])

            elif element == 'NA':
                na[gene].add(header[index])


count_per_manifestation = {}
count_per_inheritance = {}
for manifestationName in manifestations:
    manifestationName = manifestationName.strip()
    count_per_manifestation[manifestationName] = {'outlier':{}, 'not_outlier':{}, 'na':{}}
count_per_manifestation['other'] = {'outlier':{}, 'not_outlier':{}, 'na':{}}
for inheritanceName in inheritances:
    inheritanceName = inheritanceName.strip()
    count_per_inheritance[inheritanceName] = {'outlier':{} , 'not_outlier':{},'na':{}}
count_per_inheritance['other'] = {'outlier':{}, 'not_outlier':{}, 'na':{}}


s = set([])
s2 = set([])

with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/omim_enrichment/omim_carriers/OMIM_carriers_allIMPACT.txt') as input_file:
    input_file.readline()
    x = 0
    for line in input_file:
        x += 1
        if x % 100 == 0:
            print(x,'lines processed')
        line = line.strip().split('\t')
        gene, sample, snp, genotype, symbol, impact, isOmim, MAF = line[0],line[2],line[1],line[3],line[4],line[7],line[8], line[6]
        s2.add(gene)
        if gene not in outliers and gene not in not_outliers:
             s.add(gene)
             continue
        if line[2] in outliers[gene]:
            if symbol in manifestation:
                for manifestationName in manifestation[symbol]:
                    if impact not in count_per_manifestation[manifestationName]['outlier']:
                        count_per_manifestation[manifestationName]['outlier'][impact] = []
                    count_per_manifestation[manifestationName]['outlier'][impact].append([gene, snp, sample, MAF])
                for inheritanceName in inheritance[symbol]:
                    if impact not in count_per_inheritance[inheritanceName]['outlier']:
                        count_per_inheritance[inheritanceName]['outlier'][impact] = []
                    count_per_inheritance[inheritanceName]['outlier'][impact].append([gene, snp, sample, MAF])
            else:
                if impact not in count_per_manifestation['other']['outlier']:
                    count_per_manifestation['other']['outlier'][impact] = []
                    count_per_inheritance['other']['outlier'][impact] = []
                count_per_manifestation['other']['outlier'][impact].append([gene, snp, sample, MAF])
                count_per_inheritance['other']['outlier'][impact].append([gene, snp, sample, MAF])
        elif line[2] in not_outliers[gene]:
            if symbol in manifestation:
                for manifestationName in manifestation[symbol]:
                    if impact not in count_per_manifestation[manifestationName]['not_outlier']:
                        count_per_manifestation[manifestationName]['not_outlier'][impact] = []
                    count_per_manifestation[manifestationName]['not_outlier'][impact].append([gene, snp, sample, MAF])
                for inheritanceName in inheritance[symbol]:
                    if impact not in count_per_inheritance[inheritanceName]['not_outlier']:
                        count_per_inheritance[inheritanceName]['not_outlier'][impact] = []
                    count_per_inheritance[inheritanceName]['not_outlier'][impact].append([gene, snp, sample, MAF])
            else:
                if impact not in count_per_manifestation['other']['not_outlier']:
                    count_per_manifestation['other']['not_outlier'][impact] = []
                    count_per_inheritance['other']['not_outlier'][impact] = []
                count_per_manifestation['other']['not_outlier'][impact].append([gene, snp, sample, MAF])
                count_per_inheritance['other']['not_outlier'][impact].append([gene, snp, sample, MAF])
        elif line[2] in na[gene]:
            if symbol in manifestation:
                for manifestationName in manifestation[symbol]:
                    if impact not in count_per_manifestation[manifestationName]['na']:
                        count_per_manifestation[manifestationName]['na'][impact] = []
                    count_per_manifestation[manifestationName]['na'][impact].append([gene, snp, sample, MAF])
                for inheritanceName in inheritance[symbol]:
                    if impact not in count_per_inheritance[inheritanceName]['na']:
                         count_per_inheritance[inheritanceName]['na'][impact] = []
                    count_per_inheritance[inheritanceName]['na'][impact].append([gene, snp, sample, MAF])
            else:
                if impact not in count_per_manifestation['other']['na']:
                    count_per_manifestation['other']['na'][impact] = []
                    count_per_inheritance['other']['na'][impact] = []
                count_per_manifestation['other']['na'][impact].append([gene, snp, sample, MAF])
                count_per_inheritance['other']['na'][impact].append([gene, snp, sample, MAF])


disease_outfile = '/groups/umcg-bios/tmp03/projects/outlierGeneASE/omim_enrichment/carriers_per_disease/carriers_per_disease.txt'
with open(disease_outfile,'w') as out:
    out.write('disease\ttype\timpact\tgene\tsnp\tsample\tlogFC\tMAF\n')
    for m in count_per_manifestation:
        for type in ['outlier','not_outlier','na']:
            for impact in count_per_manifestation[m][type]:
                for result in count_per_manifestation[m][type][impact]:
                    if result[2] in logFC_per_sample_per_gene:
                        if result[0] in logFC_per_sample_per_gene[result[2]]:
                            out.write(m+'\t'+type+'\t'+impact+'\t'+result[0]+'\t')
                            out.write(result[1]+'\t'+result[2]+'\t')
                            out.write(logFC_per_sample_per_gene[result[2]][result[0]])
                            out.write('\t'+result[3]+'\n')
                        else:
                            print(result[0],'not in logFC_per_sample_per_gene[result[2]')
                    else:
                        print(result[2],'not in logFC_per_sample_per_gene')

inheritance_outfile = '/groups/umcg-bios/tmp03/projects/outlierGeneASE/omim_enrichment/carriers_per_inheritance/carriers_per_inheritance.txt'
with open(inheritance_outfile,'w') as out:
    out.write('inheritance\ttype\timpact\tgene\tsnp\tsample\tlogFC\tMAF\n')
    for m in count_per_inheritance:
        for type in ['outlier','not_outlier','na']:
            for impact in count_per_inheritance[m][type]:
                for result in count_per_inheritance[m][type][impact]:
                    if result[2] in logFC_per_sample_per_gene:
                        if result[0] in logFC_per_sample_per_gene[result[2]]:
                            out.write(m+'\t'+type+'\t'+impact+'\t'+result[0]+'\t')
                            out.write(result[1]+'\t'+result[2]+'\t')
                            out.write(logFC_per_sample_per_gene[result[2]][result[0]])
                            out.write('\t'+result[3]+'\n')

print('output written to',disease_outfile,'and',inheritance_outfile)
