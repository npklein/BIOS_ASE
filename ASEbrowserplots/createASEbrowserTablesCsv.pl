#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use List::Util qw(min max);
use List::Util qw(sum);
use List::Util 'shuffle';
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);


###Variables to set
my $outputSep = "\t"; #Make it \t for testing, comma for production to create csv file.
my $tabixPath="/apps/software/HTSlib/1.3.2-foss-2015b/bin/";


##Read GTF file
print "Processing GTF file..\n";
open(GTF, "< /apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf") || die "Can't open GTF file\n";
my %geneName2geneID;
my %geneID2geneName;
my $geneCount=0;
my %geneChr;
my %geneLength;
my %geneStart;
my %geneStop;
my %geneBiotype;

open(GENES, "> /groups/umcg-bios/tmp03/projects/BIOS_manuscript/ase_genes.txt") || die "Can't open file: ase_genes.txt\n";
print GENES "\"Gene_symbol\"" . $outputSep . "\"Ensembl_ID\"" . $outputSep . "\"Biotype\"" . "\n";
while (my $line=<GTF>){ #Read GTF file line by line
    chomp($line);
    if ($line !~ m/^#.+/gs) { #If not header line process further
        my @array = split("\t", $line);
        my $chr = $array[0];
        #$chr =~ s/X/23/gs;
        #$chr =~ s/Y/24/gs;
        #$chr =~ s/MT/25/gs;
        my $feature = $array[2];
        my $start = $array[3];
        my $stop = $array[4];
        my $info = $array[8];
        my $geneID = "NA";
        my $gene = "NA";
        my $gBiotype = "NA";
        my $geneLength = ($stop-$start);
        if ($feature eq "gene") {
            #if (looks_like_number($chr)) { #Check if chromosome is a number
                if ($info =~ m/gene_id "(ENSG[0-9]{1,})"; gene_name "(.+)"; gene_sourc.+gene_biotype "(.+)".+/gs) { #Extract Ensembl gene ID
                    $geneID = $1;
                    $gene = $2;
                    $gBiotype = $3;
                    $gene =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//;
                    $geneID =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//;
                    $geneName2geneID{ $gene } = $geneID;
                    $geneID2geneName{ $geneID } = $gene;
                    $geneChr{ $geneID } = $chr;
                    $geneLength{ $geneID } = $geneLength;
                    $geneStart{ $geneID } = $start;
                    $geneStop{ $geneID } = $stop;
                    $geneBiotype{ $geneID } = $gBiotype;
                    $geneCount++;
                    print GENES "\"$gene\"" . $outputSep . "\"$geneID\"" . $outputSep . "\"$gBiotype\"" . "\n";
                }
            #}
        }
    }
}
close(GENES);
close(GTF);
print "#Genes detected: $geneCount\n";
print "Done processing GTF file.\n\n";


##Process annotation tables
print "Processing annotation tables ..\n";
my @annotationFiles = glob('/groups/umcg-bios/tmp03/projects/outlierGeneASE/annotatedWith.snpEff.closest.VEP/BIOS_LLDeep_noRNAeditSites_phASER.snpEff.closest.VEP.chr*.addedExACandGONLAlleleFrequency.addedpLI.annotation.table');
#Loop through these files to extract annotation information per variant
my %snpIDs;
my %varGeneName;
my %varGeneID;
foreach my $anFile (@annotationFiles){
    open(ANNO, "< $anFile") || die "Can't open annotation file: $anFile\n";
    my @annotationFile = <ANNO>;
    close(ANNO);
    for (my $i=1; $i<=$#annotationFile; $i++){
        my $line = $annotationFile[$i];
        chomp($line);
        my @array = split("\t", $line);
        my $chr = $array[0];
        my $pos = $array[1];
        my $rsID = $array[2];
        my $ref = $array[3];
        my $alt = $array[4];
        my $annotation = $array[8];
        my $impact = $array[9];
        my $geneName = $array[10];
        my $geneID = $array[11];
        my $vepConsequence = $array[32];
        my $vepImpact = $array[33];
        my $vepFeatureType = $array[36];
        my $vepBiotype = $array[38];
        $geneName =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//;
        my $gonlAF = $array[140];
        my $exacAF = $array[141];
        my $pLI = $array[142];
        my $key = "$chr\_$pos\_$ref\_$alt";
        if ($rsID =~ m/\./gs) {
            $rsID = "$chr:$pos";
        }
        $snpIDs{ $key } = $rsID;
        $varGeneName{ $key } = $geneName;
        $varGeneID{ $key } = $geneID;
    }
}
print "Done processing annotation tables.\n\n";


##Read sample file to include only pass QC samples
print "Processing sample file..\n";
open(SAMPLE, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/samples_NOUTLIERS1000.depthFiltered.bonferroni.txt") || die "Can't open sample file\n";
my @samples;
while (my $lin = <SAMPLE>) { #Only one ID per line in the file
    chomp($lin);
    push(@samples, $lin);
}
close(SAMPLE);
print "Done processing sample file.\n\n";


##Create ase_sampleAse table
print "Processing ASEcount files..\n";
my %varIDs;
my %varSampleRefAll;
my %varSampleRefC;
my %varSampleAltC;
my @TotVarCount;
for (my $CHR=1; $CHR<=22; $CHR++){
    print "Processing files for chr$CHR ..\n";
    my @count_files;
    foreach my $sample (@samples){
        push( @count_files, "/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/readbackedPhasing/allelic_counts_mergedPerSample/chr$CHR/BIOS_LLDeep_noRNAeditSites_phASER.$sample.chr$CHR.allelic_counts.txt");
    }
    foreach my $countFile (@count_files){
        open(COUNTS, "< $countFile") || die "Can't open file: $countFile\n";
        my @file = <COUNTS>;
        close(COUNTS);
        #Extract sampleID from filename
        my $sampleID = "NA";
        if ($countFile =~ m/.+BIOS_LLDeep_noRNAeditSites_phASER.(.+).chr$CHR.allelic_counts.txt/gs) {
            $sampleID = $1;
        }
        for (my $i=1; $i <= $#file; $i++){ #Loop over all lines in file
            my $line = $file[$i];
            chomp($line);
            my @array = split("\t", $line);
            my $chr = $array[0];
            my $pos = $array[1];
            my $varID = $array[2];
            my $refAll = $array[3];
            my $altAll = $array[4];
            my $refC = $array[5];
            my $altC = $array[6];
            my $totC = $array[7];
            if ($refC >= 10 && $altC >= 10 && $totC >= 20) {
                #Push all info into hashes
                $varIDs{ $varID }++; #Number of samples carrying this variant
                $varSampleRefAll{ $varID }{ $sampleID } = $refAll;
                $varSampleRefC{ $varID }{ $sampleID } = $refC;
                $varSampleAltC{ $varID }{ $sampleID } = $altC;
                #push(@TotVarCount, $varID); #Needed for random number selection
            }
        }
        #undef(@file);
    }
    undef(@count_files);
}
print "Creating output table..\n";


#Write output
open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/BIOS_manuscript/ase_sampleAse.txt") || die "Can't open output file: ase_sampleAse.txt\n";
open(SUMMEDOUTPUT, "> /groups/umcg-bios/tmp03/projects/BIOS_manuscript/ase_ase.txt") || die "Can't open output file: ase_ase.txt\n";
print OUTPUT "\"snp_id\"" . $outputSep . "\"Ref_Counts\"" . $outputSep . "\"Alt_Counts\"" . $outputSep . "\"Chromosome\"" . $outputSep . "\"Position\"" . $outputSep . "\"ID\"" . $outputSep . "\"Pval\"" . $outputSep . "\"Bonf_corr_Pval\"" . $outputSep . "\"Gene_symbol\"" . $outputSep . "\"Ensembl_ID\"" . "\n";
print SUMMEDOUTPUT "\"Fraction_alternative_allele\"" . $outputSep . "\"Alternative_allele\"" . $outputSep . "\"Reference_allele\"" . $outputSep . "\"Samples\"" . $outputSep . "\"SNP_ID\"" . $outputSep . "\"Chr\"" . $outputSep . "\"Pos\"" . $outputSep . "\"Likelihood_ratio_test_D\"" . $outputSep . "\"P_Value\"" . $outputSep . "\"Genes\"" . "\n";
my $increment = 1;
foreach my $variant (sort keys %varIDs){ #Foreach variant
    my $sumRefC = 0;
    my $sumAltC = 0;
    my $chr = "NA";
    my $pos = "NA";
    my $ref = "NA";
    my $alt = "NA";
    my $nSamples = 0;
    my $snpID = "NA";
    foreach my $sample (shuffle keys %{ $varSampleRefAll{ $variant }}){
        my $refAll = $varSampleRefAll{ $variant }{ $sample };
        my $refC = $varSampleRefC{ $variant }{ $sample };
        my $altC = $varSampleAltC{ $variant }{ $sample };
        my $geneSymbol = $varGeneName{ $variant };
        my $geneID = $varGeneID{ $variant };
        $snpID = $snpIDs{ $variant };
        $nSamples = $varIDs{ $variant };
        chomp($variant);
        my @array = split("_", $variant);
        $chr = $array[0];
        $pos = $array[1];
        $ref = $array[2];
        $alt = $array[3];
        if ($ref ne $refAll) {
            exit("ERROR: variant $variant is not matching reference allele for variant and sample $sample\n");
        }
        print OUTPUT "\"$snpID\"" . $outputSep . "\"$refC\"" . $outputSep . "\"$altC\"" . $outputSep . "\"$chr\"" . $outputSep . "\"$pos\"" . $outputSep;
        printf( OUTPUT "\"%08d\"",$increment);
        #print "\t$sample";
        print OUTPUT $outputSep . "\"NA\"" . $outputSep . "\"NA\"" . $outputSep . "$geneSymbol" . $outputSep . "$geneID" . "\n";
        $increment++;
        
        #Collect counts for one variant
        $sumRefC = $sumRefC + $refC;
        $sumAltC = $sumAltC + $altC;
    }
    #Do calculations for one variant
    my $total = ($sumRefC+$sumAltC);
    my $fraction = ($sumAltC/$total);
    print SUMMEDOUTPUT "\"$fraction\"" . $outputSep . "\"$alt\"" . $outputSep . "\"$ref\"" . $outputSep . "\"$nSamples\"" . $outputSep . "\"$snpID\"" . $outputSep . "\"$chr\"" . $outputSep. "\"$pos\"" . $outputSep . "\"NA\"" . $outputSep . "\"NA\"" . $outputSep . "\"NA\"" . "\n";
}
close(OUTPUT);
close(SUMMEDOUTPUT);
print "Done creating output table.\n";
print "Done processing ASEcount files.\n";


