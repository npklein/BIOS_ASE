#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use List::Util qw(min max);
use List::Util qw(sum);
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);


my $tabixPath="/apps/software/HTSlib/1.3.2-foss-2015b/bin/";


#Create an array of all samples
print "Collecting all sample names..\n";
my @files = glob( "/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/geneAE_metaGenes/merged/chr1/*.chr1.geneAE.txt" );
my @samples;
foreach my $file (@files) { #Iterate over all geneAE files
    my $fileName = basename($file);
    #Extract part of filename to use for matching/coupling samples further downstream
    my $fID;
    if ($fileName =~ m/(.+).chr.+.geneAE.txt/gs) { #Extract samplename from allelic_counts filename
        $fID = $1;
        $fID =~ s/-lib1//gs;
    }
    push(@samples, $fID);
}
print "Done collecting all sample names.\n\n";


#/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/geneAE_metaGenes/merged/chr1/AD1NFNACXX-4-9.chr1.geneAE.txt
#contig	start	stop	name	aCount	bCount	totalCount	log2_aFC	n_variants	variants	gw_phased	bam
#1	169823221	169828897	ENSG00000000457	8	7	15	0.19264507794239583	1	1_169823660_C_T	1	BC43MJACXX-1-9.mdup.sorted.readGroupsAdded

my @geneAE_files = glob("/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/geneAE_metaGenes/merged/chr*/*.chr*.geneAE.txt");
print "\n# Files to process: " . ($#geneAE_files+1) . "\n\n";

my %geneSampleLogfold;
my %geneSampleCov;
my %geneSampleAcount;
my %geneSampleBcount;
my %genes;
my $count = 0;
foreach my $file (@geneAE_files) { #Iterate over all geneAE files
    my $fileName = basename($file);
    #Extract part of filename to use for matching/coupling samples further downstream
    my $fID;
    if ($fileName =~ m/(.+).chr.+.geneAE.txt/gs) { #Extract samplename from allelic_counts filename
        $fID = $1;
        $fID =~ s/-lib1//gs;
    }
    #Read file into array
    open(ALL, "< $file") || die "Can't open file: $file\n";
    my @all=<ALL>;
    close(ALL);
    #If file contains more than 1 line the alternative sampleID can be retrieved from the lines
    if ($#all > 0) {
        #More than one line, so process the file
        for (my $i=1; $i<=$#all; $i++){ #Loop over lines in allelic_counts file
            my $line = $all[$i];
            my @array = split("\t", $line);
            my $gene = $array[3];
            my $aCount = $array[4];
            my $bCount = $array[5];
            my $totCount = $array[6];
            my $logFold = $array[7];
            $geneSampleLogfold{ $gene }{ $fID } = $logFold;
            $geneSampleCov{ $gene }{ $fID } = $totCount;
            $geneSampleAcount{ $gene }{ $fID } = $aCount;
            $geneSampleBcount{ $gene }{ $fID } = $bCount;
            $genes{ $gene }++;
        }
    }
    #Print counter every 1000 lines processed
    if ($count % 1000 == 0){
        print "Processed $count files ..\n"
    }
    $count++;
    undef(@all);
}

print "\n\n\nGenerating output tables ..\n\n";
open(OUTA, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/binomialTest/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.20190124.aCount.txt") || die "Can't open outA!\n";
open(OUTB, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/binomialTest/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.20190124.bCount.txt") || die "Can't open outB!\n";

my $sampleString = join("\t", @samples);
my $outputHeader = "ENSEMBLID\t$sampleString";
print OUTA "$outputHeader\n";
print OUTB "$outputHeader\n";

#Iterate over all genes and output values per sample
foreach my $gene (sort keys %genes){
    if ($gene !~ m/.+;.+/) {
        print OUTA "$gene";
        print OUTB "$gene";
        #Iterate over all samples
        foreach my $sample (@samples){
            my $aCount = "NA";
            my $bCount = "NA";
            if (exists $geneSampleAcount{ $gene }{ $sample }) { #If logfold value is present continue
                $aCount = $geneSampleAcount{ $gene }{ $sample };
            }
            if (exists $geneSampleBcount{ $gene }{ $sample }) { #If logfold value is present continue
                $bCount = $geneSampleBcount{ $gene }{ $sample };
            }
            print OUTA "\t$aCount";
            print OUTB "\t$bCount";
        }
        print OUTA "\n";
        print OUTB "\n";
    }
}
close(OUTA);
close(OUTB);

print "Done generating output tables\n\n";
