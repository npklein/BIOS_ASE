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
use lib '/home/umcg-fvandijk/perl_modules/';
use lib '/home/umcg-fvandijk/perl_modules/lib/perl5/';
use Statistics::Basic qw(:all);



<<<<<<< HEAD
open(BINOM, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferonni.txt") || die "Can't open input BINOM file!\n";
=======
open(BINOM, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.txt") || die "Can't open input BINOM file!\n";
>>>>>>> 2a80f51c302edfb1917d2925ee56cc7cdaab0734
my @BINOMfile = <BINOM>;
close(BINOM);
my $BINOMhead = $BINOMfile[0];
chomp($BINOMhead);
my @BINOMheader = split("\t", $BINOMhead);

print "Processing BINOM file ..\n";
my %outlierHash;
my %geneHash;
for (my $i=1; $i<=$#BINOMfile; $i++){ #Process binom file line by line
    my $line = $BINOMfile[$i];
    chomp($line);
    my @array = split("\t", $line);
    #Gene is always the first element
    for (my $j=1; $j<=$#array; $j++){
        my $val = $array[$j];
        my $sample = $BINOMheader[$j];
        if ($val ne "NA") { #Not NA, so we measured something
            $geneHash{ $sample }++; #Increase count
            #If YES it was an outlier, count that too
            if ($val eq "YES") {
                $outlierHash{ $sample }++;
            }
        }
    }
}
print "Done processing BINOM file\n";

#Create output file

open(INPUT, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/phenotypeTables/AllSamples.phenotypes.txt") || die "Can't open input file!\n";
my @INPUTfile = <INPUT>;
close(INPUT);
my $INPUThead = $INPUTfile[0];
chomp($INPUThead);
my @INPUTheader = split("\t", $INPUThead);

open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/phenotypeTables/AllSamples.phenotypes.nGenesANDnOutliers.txt") || die "Can't open OUTPUT file!\n";
print OUTPUT "$INPUThead\tNGENES\tNOUTLIERS\n";
for (my $k=1; $k<=$#INPUTfile; $k++){
    my $line = $INPUTfile[$k];
    chomp($line);
    my @array = split("\t", $line);
    my $sample = $array[0];
    my $nOutliers = $outlierHash{ $sample };
    my $nGenes = $geneHash{ $sample };
    #print "$sample\t$nGenes\t$nOutliers\n";
    print OUTPUT "$line\t$nGenes\t$nOutliers\n";
}
close(OUTPUT);