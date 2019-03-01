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

my $tabixPath="/apps/software/HTSlib/1.3.2-foss-2015b/bin/";

print "Reading GoNL AF file ..\n";
open(GONL, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/geneAndVariantLists/GoNL.AF.txt") || die "Can't open file: GoNLAF!\n";
my @GONLFile = <GONL>;
close(GONL);
my $GONLHeader = $GONLFile[0];
chomp($GONLHeader);
my @GONLHeaderArray = split("\t", $GONLHeader);

my %GoNLhash;
for (my $i=1; $i<=$#GONLFile; $i++){
    my $line = $GONLFile[$i];
    chomp($line);
    my @array = split("\t", $line);
    my $var = $array[0];
    my $af = $array[3];
    $GoNLhash{ $var } = $af;
}
print "Done reading GoNL AF file\n\n";


print "Reading ExAC AF file ..\n";
open(EXAC, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/geneAndVariantLists/ExAC.r0.3.1.sites.AF.EUR.txt") || die "Can't open file: ExACAF!\n";
my @EXACFile = <EXAC>;
close(EXAC);
my $EXACHeader = $EXACFile[0];
chomp($EXACHeader);
my @EXACHeaderArray = split("\t", $EXACHeader);

my %EXAChash;
for (my $i=1; $i<=$#EXACFile; $i++){
    my $line = $EXACFile[$i];
    chomp($line);
    my @array = split("\t", $line);
    my $var = $array[0];
    my $af = $array[3];
    $EXAChash{ $var } = $af;
}
print "Done reading ExAC AF file\n\n";


#/groups/umcg-bios/tmp03/projects/outlierGeneASE/pathogenicAlleles/alleleCountPerGroupPerGene.medianSD3.merged.annotated.txt

print "Processing annotated files ..\n";
for (my $CHR=1; $CHR<=22; $CHR++){
    my $file = "/groups/umcg-bios/tmp03/projects/outlierGeneASE/annotatedWith.snpEff.closest.VEP/BIOS_LLDeep_noRNAeditSites_phASER.snpEff.closest.VEP.chr$CHR.annotation.table";
    my $outFile  = "/groups/umcg-bios/tmp03/projects/outlierGeneASE/annotatedWith.snpEff.closest.VEP/BIOS_LLDeep_noRNAeditSites_phASER.snpEff.closest.VEP.chr$CHR.addedExACandGONLAlleleFrequency.annotation.table";
    print "Processing file for chr$CHR ..\n";
    open(ANN, "< $file") || die "Can't open file: $file!\n";
    my @ANNFile = <ANN>;
    close(ANN);
    my $ANNHeader = $ANNFile[0];
    chomp($ANNHeader);
    my @ANNHeaderArray = split("\t", $ANNHeader);

    
    open(OUTPUT, "> $outFile") || die "Can't open output file: $outFile!\n";
    print OUTPUT "$ANNHeader\tGoNL_AF_extracted\tExAC_AF_extracted\n";
#
    for (my $i=1; $i<=$#ANNFile; $i++){
        my $line = $ANNFile[$i];
        chomp($line);
        my @array = split("\t", $line);
        my $chr = $array[0];
        my $pos = $array[1];
        my $var = "$chr\_$pos";
        my $afGoNL = "NA";
        my $afExAC = "NA";
        if (exists $GoNLhash{ $var }) {
            $afGoNL = $GoNLhash{ $var };
            if (exists $EXAChash{ $var }) {
                $afExAC = $EXAChash{ $var };
            }
        }
        print OUTPUT "$line\t$afGoNL\t$afExAC\n";
    }
    close(OUTPUT);
}


print "Done processing annotated files\n\n";