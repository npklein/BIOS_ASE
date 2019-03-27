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


#Read DDG2P database

open(DDG2P, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/geneAndVariantLists/DDG2P_DevelopmentDisorderGenotypePhenotypeDatabase_21_2_2018.txt") || die "Can't open file: DDG2P_DevelopmentDisorderGenotypePhenotypeDatabase_21_2_2018.txt!\n";
my @DDG2PFile = <DDG2P>;
close(DDG2P);
my $DDG2PHeader = $DDG2PFile[0];
chomp($DDG2PHeader);
my @DDG2PHeaderArray = split("\t", $DDG2PHeader);

my %DDG2P;
my %DDG2Panno;
for (my $i=1; $i<=$#DDG2PFile; $i++){ #loop over lines
    my $line = $DDG2PFile[$i];
    chomp($line);
    my @array = split("\t", $line); #split line by tab
    my $geneName = $array[0]; #extract gene name
    my $dddCategory = $array[4]; #DDD category
    my $allelicReq = $array[5]; #allelic requirement
    my $mutConsequence = $array[6]; #mutation consequence
    $geneName =~ s/^\s+|\s+$//g;
    $dddCategory =~ s/^\s+|\s+$//g;
    $allelicReq =~ s/^\s+|\s+$//g;
    $mutConsequence =~ s/^\s+|\s+$//g;
    $DDG2P{ $geneName } = $line; #push gene name as key in hash
    my $annotation = "$dddCategory\t$allelicReq\t$mutConsequence";
    $DDG2Panno{ $geneName } = $annotation;
}


#Loop over lines in alleleCounts file
open(ACF, "< alleleCountPerGroupPerGene.binom.annotated.alleleFiltered.removedCODAMandOutliers.splitOutliers.ALL.txt") || die "Can't open file: alleleCountPerGroupPerGene.medianSD3.merged.annotated.depthFiltered.removedCODAM.4outliersRemoved.V3.txt!\n";
open(OUTPUT, "> alleleCountPerGroupPerGene.binom.annotated.alleleFiltered.removedCODAMandOutliers.splitOutliers.DDG2P.txt") || die "Can't open file: alleleCountPerGroupPerGene.medianSD3.merged.annotated.depthFiltered.removedCODAM.4outliersRemoved.DDG2PgenesOnly.V3.txt!\n";

my $acfHeader = `head -1 alleleCountPerGroupPerGene.binom.annotated.alleleFiltered.removedCODAMandOutliers.splitOutliers.ALL.txt`;
chomp($acfHeader);
print OUTPUT "$acfHeader\tDDDcategory\tallelicRequirement\tmutationConsequence\n";

while (my $line = <ACF>) { #loop over file
    next if $. == 1; #skip header line
    chomp($line);
    my @array = split("\t", $line); #split line
    my $geneName = $array[24]; #extract gene  name
    if (exists $DDG2P{ $geneName }) { #if gene name exists in CGD hash, output line into new file.
        my $anno = $DDG2Panno{ $geneName };
        print OUTPUT "$line\t$anno\n";
    }
}
close(OUTPUT);
close(ACF);




