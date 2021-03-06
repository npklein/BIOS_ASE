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


#Read OMIM database

open(OMIM, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/geneAndVariantLists/OMIM_20180211.mim2gene.txt") || die "Can't open file: OMIM_20180211.mim2gene.txt!\n";
my @omimFile = <OMIM>;
close(OMIM);


my %OMIM;
for (my $i=5; $i<=$#omimFile; $i++){ #loop over lines
    my $line = $omimFile[$i];
    chomp($line);
    my @array = split("\t", $line); #split line by tab
    my $entryType = $array[1]; #entry Type, should match gene
    if ($entryType eq "gene") {
        my $geneName = $array[3];
        if (defined $geneName && $geneName ne '') {
            $OMIM{ $geneName } = $line; #push gene name as key in hash
        }
    }
}


#Loop over lines in alleleCounts file
open(ACF, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/pathogenicAlleles/alleleCountPerGroupPerGene.binom.annotated.alleleFiltered.removedCODAMandOutliers.splitOutliers.ALL.txt") || die "Can't open file: alleleCountPerGroupPerGene.medianSD3.merged.annotated.depthFiltered.removedCODAM.4outliersRemoved.V3.txt!\n";
open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/pathogenicAlleles/alleleCountPerGroupPerGene.binom.annotated.alleleFiltered.removedCODAMandOutliers.splitOutliers.OMIM.txt") || die "Can't open file: alleleCountPerGroupPerGene.medianSD3.merged.annotated.depthFiltered.removedCODAM.4outliersRemoved.OMIMgenesOnly.V3.txt!\n";

my $acfHeader = `head -1 /groups/umcg-bios/tmp03/projects/outlierGeneASE/pathogenicAlleles/alleleCountPerGroupPerGene.binom.annotated.alleleFiltered.removedCODAMandOutliers.splitOutliers.ALL.txt`;
chomp($acfHeader);
print OUTPUT "$acfHeader\n";

while (my $line = <ACF>) { #loop over file
    next if $. == 1; #skip header line
    chomp($line);
    my @array = split("\t", $line); #split line
    my $geneName = $array[24]; #extract gene  name
    if (exists $OMIM{ $geneName }) { #if gene name exists in CGD hash, output line into new file.
        print OUTPUT "$line\n";
    }
}
close(OUTPUT);
close(ACF);




