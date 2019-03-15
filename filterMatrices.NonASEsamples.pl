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

my $tabixPath="/apps/software/HTSlib/1.3.2-foss-2015b/bin/";


##Read counts file
print "Processing count matrix files..\n";
my @countFiles = glob('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.m*orAllelle.nonASEsamples.chrALL.txt');


foreach my $cFile (@countFiles){
    print "Processing file $cFile..\n";
    open(COUNT, "< $cFile") || die "Can't open count file: $cFile\n";
    my @countFile = <COUNT>;
    close(COUNT);
    my $header = $countFile[0];
    chomp($header);
    
    #Open output
    open(OUTPUT, "> $cFile.filtered.txt") || die "Can't open output file: $cFile.filtered.txt\n";
    print OUTPUT "$header\n";
    for (my $i=1; $i<=$#countFile; $i++){
        my $line = $countFile[$i];
        chomp($line);
        my @array = split("\t", $line);
        #Loop over all elements, from 1 to last
        #If a non NA occurs write the line and break the loop
        for (my $j=1; $j<=$#array; $j++){
            my $element = $array[$j];
            if ($element ne "NA") {
                print OUTPUT "$line\n";
                last;
            }
        }
    }
}
print "Done processing count matrix files\n";

