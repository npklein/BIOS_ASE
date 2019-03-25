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
use POSIX qw/ceil/;
use Scalar::Util qw(looks_like_number);


###Variables to set
my $tabixPath="/apps/software/HTSlib/1.3.2-foss-2015b/bin/";





##Read ASE file
print "Processing ASE file..\n";
my %varSampleLine;
my %snps;
#open(ASE, "< /groups/umcg-bios/tmp03/projects/BIOS_manuscript/test.test.txt") || die "Can't open ASE file\n";
open(ASE, "< /groups/umcg-bios/tmp03/projects/BIOS_manuscript/ASEbrowserTables/ase_sampleAse.all.binom.txt") || die "Can't open ASE file\n";
my @ASEfile = <ASE>;
close(ASE);

for (my $i = 1; $i<=$#ASEfile; $i++){
    my $line = $ASEfile[$i];
    chomp($line);
    my @array = split("\t", $line);
    my $snpid = $array[0];
    my $id = $array[5];
    $varSampleLine{ $snpid }{ $id } = "$line";
    $snps{ $snpid }++;
}
print "Done processing ASE file.\n\n";



print "Generating output file..\n";
my $nChunks = 100;
my @sorted_vals = (keys %snps);
my $chunk = ceil(($#sorted_vals+1)/$nChunks);

#print "$chunk\n";

for (my $k = 1; $k<=$nChunks; $k++){
    my @sub = splice(@sorted_vals, 0, $chunk);
    open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/BIOS_manuscript/ASEbrowserTables/ase_sampleAse.chunk$k.txt") || die "Can't open output file: ase_sampleAse.chunk$k.txt\n";
    my $header = $ASEfile[0];
    print OUTPUT "$header";
    foreach my $var (@sub){
        #print "$k\t\t$var\n";
        foreach my $sample (keys %{ $varSampleLine{ $var }}){
            my $line = $varSampleLine{ $var }{ $sample };
            #print "$k\t$var\t$sample\t$line\n";
            print OUTPUT "$line\n";
        }
    }
    close(OUTPUT);
}
print "Done generating output file.\n\n";




