use warnings;
use diagnostics;
use List::Util qw(min max);
use List::Util qw(sum);
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $tabixPath="/apps/software/HTSlib/1.3.2-foss-2015b/bin/";


#Read Counts file
print "Processing counts file..\n";
my %clinvar;
open(COUNTS, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.chr22.txt") || die "Can't open counts file\n";
my @COUNTSfile = <COUNTS>;
my $header = $COUNTSfile[0];
chomp($header);

open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.chr22.addedCADD.txt") || die "Can't open output file\n";
print OUTPUT "$header\tCADDRAW\tCADDPHRED\n";
for (my $m=1; $m<=$#COUNTSfile; $m++){
    my $line = $COUNTSfile[$m];
    chomp($line);
    my @array = split("\t", $line);
    my $chr = $array[1];
    my $pos = $array[2];
    my $ref = $array[4];
    my $alt = $array[5];
    #Use tabix to extract region from CADD file
    my $cmd = "$tabixPath/tabix /apps/data/CADD/whole_genome_SNVs.tsv.gz $chr:$pos-$pos";
    my $exec = `$cmd`;
    #print "$line\n$exec\n\n";
    my @results = split("\n", $exec);
    foreach my $lin (@results){
        chomp($lin);
        my @res = split("\t", $lin);
        my $rChr = $res[0];
        my $rPos = $res[1];
        my $rRef = $res[2];
        my $rAlt = $res[3];
        my $rRaw = $res[4];
        my $rPhred = $res[5];
        if ($ref eq $rRef && $alt eq $rAlt){ #Ref and Alt allele both match, so match
            print OUTPUT "$line\t$rRaw\t$rPhred\n";
            last; #Break out of loop
        }
    }
}
close(COUNTS);
print "Done processing counts file.\n\n";

