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
        my $geneLength = ($stop-$start);
        if ($feature eq "gene") {
            if (looks_like_number($chr)) { #Check if chromosome is a number
                if ($info =~ m/gene_id "(ENSG[0-9]{1,})"; gene_name "(.+)"; gene_sourc.+/gs) { #Extract Ensembl gene ID
                    $geneID = $1;
                    $gene = $2;
                    $gene =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//;
                    $geneID =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//;
                    #if ($gene eq "PPT1" || $geneID eq "ENSG00000164458") {
                    #    print "$geneID\t$gene\n";
                    #}
                    $geneName2geneID{ $gene } = $geneID;
                    $geneID2geneName{ $geneID } = $gene;
                    $geneChr{ $geneID } = $chr;
                    $geneLength{ $geneID } = $geneLength;
                    $geneStart{ $geneID } = $start;
                    $geneStop{ $geneID } = $stop;
                    $geneCount++;
                }
            }
        }
    }
}
close(GTF);
print "#Genes detected: $geneCount\n";
print "Done processing GTF file.\n\n";


#Read Genes Of Interest file
print "Processing GOI file..\n";
my %GOI;
open(GOI, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/GOI.list.CADD15.AD.txt") || die "Can't open GOI file\n";
#open(GOI, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/GOI.list.ClinVarPATHOGENIC.txt") || die "Can't open GOI file\n";
while (my $lin = <GOI>){
    chomp($lin);
    my $gene = $lin;
    if (exists $geneName2geneID{ $gene }) {
        $geneID = $geneName2geneID{ $gene };
        $GOI{ $geneID } = $geneID;
    }
}
close(GOI);
print "Done processing GOI file.\n\n";


#Read gene expression file
print "Processing gene expression file..\n";
open(EXPR, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/geneExpressionTables/AllSamplesExcludingCODAMand4outliers.geneCounts.AllGenes.TPM.txt") || die "Can't open gene expression file\n";
my @EXPRfile = <EXPR>;
close(EXPR);
my $header = $EXPRfile[0];
chomp($header);
my @headerArray = split("\t", $header);
my %expression;
for (my $m=1; $m<=$#EXPRfile; $m++){
    my $line = $EXPRfile[$m];
    chomp($line);
    my @array = split("\t", $line);
    my $geneID = $array[0];
    #Iterate over values in line
    for (my $j=1; $j<=$#array; $j++){
        my $value = $array[$j];
        my $sample = $headerArray[$j];
        $expression{ $geneID }{ $sample } = $value;
    }
}
print "Done processing gene expression file.\n\n";


#Proces matrix files
print "Processing count matrix file..\n";
open(COUNT, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/MajorMinorCountsPerVariantPerSample.txt") || die "Can't open count matrix file!\n";
my @countFile = <COUNT>;
close(COUNT);
my %major;
my %minor;
for (my $i=1; $i<=$#countFile; $i++){
    my $line = $countFile[$i];
    chomp($line);
    my @array = split("\t", $line);
    my $var = $array[0];
    my $sample = $array[1];
    $sample =~ s|\.|-|g; #Replace dots
    my $majorC = $array[2];
    my $minorC = $array[3];
    $major{ $var }{ $sample } = $majorC;
    $minor{ $var }{ $sample } = $minorC;
}
print "Done processing count matrix file\n";



#Read Counts file
print "Processing counts file..\n";
my %clinvar;
open(COUNTS, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.chr22.addedCADD.addedVKGL.txt") || die "Can't open counts file\n";
my @COUNTSfile = <COUNTS>;
close(COUNTS);

open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/geneExpressionAndMajorMinorAlleleCounts.GOI.list.CADD15.AD.txt") || die "Can't open output file\n";
print OUTPUT "VARIANT\tCHR\tPOS\tREF\tALT\tGENEID\tGENENAME\tSAMPLE\tGENEEXPRESSION\tIMPACT\tMAJOR\tMINOR\tMINORRATIO\tLABEL\n";
for (my $k=1; $k<=$#COUNTSfile; $k++){
    my $line = $COUNTSfile[$k];
    chomp($line);
    my @array = split("\t", $line);
    my $var = $array[0];
    my $chr = $array[1];
    my $pos = $array[2];
    my $ref = $array[4];
    my $alt = $array[5];
    my $impact = $array[7];
    my $gene = $array[8];
    my $geneID = "NA";
    if (exists $geneName2geneID{ $gene }) {
        $geneID = $geneName2geneID{ $gene };
    }
    my $key = "$chr\_$pos\_$ref\_$alt";
    #If gene of interest continue
    if (exists $GOI{ $geneID }) {
        for (my $p=1; $p<=$#headerArray; $p++){ #Iterate over all samples
            my $sample = $headerArray[$p];
            my $expression = $expression{ $geneID }{ $sample };
            my $ratio = "NA";
            my $major = "NA";
            my $minor = "NA";
            my $label = "";
            if (exists $major{ $var }{ $sample }) {
                $major = $major{ $var }{ $sample };
                $minor = $minor{ $var }{ $sample };
                $ratio = ($minor/($major+$minor));
                $label = "$sample";
            }
            print OUTPUT"$var\t$chr\t$pos\t$ref\t$alt\t$geneID\t$gene\t$sample\t$expression\t$impact\t$major\t$minor\t$ratio\t$label\n"
        }
    }
}
close(OUTPUT);
print "Done processing counts file.\n\n";






