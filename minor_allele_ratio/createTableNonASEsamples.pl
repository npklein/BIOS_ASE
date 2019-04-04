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


#Read CGD file
print "Processing CGD file..\n";
open(CGD, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/geneAndVariantLists/CGD.20171220.txt") || die "Can't open CGD file\n";
my @CGDfile = <CGD>;
my %inheritance;
my %ageGroup;
my %manifestation;
my %CGDgenes;
for (my $j=1; $j<=$#CGDfile; $j++){
    my $line = $CGDfile[$j];
    chomp($line);
    my @array = split("\t", $line);
    my $gene = $array[0];
    my $inheritance = $array[4];
    $inheritance =~ s/ //gs;
    my $ageGroup = $array[5];
    my $manifestation = $array[7];
    $inheritance{ $gene } = $inheritance;
    $ageGroup{ $gene } = $ageGroup;
    $manifestation{ $gene } = $manifestation;
    $CGDgenes{ $gene } = $gene;
}
close(CGD);
print "Done processing CGD file.\n\n";


#Read CLINVAR file
print "Processing CLINVAR file..\n";
my %clinvar;
open(CLNVR, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/clinvar/clinvar_data/clinvar_20180128.criteriaProvided_multiple_submitters_no_conflicts.txt") || die "Can't open CLINVAR file\n";
my @CLNVRfile = <CLNVR>;
for (my $m=1; $m<=$#CLNVRfile; $m++){
    my $line = $CLNVRfile[$m];
    chomp($line);
    my @array = split("\t", $line);
    my $chr = $array[1];
    my $pos = $array[2];
    my $ref = $array[3];
    my $alt = $array[4];
    my $clnsig = $array[15];
    my $variant = "$chr\_$pos\_$ref\_$alt";
    $clinvar{ $variant } = $clnsig;
}
close(CLNVR);
print "Done processing CLINVAR file.\n\n";

#Genes with high pLI scores (pLI <B3> 0.9) are extremely LoF intolerant, whereby genes with low pLI scores (pLI <B2> 0.1) are LoF tolerant.

print "Processing annotation tables ..\n";
my @annotationFiles = glob('/groups/umcg-bios/tmp03/projects/outlierGeneASE/annotatedWith.snpEff.closest.VEP/BIOS_LLDeep_noRNAeditSites_phASER.snpEff.closest.VEP.chr*.addedExACandGONLAlleleFrequency.addedpLI.annotation.table');
#Loop through these files to extract annotation information per variant
my %annotations;
my %impacts;
my %exacs;
my %gonls;
my %plis;
my %geneplis;
my %geneImpacts;
my %geneImpactHigh;
my %geneImpactModerate;
my %geneImpactModifier;
my %geneImpactLow;
my %genes;
my %varToGenes;
my $totalGenes=0;
my $tolGenes=0;
my $intolGenes=0;
my $nonGenes=0;
my $naGenes=0;
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
        my $gene = $array[10];
        my $vepConsequence = $array[32];
        my $vepImpact = $array[33];
        my $vepFeatureType = $array[36];
        my $vepBiotype = $array[38];
        $gene =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//;
        my $gonlAF = $array[140];
        my $exacAF = $array[141];
        my $pLI = $array[142];
        my $key = "$chr\_$pos\_$ref\_$alt";
        $annotations{ $key } = "$chr\t$pos\t$rsID\t$ref\t$alt\t$annotation\t$impact\t$gene\t$vepConsequence\t$vepImpact\t$vepFeatureType\t$vepBiotype\t$gonlAF\t$exacAF\t$pLI";
        $impacts{ $key } = $impact;
        $exacs{ $key } = $exacAF;
        $gonls{ $key } = $gonlAF;
        $plis{ $key } = $pLI;
        $geneplis{ $gene } = $pLI;
        if (not exists $genes{ $gene }) {
            $totalGenes++;
            if (looks_like_number($pLI)) {
                if ($pLI >= 0.9) {
                    $intolGenes++;
                }elsif ($pLI <= 0.1){
                    $tolGenes++;
                }else{
                    $nonGenes++;
                }
            }else{
                $naGenes++;
            }
        }
        $varToGenes{ $key } = $gene;
        $genes{ $gene } = $gene;

        $geneImpacts{ $gene }{ $impact }++; #Counts number of impacts per gene
        #print "$key\t$annotation\t$impact\n";
    }
}
print "#Total genes: $totalGenes\n";
print "#Tolerant genes: $tolGenes\n";
print "#Intolerant genes: $intolGenes\n";
print "#Neither tolerant or intolerant genes: $nonGenes\n";
print "#NA genes: $naGenes\n";

print "Done processing annotation tables.\n\n";


#Read AF genotype data file
print "Processing AF file..\n";
open(AF, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/annotatedWith.snpEff.closest.VEP/chrALL.AFsFromData.txt") || die "Can't open AF file\n";
my @AFfile = <AF>;
my %afs;
for (my $j=0; $j<=$#AFfile; $j++){
    my $line = $AFfile[$j];
    chomp($line);
    my @array = split("\t", $line);
    my $key = $array[0];
    my $af = $array[1];
    $afs{ $key } = $af;
}
close(AF);
print "Done processing AF file.\n\n";


open(BINOM, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.samplesFILTERED.txt") || die "Can't open input BINOM file!\n";
my @BINOMfile = <BINOM>;
close(BINOM);
my $BINOMhead = $BINOMfile[0];
chomp($BINOMhead);
my @BINOMheader = split("\t", $BINOMhead);

print "Processing BINOM file ..\n";
my %outlierHash;
my %nonoutlierHash;
my %binomStatusHash;
my %geneHash;
for (my $i=1; $i<=$#BINOMfile; $i++){ #Process binom file line by line
    my $line = $BINOMfile[$i];
    chomp($line);
    my @array = split("\t", $line);
    #Gene is always the first element
    my $gene = $array[0];
    for (my $j=1; $j<=$#array; $j++){
        my $val = $array[$j];
        my $sample = $BINOMheader[$j];
        if ($val ne "NA") { #Not NA, so we measured something
            $binomStatusHash{ $sample }{ $gene } = $val;
            $geneHash{ $sample }++; #Increase count
            #If YES it was an outlier, count that too
            if ($val eq "YES") {
                $outlierHash{ $sample }{ $gene } = $val;
            }elsif ($val eq "NO"){
                $nonoutlierHash{ $sample }{ $gene } = $val;
            }
            #print "$sample\t$gene\t$val\n";
        }
    }
}
print "Done processing BINOM file\n\n";





#Process haplotype files
print "Processing haplotypic count files ..\n";
#my $hapFile = "./chr$CHR.TEST.txt";
my %varRefCount;
my %varAltCount;
my %varStatusCount;
my %varSampleRefCount;
my %varSampleAltCount;
for (my $CHR=1; $CHR<=22; $CHR++){
    my $hapFile = "/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/geneAE_metaGenes/merged/chr$CHR.txt";
    print "Processing file for chr$CHR ..\n";
    open(HAP, "< $hapFile") || die "Can't open file: $hapFile!\n";
    my @HAPfile = <HAP>;
    close(HAP);
    
    my $HAPHeader = $HAPfile[0];
    chomp($HAPHeader);
    my @HAPHeaderArray = split("\t", $HAPHeader);

    for (my $i=1; $i<=$#HAPfile; $i++){ #Iterate over haplotype file lines
        my $line = $HAPfile[$i];
        chomp($line);
        my @array = split("\t", $line);
        my $geneID = $array[3];
        $geneID =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//;
        #$hapGeneIDs{ $geneID } = $geneID;
        my $aCount = $array[4];
        my $bCount = $array[5];
        my $totCount = $array[6];
        #If total coverage is above 20X and both haplotypes have > 0 reads covering it
        if ($totCount >= 20 && $aCount > 0 && $bCount > 0) {
            my $nVariants = $array[8];
            my $variants = $array[9];
            my $sample = $array[13];
            my @allDirAltAlls = split(",", $array[14]);
            my @allHapA = split(",", $array[17]);
            my @allHapB = split(",", $array[18]);
            
            #If gene exists in genes hash continue, otherwise don't process it at all
            #print "Gene: $geneID\n";
            #if (exists $genes{ $geneID }) {
                my @variantsToProcess = split(",", $variants);    
                #Iterate over all variants on a line
                for (my $j=0; $j<=$#variantsToProcess; $j++){
                    my $variant = $variantsToProcess[$j];
                    my $dirAltAll = $allDirAltAlls[$j];
                    #my $varAltAll = $variantAlts[$j];
                    my $allHapA = $allHapA[$j];
                    my $allHapB = $allHapB[$j];
                    my @var = split("_", $variant);
                    my $ref = $var[2];
                    my $alt = $var[3];
                    
                    my $refC=0;
                    my $altC=0;
                    if (exists $binomStatusHash{ $sample }{ $geneID }) {
                        my $binomStatus = $binomStatusHash{ $sample }{ $geneID };
                        if (exists $varRefCount{ $variant }{ $binomStatus }) { #If variant already observed
                            $refC = $varRefCount{ $variant }{ $binomStatus }; #retrieve ref count already present in hash
                            $altC = $varAltCount{ $variant }{ $binomStatus };
                        }
                        
                        my $refCtoAdd;
                        my $altCtoAdd;
                        if ($ref eq $allHapA) { #if reference allele matches hapA allele, add it here, otherwise add it to alt count
                            $refCtoAdd = $refC + $aCount;
                            $altCtoAdd = $altC + $bCount;
                        }else{ #
                            $refCtoAdd = $refC + $bCount;
                            $altCtoAdd = $altC + $aCount;
                        }
                        #Debug
                        #if ($variant eq "22_31371578_C_G"){                    
                        #    print "$variant\t$binomStatus\t$refC\t$refCtoAdd\t$altC\t$altCtoAdd\n";
                        #}
                        
                        #Add new numbers to hashes
                        $varRefCount{ $variant }{ $binomStatus } = $refCtoAdd;
                        $varAltCount{ $variant }{ $binomStatus } = $altCtoAdd;
                        $varStatusCount{ $variant }{ $binomStatus }++;
                        
                        #Check if sample exists in nonOutlierHash, if so add ref and alt counts to new hash.
                        #This new hash will later be used to construct two seperate tables
                        if (exists $nonoutlierHash{ $sample }{ $geneID }) {
                            if ($ref eq $allHapA) { #if reference allele matches hapA allele, add it here, otherwise add it to alt count
                                $varSampleRefCount{ $variant }{ $sample } = $aCount;
                                $varSampleAltCount{ $variant }{ $sample } = $bCount;
                            }else{
                                $varSampleRefCount{ $variant }{ $sample } = $bCount;
                                $varSampleAltCount{ $variant }{ $sample } = $aCount;
                            }
                        }
                    }
                }
            #}
        }
    }

    print "Done processing haplotypic count files\n\n";
    
    
    print "\nGenerating output file ..\n";
    
    foreach my $key (sort keys %annotations){
        my $variant = $key; #Variant
        my $af = $afs{ $variant }; #Allele frequency from data
        my $exacAF = $exacs{ $variant };
        my $gonlAF = $gonls{ $variant };
        my $clinvar = "NA";
        if (exists $clinvar{ $variant }) { #Extract Clinvar status
            $clinvar = $clinvar{ $variant };
        }
        my $annotation = $annotations{ $variant }; #Row of annotations
        my $pLI = "NA";
        my $pliStat = "NA";
        #Set pLI status
        if (exists $plis{ $variant }) {
            $pLI = $plis{ $variant }; #pLI score
            if (looks_like_number($pLI)) {
                if ($pLI >= 0.9) {
                    $pliStat = "INTOLERANT";
                }elsif ($pLI <= 0.1){
                    $pliStat = "TOLERANT";
                }else{
                    $pliStat = "NEITHER";
                }
            }else{
                $pliStat = "NA";
            }
        }
        #Set CGD status
        my $geneID = $varToGenes{ $key };
        my $inheritance = "NA";
        my $ageGroup = "NA";
        my $manifestation = "NA";
        my $geneName = "NA";
        #print "geneID; $geneID\n";
        if (exists $varToGenes{ $key }) {
            $geneName = $varToGenes{ $key };
            #print "geneName: $geneName\n";
            if (exists $CGDgenes{ $geneName }) {
                $inheritance = $inheritance{ $geneName };
                $ageGroup = $ageGroup{ $geneName };
                $ageGroup =~ s|N/A|NA|g;
                $manifestation = $manifestation{ $geneName };
            }
        }
        if (exists $varRefCount{ $variant }{ YES }) { #Counts for this variant exist, so print it to output
            my $nonoutlierVarRefCount = "NA";
            my $nonoutlierVarAltCount = "NA";
            my $nonoutlierVarStatusCount = "NA";
            if (exists $varStatusCount{ $variant }{ NO }) {
                $nonoutlierVarRefCount = $varRefCount{ $variant }{ NO };
                $nonoutlierVarAltCount = $varAltCount{ $variant }{ NO };
                $nonoutlierVarStatusCount = $varStatusCount{ $variant }{ NO };
            }
            my $outlierVarRefCount = $varRefCount{ $variant }{ YES };
            my $outlierVarAltCount = $varAltCount{ $variant }{ YES };
            my $outlierVarStatusCount = $varStatusCount{ $variant }{ YES };
            
            #If AF is > 0.5, swap counts
            my $updNonoutlierVarRefCount = $nonoutlierVarRefCount;
            my $updNonoutlierVarAltCount = $nonoutlierVarAltCount;
            my $updOutlierVarRefCount = $outlierVarRefCount;
            my $updOutlierVarAltCount = $outlierVarAltCount;
            if ($af > 0.5) {
                $updNonoutlierVarRefCount = $nonoutlierVarAltCount;
                $updNonoutlierVarAltCount = $nonoutlierVarRefCount;
                $updOutlierVarRefCount = $outlierVarAltCount;
                $updOutlierVarAltCount = $outlierVarRefCount;
            }
            
        }
    }
    print "Done generating output file for chr$CHR.\n\n";

}

#Create matrices containing counts for major and minor alleles
print "Generating output major/minor allele count files ..\n";
open(MAJOR, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.majorAllelle.nonASEsamples.chrALL.txt") || die "Can't open output file for major allele!\n";
open(MINOR, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.minorAllelle.nonASEsamples.chrALL.txt") || die "Can't open output file for minor allele!\n";
shift(@BINOMheader);
my $header = join("\t", @BINOMheader);
chomp($header);
print MAJOR "VARIANT\t$header\n";
print MINOR "VARIANT\t$header\n";
foreach my $key (sort keys %annotations){
    my $variant = $key; #Variant
    print MAJOR "$variant";
    print MINOR "$variant";
    for (my $s=0; $s<=$#BINOMheader; $s++){
        my $sample = $BINOMheader[$s];
        my $refC = "NA";
        my $altC = "NA";
        my $af = $afs{ $variant }; #Allele frequency from data
        if (exists $varSampleRefCount{ $variant }{ $sample }) {
            if ($af > 0.5) { # Allele frequency above 0.5, swap major and minor allele
                $refC = $varSampleAltCount{ $variant }{ $sample };
                $altC = $varSampleRefCount{ $variant }{ $sample };
            }else{
                $refC = $varSampleRefCount{ $variant }{ $sample };
                $altC = $varSampleAltCount{ $variant }{ $sample };
            }
        }
        print MAJOR "\t$refC";
        print MINOR "\t$altC";
    }
    print MAJOR "\n";
    print MINOR "\n";
}
close(MAJOR);
close(MINOR);
print "Done generating output major/minor allele count files\n";
