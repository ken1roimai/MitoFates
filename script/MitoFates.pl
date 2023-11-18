#!/usr/bin/perl
#  Author:  Yoshinori Fukasawa and Kenichiro Imai
#  Organizations:  Department of Computational Biology, University of Tokyo
#  Copyright (C) 2014, Yoshinori Fukasawa and Kenichiro Imai, All rights reserved.

BEGIN{
    use FindBin qw ($Bin);
    use lib "$Bin/bin/modules";
}

use strict;
use warnings;
use MotifPosition;
use MotifHelixPosition;

my $VERSION = "v1.1";

######### Initialize variables ##########
my $scriptDir = $Bin;
####################

if (@ARGV<2){
    die "Usage: [MitoFates.pl] [MultiFastaFile] [Organism Flag: fungi, metazoa or plant]\n";
}

my $osFlag = "$ARGV[1]";

####################
my @Tom20Motif_array;
my @Helix_array;
my @MTSMotif_array;
#to print out motifs in the header, sorted motif list is required.
my @motifs = sort(MotifHelixPosition->new()->getMotifs());
#################### 

# Get Positions for seqs.
open my $fastain, "<", $ARGV[0] || die "Given fasta file cannot be opend";
{
    local $/ ="\n>";

    while(my $line = <$fastain>){

    # Turn $/ into original.
	local $/ = "\n";

	chomp($line);
	my @head_seq = split(/\n/, $line);
	my $id = shift @head_seq;
	my $sequence = join("", @head_seq);
	chomp($sequence);
	$id =~ s/>//g; #for the first header

	my @posArray = findMotifPosition(Seq => $sequence, Length=>100, Pattern=>"pcbp2");
	push @Tom20Motif_array, \@posArray;

	my $searchSpace = 30;
	my $mhp = MotifHelixPosition->new();
	$mhp->calcHmoment(substr($sequence,0, $searchSpace), 96, 8.5);
	push @Helix_array, join(
	    "-",
	    $mhp->getPos,
	    $mhp->getPos+$mhp->getWindowSize-1,
	    $mhp->getMoment >= 2 ? "high" : "low"
	    );
	my $ref = $mhp->searchAndGetPositions($sequence, $searchSpace);

	if($ref){
	    push @MTSMotif_array, $ref;
	}
	##### <<

    }
}

## Prediction of Presequence
my @resultPS = `perl $Bin/bin/predictPreSeq.pl $ARGV[0]`;


## Prediction of CleavageSite
my @resultCS = `perl $Bin/bin/cleavage.pl --gamma --svm --$osFlag $ARGV[0]`;

## Print header
print "Sequence ID\tProbability of presequence\tPrediction\tCleavage site (processing enzyme)\tNet charge\tPositions for TOM20 recognition motif (deliminated by comma)\tPosition of amphypathic alpha-helix\t";

print join("\t", @motifs);
print "\n";

for(0..$#resultPS){
    chomp($resultPS[$_]);
    chomp($resultCS[$_]);
    my ($id,$PreseqProb,$preSeqLabel) = split(/\t/,$resultPS[$_]);
    my ($id2, $mppProb, $preposi, $netcharge, $fragment, $oct1PWMscore, $Icp55PWMscore) = split(/\t/,$resultCS[$_]);

    # Prediction result
    printf "$id\t%0.3f\t$preSeqLabel\t$preposi\t%0.3f\t",$PreseqProb, $netcharge;
    # Tom20 motif
    print join(",", @{$Tom20Motif_array[$_]}), "\t";
    # helix
    print $Helix_array[$_], "\t";
    # MTS motifs
    printMTSMotif($MTSMotif_array[$_]);
    print "\n";
}

sub printMTSMotif
{
    if(@_ != 1){
	print STDERR "\tError\n";
	return 0;
    }

    my $ref = shift;
    foreach my $motif (sort keys %{$ref}){
	if(@{$ref->{$motif}}){
	    print join(",", @{$ref->{$motif}});
	} else {
	    print "-";
	}
	print "\t";
    }
}
