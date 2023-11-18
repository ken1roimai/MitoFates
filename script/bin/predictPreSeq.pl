#!/usr/bin/perl -w
use strict;
use warnings;
use File::Temp qw/ tempfile tempdir /;

BEGIN{
    use FindBin qw ($Bin);
    push @INC, "$FindBin::Bin/bin";
}

if (@ARGV>0) {
    open(IN,"$ARGV[0]")||die "Cannot open $ARGV[0] \n";
}else {
    die "Usage: [predictPreSeq.pl] [MultiFastaFile]\n";
}

my @ids;
while(<IN>){
    chomp($_);
    if(/^>/){
	$_ =~ s/>//g;
	push(@ids,$_);
    }
}

# making a temporary directory
my $t_dir = tempdir( CLEANUP=> 1);
# making a temp file for feature vectors
my ($fh_svm, $svm_file) = tempfile( DIR => $t_dir, SUFFIX => '.svm');
# making a temp file for scaled feat vec file
my ($fh_scl, $scaled_file) = tempfile( DIR => $t_dir, SUFFIX => '.scaled');
# making a temp file for a result
my ($fh_res, $result_file) = tempfile( DIR => $t_dir, SUFFIX => '.result');

#Get SVM features
system("perl $Bin/computeSVMFeatureForPresequence.pl -c 1 $ARGV[0] > $svm_file");

#scaling
system("svm-scale -r $Bin/PreSeq.scaling $svm_file > $scaled_file");

#Predict Preseq
`svm-predict -b 1 $scaled_file $Bin/PreSeq.model $result_file`;

open(RS,"$result_file")||die "Cannot open input.result \n";
my @result = <RS>;
for(0..$#ids){
    my $id = $ids[$_];
    my ($label,$prob1,$prob2) = split(/\s/,$result[$_+1]);
    my $tag;
    if($prob1 >= 0.385){
	$tag = "Possessing mitochondrial presequence";
    }else{
	$tag = "No mitochondrial presequence";
    }
    print "$id\t$prob1\t$tag\n";
}
