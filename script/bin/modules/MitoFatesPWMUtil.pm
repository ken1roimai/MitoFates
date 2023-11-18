#  Author:  Yoshinori Fukasawa
#  Organizations:  Department of Computational Biology, University of Tokyo
#  Copyright (C) 2011, Yoshinori Fukasawa, All rights reserved.
#  Creation Date:  2011/10/11

package MitoFatesPWMUtil;
require Exporter;
our @ISA = qw (Exporter);
our @EXPORT = qw (ret_gammaParams Log2 seq GammaWeighting findMax cutter netcharge chargedResidue Calc_mean_hydrophobicity
                  Calc_aa_compositions mean stddev charge_around_Cleavage
                  Hmoment Calc_aa_composition_DirichletPosterior);
BEGIN{
    use FindBin qw ($Bin);
    use lib "$Bin";
}

use strict;
use warnings;
use diagnostics;
use Data::Dumper;
use Math::Cephes qw(gamma gdtr);
use DirichletRegulator_fast;


# Some mathmatical sub-routines ######################
## Ceiling Function
## Argument:
## $X value
sub CEILINGFUNCTION{
    my ($X) = @_;
    my $CeilingFunction = 0;
    my $Temp = 0;

    if($X == int($X)){
        $Temp = $X;
    }else {
        if($X > 0){
            $Temp = int($X) + 1;
        }else {
            $Temp = int($X);
        }
    }
    $CeilingFunction = $Temp;
    return $CeilingFunction;
}

sub HarmonicAverage{
  my ($a, $b) = @_;
  my $rtnValue = (2*$a*$b)/($a+$b);
  return($rtnValue);
}

sub Log2{
    my $x = shift;
    return log($x) / log(2);
}

sub mean{
    my @values = @_;
    return 0 unless(@values);
    my $x=0;
    foreach my $value (@values){
        $x += $value;
    }
    return $x/(scalar @values);
}

sub stddev{
    my @values = @_;
    return 0 unless(@values);
    my $mean = &mean(@values);
    my $x=0;
    foreach my $xi (@values){
        $x += ($xi - $mean)**2;
    }

    return sqrt($x/((scalar @values)-1));
}

###################################################################################################

# Miscellaneous sub-routines #########################
#

# Return a list which containing each windows
sub cutter{
    my %args = (
        sequence=> "",
        length => 100,
        windowSize=>10,
        @_
    );
    die "==\tGive a sequence.\n" unless(exists($args{sequence}));
    my $len = length($args{sequence});
    my @array;
    my $threshold = $args{length}-$args{windowSize}+2;
    for(my $i = 0; $i < $threshold; $i++){
        if($i <= $len-$args{windowSize}){
            if($i <= $args{length}-$args{windowSize}){
                my $negative = substr($args{sequence}, $i, $args{windowSize});
                $array[$i] = $negative;
            } elsif ($i > $len-$args{windowSize}){
                last;
            }
        }
    }
    return @array;
}

## Calculate 20 amino acid compositions
## Arguments:
## seq return 20 compositions from this gived sequence
sub Calc_aa_compositions{
    my $seq = shift;
    my $length = length($seq);
    my %hash;
    my $total=0;
    my @AA = (
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
        'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
    );

    foreach my $amino (@AA){
        my $temp = grep (/$amino/, split(//, $seq));
        $hash{$amino} = $temp/$length;
        $total += $temp;
    }

    if($length == $total){
        return \%hash;
    }
    else {
        print STDERR "Error Length: $length, Counted: $total.\n";
        print STDERR "Seems to include weird amino acid character.\n";
    }
}

sub Calc_aa_composition_DirichletPosterior{
    if(@_ != 2){
        print STDERR "\tError occured at Calc_aa_composition_DirichletPosterior\n";
        print STDERR "\tThe number of arguments is incorrect\n";
        exit(1);
    }

    my $seq = shift;
    my $struct = shift;
    my $length = length($seq);
    my %hash;
    my $total=0;
    my @AA = (
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
        'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
    );

    $DB::single = 1;

    my @aacounts;
    foreach my $amino (@AA){
        my $temp = grep (/$amino/, split(//, $seq));
        push @aacounts, $temp;
        $total += $temp;
    }

    $DB::single = 1;

    my $pseudoCounts = DirichletPseudoCount(\@aacounts, $struct);

    if($length == $total){
        return $pseudoCounts;
    }
    else {
        print STDERR "Error Length: $length, Counted: $total\n";
    }
}

## Similar function of "seq" in R language.
## But code iself is not copied from R. Need more sophisticated implementation 
sub seq{
    my ($start, $end, $by) = @_;
    my @array;
    for(my $i=$start; $i<$end; $i+=$by){
        push @array, $i;
    }
    return @array;
}

## This routine searches maximum and second maximum values.
## Arguments: Score array to be searched
sub findMax{
    my $maxsite_pos = 0;
    my $secondsite_pos = -1;
    my @score_array = @_;

    for(my $i=0; $i<@score_array; $i++){
        if($score_array[$maxsite_pos]<$score_array[$i]){
            $maxsite_pos = $i;
        }
    }

    for(my $i=0; $i<@score_array; $i++){
        if($score_array[$maxsite_pos] < $score_array[$i]){
            $maxsite_pos = $i;
        }
        if($score_array[$maxsite_pos] > $score_array[$i]){
            if($secondsite_pos != -1){
                $secondsite_pos = $i if($score_array[$secondsite_pos] < $score_array[$i]);
            } else {
                $secondsite_pos = $i;
            }
        }
    }

    return($maxsite_pos, $secondsite_pos);
}

###################################################################################################

# Sub-routines related to Gamma distribution #########
#

## Return a construct which contains parameters for Gamma mixture
sub ret_gammaParams{
    my $filePath = shift;
    my @params;
    open my $handle, "<" ,$filePath || die "Cannot open Gammma Distribution file: $!\n";
    my @lines = <$handle>; #each line contains attribute, shape and scale.
    close $handle;
    for(my $i=0;$i<@lines;$i++){
        my %array;
        chomp $lines[$i];
        my @temp = split(/\t/, $lines[$i]);
        $array{"lambda"} = $temp[0]; # lambda
        $array{"shape"} = $temp[1]; # shape
        $array{"scale"} = $temp[2]; # scale
        $params[$i] = \%array;
    }
    return @params;
}

## Calculate Gamma Distribution's not cdf but pdf
## Math::Cephes requires rate parameter, but we have scale param.
## These routines are wrapper to avoid bugs as to the above difference.

## Is this required in this script?
sub gammaPDF{
    my ($x, $shape, $scale) = @_;
    my $beta = 1/$scale;
    my $y = ($beta**$shape/gamma($shape))*($x**($shape-1)*exp(-1*$beta*$x));
    return $y;
}

## Calculate Gamma Distribution's CDF
sub gammaCDF{
    my ($x, $shape, $scale) = @_;
    my $beta = 1/$scale; #a.k.a rate parameter
    my $y = gdtr($x, $shape, $beta);
    return $y;
}

###################################################################################################

# Sub-routines for weighting by given Gamma mixture###
#

## To weight logarithmic ratio of hmmer by gamma mixture's pdf
## Arguments:
## position same as the above
## scoreRef reference to hmmer score

## Define P($position) = \int_{$position-1}^{$position} f_{X}(x)dx
## => F_{X}($position) - F_{X}($position-1)
sub GammaWeighting{
    my ($position, $scoreRef, $paramsRef) = @_;
    my $score = $$scoreRef;
    my @params = @{$paramsRef};

    # Calculate probability
    my $y=0;
    for(my $i=0; $i<@params; $i++){
        $y += $params[$i]->{"lambda"}*(
          gammaCDF($position, $params[$i]->{"shape"}, $params[$i]->{"scale"})
              -gammaCDF($position-1, $params[$i]->{"shape"}, $params[$i]->{"scale"})
          );
    }
    $$scoreRef += &Log2($y);
}

###################################################################################################

###################################################################################################

# Sub-routines for physico-chemical features #########
#

## Calculate averaged hydrophobicity by KD indices
## Arguments:
## seq a sequence to be calclated its hydrophobicity
sub Calc_mean_hydrophobicity{
    my %args = (
        type => "KTDL",
        sequence => "",
        @_
    );

    my $KTDL = {
        # Kyte, J. & R. F. Doolittle. 1982. J Mol Biol 157:105-132.
        'S' => '-0.8', 'F' => '2.8', 'T' => '-0.7',
        'N' => '-3.5', 'K' => '-3.9', 'Y' => '-1.3',
        'E' => '-3.5', 'V' => '4.2', 'Q' => '-3.5',
        'M' => '1.9', 'C' => '2.5', 'L' => '3.8',
        'A' => '1.8', 'W' => '-0.9', 'P' => '-1.6',
        'H' => '-3.2', 'D' => '-3.5', 'R' => '-4.5',
        'I' => '4.5', 'G' => '-0.4',
    };
    my $ABODR = {
        # Aboderin. 1971.Int J Biochem 2:537-544.
        # J. CORNETTE, et al. 1987. J Mol Biol 195:659-685
        'S' => '-0.83', 'F' => '4.56', 'T' => '-0.50',
        'N' => '-2.90', 'K' => '-2.32', 'Y' => '3.23',
        'E' => '-1.91', 'V' => '3.65', 'Q' => '-2.24',
        'M' => '3.81', 'C' => '0.00', 'L' => '4.89',
        'A' => '0.83', 'W' => '4.23', 'P' => '0.66',
        'H' => '-2.07', 'D' => '-2.82', 'R' => '-1.74',
        'I' => '4.31', 'G' => '0.00',
    };

    unless(exists($args{sequence})){
        print STDERR "==\tCan't calculate hydrophobicity. Please give a sequence.\n";
        return 0;
    }
    my %hydrophobicity;
    if($args{type} eq "KTDL"){
        %hydrophobicity = %{$KTDL};
    } elsif($args{type} eq "ABODR"){
        %hydrophobicity = %{$ABODR};
    } else {
        print STDERR "==\tCan't calculate hydrophobicity. Please specify the index of hydrophobicity.\n";
        return 0;
    }

    chomp($args{sequence});
    my $length = length($args{sequence});
    my @Seq_splitted = split(//, $args{sequence});
    my $total=0;
    foreach my $amino (@Seq_splitted){
        if($hydrophobicity{$amino}){
            $total+= $hydrophobicity{$amino};
        } else {
            print STDERR "Undefined amino acid is used: $amino\n";
        }
    }

    return $total/$length;
}

## Calculate averaged net charge
## Arguments:
## seq_temp a query sequence
## clsite position of (predicted) cleavage site, which starts from 0 like an array
sub netcharge{
    my ($seq_temp, $clsite) = @_; # +1: clevage position from zero start.
    chomp($seq_temp);
    my $preseq = substr( $seq_temp, 0, $clsite+1 );
    my @array = split(//, $preseq);
    my $charge = 0;
    foreach my $aa (@array){
        if($aa eq "D" || $aa eq "E"){
            $charge = $charge - 1;
        } elsif($aa eq "R" || $aa eq "K"){
            $charge = $charge + 1;
        }
    }
    $charge = $charge/length($preseq);
    return($charge);
}

## To return number of charged residues, counts them
## Arguments:
## seq_temp same as the above
## clsite same as the above
sub chargedResidue{
    my ($seq_temp, $clsite) = @_;
    chomp($seq_temp);
    my $preseq = substr($seq_temp, 0, $clsite+1); # +1: clevage position from zero start point.
    my @Pre = split(//, $preseq);
    my @ChaResCount;
    $ChaResCount[0] = grep /R/, @Pre;
    $ChaResCount[1] = grep /K/, @Pre;
    $ChaResCount[2] = grep /H/, @Pre;
    $ChaResCount[3] = grep /D/, @Pre;
    $ChaResCount[4] = grep /E/, @Pre;

    for(my $i=0; $i<5; $i++){
        $ChaResCount[$i] /= ($clsite+1);
    }

    return(@ChaResCount);

}
######################################################

1;
