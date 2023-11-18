#  Author:  Yoshinori Fukasawa
#  Organizations:  Department of Computational Biology, University of Tokyo
#  Copyright (C) 2014, Yoshinori Fukasawa, All rights reserved.
#  Creation Date:  2014/7/8

package DirichletRegulator;
use Exporter;
@ISA = (Exporter);
@EXPORT = qw(DirichletPseudoCount DirichletParser);
use strict;
use warnings;

sub DirichletPseudoCount
{
  ## Arguments
  ## @vec : Count Vector $#vec == $AlphabetSize. This must be given by user(s). AA frequencies at each columns.
  ## @eq  : Prior Mixture Prob. A.K.A coefficients of mixture distribution.
  ## @e   : An array which stores alpha terms in each Dirichlet components. @e[$numDirichletComponents][$AlphabetSize]
  ##
  my ($vecRef, $Struct) = @_;
  my @vec = @{$vecRef};
  my @mix; ## Posterior Proberbility
  my $AlphabetSize = 20; ## constant

  ## Mixture coefficients of default Dirichlet mixture 9
  my @eq = @{$Struct->{eq}};
  my @e = @{$Struct->{e}};

  ## Calculating Positerior Probability of each component given counted vector.
  ## $q is index for Dirichlet Comp.
  my $q;
  my $numDirichletComponents = $#e+1;

  $mix[0] = 1;
  for($q = 0; $q < $numDirichletComponents; $q++){
    $mix[$q] = $eq[$q] > 0 ? log($eq[$q]) : -999 ;
    $mix[$q] += &Logp_cvec($vecRef,$AlphabetSize,$e[$q]);
  }

  #$NormMix[$q], posterior prob. of each comp., is P(component_q | n)
  my @NormMix = &LogNorm(\@mix, $numDirichletComponents);

  # Convert the counts to prob., following Sjolander(1996)
  my $totc = &Sum($vecRef);
  for(my $i = 0; $i < $AlphabetSize; $i++){
    my $xi = 0;
    for($q = 0; $q < $numDirichletComponents; $q++){
      my $tota = &Sum($e[$q]);
      $xi += $NormMix[$q] * ($vec[$i] + $e[$q][$i]) / ($totc + $tota);
    }
    $vec[$i] = $xi;
  }

  &Norm($vecRef, $AlphabetSize);
  return \@vec;
}

# Sub-routine: prob2score and score2prob
#
# Purpose: Each routines return score and probability, respectively.
# Formula is following HMMER 2.3.2
# These routines used for test

sub prob2score
{
  my $score = shift;
  my $nullProb = shift;

  return $score =~ /\*/ ? 0 : $nullProb*exp(log(2)*$score/1000);
}

sub score2prob
{
  my $prob = shift;       # convert this probability to a score
  my $nullProb = shift;   # given this null probability

  return $prob != 0 ? int( 0.5+1000*log($prob/$nullProb)/log(2.0) ) : "*";
}

# Sub-routine: Logp_cvec()
#
# Purpose:  Calculates ln P(cvec|dirichlet), the log probability of a
#           count vector given a Dirichlet distribution. Adapted
#           from an implementation by Graeme Mitchison.
#
# Args:     cvec  - count vector
#            n     - length of cvec
#           alpha - Dirichlet alpha terms
#
# Return:   log P(cvec|dirichlet)
# From: math_support.c

sub Logp_cvec
{
  my ($cvecRef, $AlphabetSize, $alphaRef) = @_;
  my @cvec = @{$cvecRef};
  my @alpha = @{$alphaRef};

  my $lnp; #log likelihood of P(cvec | Dirichlet)
  my ($sum1,$sum2,$sum3);
  my $i;

  $lnp = $sum1 = $sum2 = $sum3 = 0;
  for($i = 0; $i < $AlphabetSize; $i++){
    $sum1 += $cvec[$i] + $alpha[$i];
    $sum2 += $alpha[$i];
    $sum3 += $cvec[$i];
    $lnp += &Gammln($alpha[$i] + $cvec[$i]);
    $lnp -= &Gammln($cvec[$i] + 1);
    $lnp -= &Gammln($alpha[$i]);
  }
  $lnp -= &Gammln($sum1);
  $lnp += &Gammln($sum2);
  $lnp += &Gammln($sum3 + 1);
  return $lnp;
}


# Sub-routine: LogNorm()
#
# Purpose:  Normalize a vector of log likelihoods, changing it
#           to a probability vector. Be careful of overflowing exp().
#           Implementation adapted from Graeme Mitchison.
#
# Args:     posterior - vector destined to become log probabilities
#           n   - length of vec
#
# Return: Normalized Posterior Proberbilities
# From: mathsupport.c

sub LogNorm
{
  if(@_ != 2){
    die "Argument Number is weird.\n";
  }
  my @posterior;
  my ($tempRef,$n) = @_;
  @posterior = @{$tempRef};

  my ($i,$max);
  $max = -1E+30;
  my $denom = 0;

  for( $i = 0; $i < $n; $i++){
    $max = $posterior[$i] if($posterior[$i] > $max);
  }
  # Calculating Summation of all log(coefficient * Z(n+alpha)/Z(alpha))
  for( $i = 0; $i < $n; $i++){
    $denom += exp($posterior[$i] - $max)if($posterior[$i] > $max -50);
  }
  for( $i = 0; $i < $n; $i++){
    if($posterior[$i] > $max -50){
      $posterior[$i] = exp($posterior[$i] - $max) / $denom;
    } else {
      $posterior[$i] = 0;
    }
  }
  return @posterior;
}

# Sub-routine: Sum()
#
# Args:    vec - vector array
#          n - length of vector
#
# Return:  Summation of given vector
# From: vectorops.c

sub Sum
{
  my $ref = shift;
  my @vec = @{$ref};
  my $sum = 0;
  for(my $i = 0; $i < @vec ; $i++){
    $sum += $vec[$i];
  }
  return $sum;
}

# Sub-routine: Norm()
# Args:    refVec - reference to count vector array
#          n - length of vector
#
# Return:  Summation of given vector
# From: vectorops.c

sub Norm
{
  my $refVec; #Reference
  my $n;
  my $sum;
  ($refVec,$n) = @_;
  $sum = &Sum($refVec);

  if($sum != 0){
    for(my $i = 0; $i < $n ; $i++){
      $$refVec[$i] /= $sum;
    }
  } else {
    for(my $i = 0; $i < $n ; $i++){
      $$refVec[$i] = 1 / $n;
    }
  }
}

# Sub-routine: Gammln()
#
# Returns the natural log of the gamma function of $x.
# $x is > 0.0.
#
# Adapted from a public domain implementation in the
# NCBI core math library. Thanks to John Spouge and
# the NCBI. (According to the NCBI, that's Dr. John
# "Gammas Galore" Spouge to you, pal.)
# From: sre_math.c

sub Gammln
{
  my $x = shift;
  my $i;
  my ($xx, $tx);
  my ($tmp, $value);
  my @cof = (
    4.694580336184385E+04,
    -1.560605207784446E+05,
    2.065049568014106E+05,
    -1.388934775095388E+05,
    5.031796415085709E+04,
    -9.601592329182778E+03,
    8.785855930895250E+02,
    -3.155153906098611E+01,
    2.908143421162229E-01,
    -2.319827630494973E-04,
    1.251639670050933E-10
  );

  # Protect against $x=0. We see this in Dirichlet code,
  # for terms alpha = 0. This is a severe hack but it is effective
  # and (we think?) safe. (due to GJM)

  return 999999 if ($x <= 0);

  $xx = $x - 1;
  $tx = $tmp = $xx + 11;
  $value    = 1;
  for ($i = 10; $i >= 0; $i--){ # sum least significant terms first
    $value += $cof[$i] / $tmp;
    $tmp   -= 1;
  }
  $value  = log($value);
  $tx    += 0.5;
  $value += 0.918938533 + ($xx+0.5)*log($tx) - $tx;
  return $value;
}

# Sub-routin: DirichletParser
# Return: Structure of Dirichlet Mixture
# Argument: filepath to open SAM developer's Dirichlet Ascii file


sub DirichletParser
{
  my $filepath = shift;
  open(my $fileHandle, "<", $filepath);
  my @FILE = <$fileHandle>;
  close $fileHandle;

  my $NumOfComp;
  my $name;
  my $compCounter=-1;
  my @eq;
  my @e;
  my $struct = {
      name => '',
      eq => \@eq,
      e => \@e
  };

  foreach(@FILE){
    my $string;
    chomp($_);
    if($_ =~ /^Name = (.+)$/){
      my $temp = $1;
      $name = $temp;
      $struct->{name} = $name;
    } elsif($_ =~ /^NumDistr= ([0-9]+)/){
      $NumOfComp = $1;
    } elsif($_ =~ /^Number=/){
      $compCounter++;
    } elsif($_ =~ /^Mixture=/){
      $string = &substitution("Mixture", '', $_);
      push(@{$struct->{eq}}, $string);
    } elsif($_ =~ /^Alpha=/){
      $string = &substitution("Alpha", '', $_);
      my @array = split(/\s/, $string);
      shift @array; #first value is summation. sigma = alpha_i
      $struct->{e}->[$compCounter] = \@array;
    } elsif($_ =~ /^Comment=/){
      $string = &substitution("Comment", "\# ", $_);
    } elsif($_ eq ''){
    } elsif($_ =~ /^EndClassName/ && $compCounter+1 != $NumOfComp){
      print $compCounter, "\n";
      print "Error. Making construct failed.\n";
    }
  }
  return $struct;
}

sub substitution
{
  (my $target, my $tr, my $string) = @_;
  $string =~ s/$target//;
  if($string !~ /^= .+/){
    print "line $string doesn't match with $target.\n ";
  } else {
    $string =~ s/^= //;
  }
  $string = "$tr"."$string";
  return $string;
}

1;
