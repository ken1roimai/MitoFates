#  Author:  Yoshinori Fukasawa
#  Organizations:  CBRC, AIST
#  Copyright (C) 2014, Yoshinori Fukasawa, All rights reserved.
#  Creation Date:  2014/10/22

package DirichletRegulator_fast;
use Exporter;
@ISA = (Exporter);
@EXPORT = qw(DirichletPseudoCount DirichletParser);
use strict;
use warnings;
use FindBin qw($Bin);
use Inline (Config => 
	      DIRECTORY => "$Bin/modules/_Inline/",
            );
use Inline 'C';

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

__DATA__
__C__

// To pass arrayref is smart enough?
double
_Sum(SV *ref)
{
  AV *array;
  SV **value;
  double sum;
  int i;

  array = (AV *)SvRV(ref);

  sum = 0;
  for(i = 0; i < av_len(array)+1 ; i++){
    value = av_fetch(array, i, 0);
    sum += SvNV(*value);
  }

  return sum;
}

double
_Gammln(double x)
{
  int i;
  double xx, tx;
  double tmp, value;
  double cof[] = {
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
  };

  if (x <= 0)
    return 999999;

  xx = x - 1;
  tx = tmp = xx + 11;
  value    = 1;
  for (i = 10; i >= 0; i--){
    // sum least significant terms first
    value += cof[i] / tmp;
    tmp   -= 1;
  }
  value  = log(value);
  tx    += 0.5;
  value += 0.918938533 + (xx+0.5)*log(tx) - tx;
  return value;

}

double
_Logp_cvec(SV *cvec, int D, SV *alpha)
{
  AV *cvec_a  = (AV*) SvRV(cvec);
  AV *alpha_a = (AV*) SvRV(alpha);

  double lnp; // log likelihood of P(cvec | Dirichlet)
  double cvec_value, alpha_value;
  double sum1,sum2,sum3;
  int i;

  lnp = sum1 = sum2 = sum3 = cvec_value = alpha_value = 0;

  for(i = 0; i < D; i++){
    cvec_value  = SvNV(*av_fetch(cvec_a, i, 0));
    alpha_value = SvNV(*av_fetch(alpha_a, i,0));
    sum1 += cvec_value + alpha_value;
    sum2 += alpha_value;
    sum3 += cvec_value;
    lnp += _Gammln(cvec_value + alpha_value);
    lnp -= _Gammln(cvec_value + 1.0);
    lnp -= _Gammln(alpha_value);
  }

  lnp -= _Gammln(sum1);
  lnp += _Gammln(sum2);
  lnp += _Gammln(sum3 + 1.0);

  return lnp;
}

void
_Norm(SV *array, int n)
{
  double sum;
  sum = _Sum(array);

  AV *a = (AV*) SvRV(array);

  int i;
  double temp;
  if(sum != 0){
    for(i = 0; i < n ; i++){
      temp = SvNV(*av_fetch(a,i,0));
      temp /= sum;
      sv_setnv(*av_fetch(a,i,TRUE), temp);
    }
  } else {
    for(i = 0; i < n ; i++){
      temp = SvNV(*av_fetch(a,i,0));
      temp /= (double)n;
      sv_setnv(*av_fetch(a,i,TRUE), temp);
    }
  }
}

void
_LogNorm(SV *post, int n)
{
  AV *posterior = (AV*) SvRV(post);

  int i;
  double max = -1E+30;
  double denom = 0;

  for(i = 0; i<n; i++){
    if(SvNV(*av_fetch(posterior, i, 0)) > max)
      max = SvNV(*av_fetch(posterior, i, 0));
  }

  // Calculating Summation of all log(coefficient * Z(n+alpha)/Z(alpha))
  for(i = 0; i<n; i++){
    if(SvNV(*av_fetch(posterior, i, 0)) > max - 50)
      denom += exp(SvNV(*av_fetch(posterior, i, 0)) - max);
  }

  for(i = 0; i<n; i++){
    if(SvNV(*av_fetch(posterior, i, 0)) > max -50){
      double temp = exp(SvNV(*av_fetch(posterior, i, 0)) - max) / denom;
      sv_setnv(*av_fetch(posterior,i,TRUE), temp);
    } else {
      sv_setnv(*av_fetch(posterior,i,TRUE), (double)0);
    }
  }

  post = newRV_noinc((SV *)posterior);
}

// Below function intends to be called by only native DirichletPseudoCount().
// Purpose: to avoid many conversion of double[] to AV. waste of memory and time.
// for convinience, keep_LogNorm for calling from Perl

static void
_c_LogNorm(double *post, int n)
{
  int i;
  double max = -1E+30;
  double denom = 0;

  for(i = 0; i<n; i++){
    if(post[i] > max)
      max = post[i];
  }

  // Calculating Summation of all log(coefficient * Z(n+alpha)/Z(alpha))
  for(i = 0; i<n; i++){
    if(post[i] > max - 50)
      denom += exp(post[i] - max);
  }

  for(i = 0; i<n; i++){
    if(post[i] > max -50){
      post[i] = exp(post[i] - max) / denom;
    } else {
      post[i] = 0.0;
    }
  }
}

SV*
DirichletPseudoCount(SV *vec, SV *st)
{
  AV *a = (AV*) SvRV(vec);
  HV *b = (HV*) SvRV(st);

  int AlphabetSize = 20; // constant

  // Mixture coefficients of default Dirichlet mixture 9
  SV** eq_ref = hv_fetch(b, "eq", strlen("eq"),0);
  SV** e_ref  = hv_fetch(b, "e", strlen("e"),0);

  AV* eq = (AV*) SvRV(*eq_ref);
  AV* e  = (AV*) SvRV(*e_ref);

  // Calculating Positerior Probability of each component given counted vector.
  // q is index for Dirichlet Comp.
  int q;
  int numDirichletComponents =  av_len(e)+1;

  double mix[numDirichletComponents]; // Posterior Proberbility

  mix[0] = 1.0;
  for(q = 0; q < numDirichletComponents; q++){
    mix[q] = SvNV(*av_fetch(eq, q,0)) > 0 ? log(SvNV(*av_fetch(eq, q,0))) : -999 ;
    mix[q] += _Logp_cvec(vec, AlphabetSize, *av_fetch(e, q, 0) );
  }

  // NormMix[q], posterior prob. of each comp., is P(component_q | n)
  // Abandoned NormMix. Instead, Normalize an argument directly
  _c_LogNorm(mix, numDirichletComponents);

  // Convert the counts to prob., following Sjolander(1996)
  double totc = _Sum(vec);
  int i;
  for(i = 0; i<AlphabetSize; i++){
    double xi = 0;
    for(q = 0; q<numDirichletComponents; q++){
      double tota = _Sum(*av_fetch(e,q,0));
      AV *elm = (AV*)SvRV(*av_fetch(e, q, 0));
      double e_iq = SvNV(*av_fetch(elm, i, 0));
      xi = xi + (mix[q] * (SvNV(*av_fetch(a,i,0)) + e_iq ) / (totc + tota));
    }
    sv_setnv(*av_fetch(a,i,TRUE), xi);
  }

  _Norm(newRV_noinc((SV*)a), AlphabetSize);
  return newRV_inc( (SV*)a );
}
