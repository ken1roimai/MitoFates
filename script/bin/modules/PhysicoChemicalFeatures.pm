#!/usr/bin/perl -w
#  Author: Kenichiro Imai
#  Organization: Computational Biology Research Center, AIST, Japan
#  Copyright (C),  Kenichiro Imai, All rights reserved.

#  Purpose: Pick up the Feature of N-terminal Signal peptide
#
package PhysicoChemicalFeatures;
use strict;
use warnings;

sub new{
	my ( $class, $paramHR ) = @_;
	my $self = bless $paramHR, ref($class) || $class;
	return $self;
}

sub divisionNumber{  $_[0]->{divisionNumber};  }
sub regionLength{  $_[0]->{regionLength};      }


sub h_propensity{
    my ( $self, $residue ) = @_;
    $self->{h_propensity}->value($residue);
}
sub pc_propensity{
    my ( $self, $residue ) = @_;
    $self->{pc_propensity}->value($residue);
}
sub nc_propensity{
    my ( $self, $residue ) = @_;
    $self->{nc_propensity}->value($residue);
}
sub am_propensity{
    my ( $self, $residue ) = @_;
    $self->{am_propensity}->value($residue);
}
sub g_propensity{
    my ( $self, $residue ) = @_;
    $self->{g_propensity}->value($residue);
}
sub t_propensity{
    my ( $self, $residue ) = @_;
    $self->{t_propensity}->value($residue);
}
sub s_propensity{
    my ( $self, $residue ) = @_;
    $self->{s_propensity}->value($residue);
}
sub p_propensity{
    my ( $self, $residue ) = @_;
    $self->{p_propensity}->value($residue);
}
sub ar_propensity{
    my ( $self, $residue ) = @_;
    $self->{ar_propensity}->value($residue);
}



sub getPhyscoChemicalFeature{
	my ( $self, $seq ) = @_;
	my $divisionNumber = $self->divisionNumber();
	my $regionLength = $self->regionLength();
	my @seq = split(//,$seq);
	my @hydro = map { $self->h_propensity($_) }  @seq;
	my @pcharge = map { $self->pc_propensity($_) }  @seq;
	my @ncharge = map { $self->nc_propensity($_) }  @seq;
	my @amphi = map { $self->am_propensity($_) }  @seq;
	my @aroma = map { $self->ar_propensity($_) }  @seq;
	my @gly = map { $self->g_propensity($_) }  @seq;
	my @ser = map { $self->s_propensity($_) }  @seq;
	my @thr = map { $self->t_propensity($_) }  @seq;
	my @pro = map { $self->p_propensity($_) }  @seq;
	
	my @hmoment = $self->calcHelicalMoment(\@hydro);
	my @smoment = $self->calcStrandMoment(\@hydro);
	
	#physicochemical features for each segment
	my @hmomentFeatures = $self->calcSegment($regionLength,$divisionNumber,\@hmoment);
	my @smomentFeatures = $self->calcSegment($regionLength,$divisionNumber,\@smoment);
	my @hydroFeatures = $self->calcSegment($regionLength,$divisionNumber,\@hydro);
	my @pchargeFeatures = $self->calcSegment($regionLength,$divisionNumber,\@pcharge);
	my @nchargeFeatures = $self->calcSegment($regionLength,$divisionNumber,\@ncharge);
	my @amphiFeatures = $self->calcSegment($regionLength,$divisionNumber,\@amphi);
	my @aromaFeatures = $self->calcSegment($regionLength,$divisionNumber,\@aroma);
	my @glyFeatures = $self->calcSegment($regionLength,$divisionNumber,\@gly);
	my @serFeatures = $self->calcSegment($regionLength,$divisionNumber,\@ser);
	my @thrFeatures = $self->calcSegment($regionLength,$divisionNumber,\@thr);
	my @proFeatures = $self->calcSegment($regionLength,$divisionNumber,\@pro);
	
	#physicochemical features for Whole sequence
	my $hmomentWhole = $self->calcWhole(\@hmoment);
	my $smomentWhole = $self->calcWhole(\@smoment);
	my $hydroWhole = $self->calcWhole(\@hydro);
	my $pchargeWhole = $self->calcWhole(\@pcharge);
	my $nchargeWhole = $self->calcWhole(\@ncharge);
	my $amphiWhole = $self->calcWhole(\@amphi);
	my $aromaWhole = $self->calcWhole(\@aroma);
	my $glyWhole = $self->calcWhole(\@gly);
	my $serWhole = $self->calcWhole(\@ser);
	my $thrWhole = $self->calcWhole(\@thr);
	my $proWhole = $self->calcWhole(\@pro);
	
	my @results;
	push(@results,@hmomentFeatures);
	push(@results,@smomentFeatures);
	push(@results,@hydroFeatures);
	push(@results,@pchargeFeatures);
	push(@results,@nchargeFeatures);
	push(@results,@amphiFeatures);
	push(@results,@aromaFeatures);
	push(@results,@glyFeatures);
	push(@results,@serFeatures);
	push(@results,@thrFeatures);
	push(@results,@proFeatures);
	
	push(@results,$hmomentWhole);
	push(@results,$smomentWhole);
	push(@results,$hydroWhole);
	push(@results,$pchargeWhole);
	push(@results,$nchargeWhole);
	push(@results,$amphiWhole);
	push(@results,$aromaWhole);
	push(@results,$glyWhole);
	push(@results,$serWhole);
	push(@results,$thrWhole);
	push(@results,$proWhole);
	
	return @results;
}

sub calcSegment{
	my ($self,$segmentTotal,$segmentNumber,$tmpValue) = @_;
	my @value = @$tmpValue;
	my $count = 0;
	if($segmentTotal%$segmentNumber == 0){
		$count = $segmentTotal/$segmentNumber;
	}else{
		print STDERR "Segment error!\n";
	}
	my (@result) = ();
	my $st=0; my $ed=0; my $cnt=0;
	for(my $i=0;$i<$segmentNumber;$i++){ 
		$ed = $st + $count;
		last if($ed > scalar(@value));
		$result[$i] = $self->calcPropertyInRegion($st, $ed, \@value);
		$st = $ed;
		$cnt++;
	}
	return (@result) if($segmentNumber == $cnt);
	for(my $i=0;$i< $segmentNumber-$cnt;$i++){
		push(@result,0);
	}
	return (@result);
}

sub calcWhole{
	my ($self,$tmpValue) = @_;
	my @value = @$tmpValue;
	my $result=0;
	for(my $i=0;$i<=$#value;$i++){
		$result += $value[$i];
	}
	return $result/scalar(@value);
	
}

sub calcPropertyInRegion{
	my ($self,$st,$ed,$tmpvalue) = @_;
	my $cnt = 0;
	my @value = @$tmpvalue;
	my $result = 0.0;
	for(my $i=$st;$i<$ed;$i++){
		$result += $value[$i];
		$cnt++;
	}
	return $result/$cnt;
}


sub calcHelicalMoment{
	my ($self, $hydro) = @_;
	my @wheel;
	my @moment;
	my $end = @$hydro-1;
### treatment of N-C terminal 3 redidues
	for my $a_h (0..2){
		$wheel[$a_h] = ($hydro->[$a_h] + $hydro->[$a_h+3]+ $hydro->[$a_h+4])/3.0;
	}
	{
		my $a2 = 3;
		$wheel[$a2] = ($hydro->[$a2-3]+ $hydro->[$a2] + $hydro->[$a2+3]+ $hydro->[$a2+4])/4.0;
	}
	for my $bb ($end-2..$end){
		$wheel[$bb] = ($hydro->[$bb] + $hydro->[$bb-3]+ $hydro->[$bb-4])/3.0;
	}
	for my $b2 ( $end-3..$end-3){
		$wheel[$b2] = ($hydro->[$b2+3]+$hydro->[$b2] + $hydro->[$b2-3]+ $hydro->[$b2-4])/4.0;
	}
	for my $c (4..$end-4){
		$wheel[$c] = ($hydro->[$c-4]+ $hydro->[$c-3]+ $hydro->[$c] + $hydro->[$c+3]+ $hydro->[$c+4]) / 5.0;
	}
	my $check1;
	my $check2;
	my $check3;
	$moment[0] = 0;
	$moment[$#wheel-1] = 0;
	$moment[$#wheel] = 0;
	my $value1;
	my $value2;
	my $value3;

	for my $d ( 1..$#wheel-2 ){
		$value1 = $wheel[$d-1]-$wheel[$d];
		$value2 = $wheel[$d+1]-$wheel[$d];
		$value3 = $wheel[$d+2]-$wheel[$d];
		$check1 = abs($value1);
		$check2 = abs($value2);
		$check3 = abs($value3);
		if($check1 >= $check2 && $check1 >= $check3){
			$moment[$d] = $check1;
		}elsif($check2 >= $check1 && $check2 >= $check1){
			$moment[$d] = $check2;
		}elsif($check3 >= $check1 && $check3 >= $check2){
			$moment[$d] = $check3;
		}
	}
	return @moment;
}

sub calcStrandMoment{
	my ($self, $hydro) = @_;
	my $end = @$hydro-1;
	my @period_array;
	for(my $a_hs2 = 0; $a_hs2 < 2; $a_hs2++){
		$period_array[$a_hs2] = $hydro->[$a_hs2]*(($hydro->[$a_hs2+2] - $hydro->[$a_hs2+1])/2.0);
	}
	for(my $b_hs2 = $end-1; $b_hs2 < $end+1; $b_hs2++){
		$period_array[$b_hs2] = $hydro->[$b_hs2]*(($hydro->[$b_hs2-2] - $hydro->[$b_hs2-1])/2.0);
	}
	for(my $c_hs2 = 2; $c_hs2 < $end-1; $c_hs2++){
		$period_array[$c_hs2] = $hydro->[$c_hs2]*(($hydro->[$c_hs2+2] + $hydro->[$c_hs2-2] - $hydro->[$c_hs2-1] - $hydro->[$c_hs2+1]) / 4.0);
	}
	return @period_array;
}

1;

