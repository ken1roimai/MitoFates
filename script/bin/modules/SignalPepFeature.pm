#!/usr/bin/perl -w
#  Author: Kenichiro Imai
#  Organization: Computational Biology Research Center, AIST, Japan
#  Copyright (C),  Kenichiro Imai, All rights reserved.
#
#  Purpose: Pick up the Feature of N-terminal Signal peptide
#
package SignalPepFeature;
use strict;
use warnings;

sub new{
	my ( $class, $paramHR ) = @_;
	my $self = bless $paramHR, ref($class) || $class;
	return $self;
}

sub halfwindowSize{  $_[0]->{halfWindowSize};  }
sub regionLength{  $_[0]->{regionLength};      }
sub hydroThresHold{  $_[0]->{hydroThresHold};  }

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
sub cs_propensity{
    my ( $self, $residue ) = @_;
    $self->{cs_propensity}->value($residue);
}

sub getSpFeature{
	my ( $self, $seq ) = @_;
	
	my $halfWindowSize = $self->halfwindowSize();
	my $regionLength = $self->regionLength();
	my $hydroThresHold = $self->hydroThresHold();
	
	my @seq = split(//,$seq);
	my @hydro = map { $self->h_propensity($_) }  @seq;
	my @pcharge = map { $self->pc_propensity($_) }  @seq;
	my @ncharge = map { $self->nc_propensity($_) }  @seq;
	my @clevageSitePreference = map { $self->cs_propensity($_) }  @seq;
	
	##*****Pick Up Candidate of Signal peptide Region****##
	#######################################################
	
	#**Search for First hydrophobic Peaks**#
	
	#**Calc double average of hydrophobicity with  slidingWindowSize = 21 **#
	my @aveHydro_w21 = $self->averageBySlidingWindow(\@hydro,$halfWindowSize);
	my @doubleAveHydro_w21 = $self->averageBySlidingWindow(\@aveHydro_w21,$halfWindowSize);
	
	#**Search First hydrophobic Peaks **#
	my $hydroMaxPosiAndValue = $self->searchFirstHydrophobicPeak(\@doubleAveHydro_w21, $regionLength, $hydroThresHold);
	my ($hydroMaxPosi,$hydroMaxValue) = split(/\t/,$hydroMaxPosiAndValue);
	
	#**define the 30aa region around firstHydroPhobicPeak **#
	my $signalRegionStart = $hydroMaxPosi-14;
	my $signalRegionEnd = $hydroMaxPosi+15;
	
	#**define the sub region around firstHydroPhobicPeak **#
	my $signalRegionPreionEnd = $signalRegionStart+4;
	my $signalRegionHregionStart = $hydroMaxPosi-9;
	my $signalRegionHregionEnd = $signalRegionStart+19;
	my $signalRegionCregionStart = $signalRegionStart+20;
	my $extendedCregionEnd = $signalRegionEnd+5;
	
	#**define the sub region for HydrophobicPeak Near Nterminal**#
	if($signalRegionStart < 0){
		$signalRegionStart = 0;
		$signalRegionEnd = $hydroMaxPosi+12;
		$signalRegionPreionEnd = $hydroMaxPosi;
		$signalRegionHregionStart = 0;
		$signalRegionHregionEnd = $hydroMaxPosi+7;
		$signalRegionCregionStart = $hydroMaxPosi+8;
	}
	
	#**define the sub region for small protein **#
	$signalRegionEnd = $#seq if($signalRegionEnd > $#seq);
	if($signalRegionCregionStart > $#seq){
		if($#seq > 3){
			$signalRegionCregionStart = $#hydro-4;
		}else{
			$signalRegionCregionStart = 0;
		}
	}
	$signalRegionPreionEnd = $#seq if($signalRegionPreionEnd > $#seq);
	$signalRegionHregionEnd = $#seq if($signalRegionHregionEnd > $#seq);
	$extendedCregionEnd = $#seq if($extendedCregionEnd > $#seq);
	
	#** Pick up Signal peptide Feature **#
	######################################
	
	# N-terminal Positive Charge Density #
	my $positiveChargeDensityNterm = $self->calcSegmentAverage($signalRegionStart,$signalRegionPreionEnd,\@pcharge);
	# N-terminal Negative Charge Density #
	my $negativeChargeDensityNterm = $self->calcSegmentAverage($signalRegionStart,$signalRegionPreionEnd,\@ncharge);
	# Negative Charge Density #
	my $negativeChargeDensity = $self->calcSegmentAverage($signalRegionStart,$signalRegionHregionEnd,\@ncharge);
	# Average of Hydrophobicity in Hydrophobic region #
	my $hydroAverageHregion = $self->calcSegmentAverage($signalRegionHregionStart,$signalRegionHregionEnd,\@hydro);
	# Average of clevage Site Preference redidues: Ala,Gly,Ser after Hydrophobic Peak #
	my $csPreferenceAverageAfterPeak = $self->calcSegmentAverage($hydroMaxPosi,$signalRegionEnd,\@clevageSitePreference);
	
#	# Negative Charge Density in Cleavage region #
#	my $negativeChargeDensityAfterClevageSite = $self->calcSegmentAverage($signalRegionEnd,$extendedCregionEnd,\@ncharge);
	
	return ($positiveChargeDensityNterm,$negativeChargeDensityNterm,$negativeChargeDensity,$hydroAverageHregion,$csPreferenceAverageAfterPeak);
}

sub calcSegmentAverage{
	my ($self, $start, $end, $value) = @_;
	my $sum = 0;
	my $count = 0;
	for($start..$end){
		$sum+= $value->[$_];
		$count++;
	}
	my $average = $sum/$count;
}

sub searchFirstHydrophobicPeak{
	my ($self, $aveHydro, $regionLength, $thresHold) = @_;
	my $seqLeng = @$aveHydro-1;
	$regionLength = $seqLeng if($regionLength > $seqLeng);
	my $hydroClusterFlag = 0;
	my $start = -1;
	my $end = -1;
	my @hydroRegions;
	#*** pick up Candidate region of first peak and hydrophobic value **#
	#*** start is 10 becuse of considering TransMembrane region **#
	for(5..$regionLength){
		if($aveHydro->[$_] >= $thresHold){
			if($hydroClusterFlag == 0){
				$start = $_;
			}
			if($hydroClusterFlag == 1 && $_ == $regionLength){
				$end = $_;
				my ($max,$maxposi) = $self->maxValueAndMaxPosi($start,$end,$aveHydro);
				push(@hydroRegions, $maxposi."\t".$max);
			}
			$hydroClusterFlag = 1 ;
		}elsif($aveHydro->[$_] < $thresHold){
			if($hydroClusterFlag == 1){
				$end = $_;
				my ($max,$maxposi) = $self->maxValueAndMaxPosi($start,$end,$aveHydro);
				push(@hydroRegions, $maxposi."\t".$max);
			}
			$hydroClusterFlag = 0;
			$start = -1;
			$end = -1
		}
	}
	#*** Search hydrophobic first Peak and Value**#
	my $maxPosiAndValue;
	if($#hydroRegions == -1){
		$start = 5;
		my ($max,$maxposi) = $self->maxValueAndMaxPosi($start,$regionLength,$aveHydro);
		$maxPosiAndValue = $maxposi."\t".$max;
	}elsif($#hydroRegions == 0){
		$maxPosiAndValue = $hydroRegions[0];
	}else{
		my ($fistMaxPosi,$firstMax) = split(/\t/,$hydroRegions[0]);
		my %peakCandidate = map { $_->[0] => $_->[1] } map {[ split /\t/ ] } @hydroRegions;
		my @keys = sort {
			$peakCandidate{$b} cmp $peakCandidate{$a}
		}keys %peakCandidate;
		if($keys[0] != $fistMaxPosi){
			if(($keys[0] - $fistMaxPosi) > 5){
				$maxPosiAndValue = $hydroRegions[0];
			}else{
				$maxPosiAndValue = $keys[0]."\t".$peakCandidate{$keys[0]};
			}
		}else{
			$maxPosiAndValue = $hydroRegions[0];
		}
	}
	return $maxPosiAndValue;
}

sub maxValueAndMaxPosi{
	my ($self, $start, $end, $value) = @_;
	my %segment;
	for($start..$end){
		$segment{$_} = $value->[$_];
	}
	my @keys = sort {
		$segment{$b} cmp $segment{$a}
	}keys %segment;
	my $max = $segment{$keys[0]};
	return ($max,$keys[0]);
}

sub averageBySlidingWindow{
	my ($self, $value, $halfwindowSize) = @_;
	my $end = @$value;
	my @average;
	$average[0] = 0;
	$average[$end-1] = 0;
	for(my $k = 1; $k <= $end-2; $k++){
		my $sum = 0;
		my $count = 0;
		my $window = $halfwindowSize;
		$window = $k if($k-$window < 0);
		$window = $end-$k-1 if($k+$window >= $end);
		for(my $n = $k-$window; $n <= $k+$window; $n++){
			$sum+= $value->[$n];
			$count++;
		}
		$average[$k] = $sum/$count;
	}
	return @average;
}

1;

