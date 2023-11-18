#!/usr/bin/perl -w
#  Author: Kenichiro Imai
#  Copyright (C) 2013, Kenichiro Imai, All rights reserved.

package MotifMatch6RScoreTop14;
use strict;
use warnings;

sub new{
	my ( $class, $paramHR ) = @_;
	my $self = bless $paramHR, ref($class) || $class;
	return $self;
}

sub regionLength{  $_[0]->{regionLength};      }

sub getFeature{
	my ( $self, $seq ) = @_;
	my $regionLength = $self->regionLength();
	my @seq = split(//,$seq);
	my $regionEnd = $regionLength-1;
	$regionEnd = $#seq if ($regionEnd > $#seq);
	##**Segmentation of Candidate Region
	my $segment;
	for(0..$regionEnd){
		$segment.=$seq[$_];
	}
	#**Search for Motif in Candidate Region **#
	
	my @motifAScores = $self->searchMTSmotif(0,$regionEnd,\@seq);
	return ($segment,@motifAScores);
}

sub searchMTSmotif{
	my ($self, $strat, $end, $seq) = @_;
	my @seg = @{$seq};
	my @scores;
	for(0..14){
		$scores[$_] = 0;
	}
	for($strat..($end-5)){
		my $pattern = $seg[$_].$seg[$_+1].$seg[$_+2].$seg[$_+3].$seg[$_+4].$seg[$_+5];
		my $physicalPattern = $self->changeAminoAcidToPhysicalPattern($pattern);
		if($physicalPattern =~ /HHPBHH/){
			$scores[0] = 12.24;
			$scores[14]+=12.24;
		}elsif($physicalPattern =~ /HHBPHH/){
			$scores[1] = 10.93;
			$scores[14]+=10.93;
		}elsif($physicalPattern =~ /HHHBPH/){
			$scores[2] = 8.97;
			$scores[14]+=8.97;
		}elsif($physicalPattern =~ /PHHBPH/){
			$scores[3] = 5.91;
			$scores[14]+=5.91;
		}elsif($physicalPattern =~ /HHBPHB/){
			$scores[4] = 8.74;
			$scores[14]+=8.74;
		}elsif($physicalPattern =~ /HHBHHB/){
			$scores[5] = 6.05;
			$scores[14]+=6.05;
		}elsif($physicalPattern =~ /HBHHBb/){
			$scores[6] =8.21;
			$scores[14]+=8.21;
		}elsif($physicalPattern =~ /BHHPPP/){
			$scores[7] = 7.03;
			$scores[14]+= 7.03;
		}elsif($physicalPattern =~ /HHHBBH/){
			$scores[8] = 6.92;
			$scores[14]+= 6.92;
		}elsif($physicalPattern =~ /HPBHHP/){
			$scores[9] = 6.01;
			$scores[14]+= 6.01;
		}elsif($physicalPattern =~ /HBHHbB/){
			$scores[10] = 5.76;
			$scores[14]+= 5.76;
		}elsif($physicalPattern =~ /HHHHBB/){
			$scores[11] = 5.31;
			$scores[14]+= 5.31;
		}elsif($physicalPattern =~ /HHBPHP/){
			$scores[12] = 5.20;
			$scores[14]+= 5.20;
		}elsif($physicalPattern =~ /BPHBHH/){
			$scores[13] = 5.03;
			$scores[14]+= 5.03;
		}
	}
	return @scores;
}

sub changeAminoAcidToPhysicalPattern{
	my ($self, $seg) = @_;
	my @segment = split(//,$seg);
	my $changedPattern;
	for(@segment){
		my $r;
		if($_=~/L|F|I|V|W|Y|M|C|A/){
			$r = "H";
		}elsif($_=~/R|K|H/){
			$r = "B";
		}elsif($_=~/E|D/){
			$r = "A";
		}elsif($_=~/S|T|N|Q/){
			$r = "P";
		}elsif($_=~/P|G/){
			$r = "b";
		}else{
			$r = "X";
		}
		$changedPattern.=$r;
	}
	return $changedPattern;
}

1;

