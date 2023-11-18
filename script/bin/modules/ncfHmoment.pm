#  Author:  Yoshinori Fukasawa
#  Organizations:  Department of Computational Biology, University of Tokyo
#  Copyright (C) 2013, Yoshinori Fukasawa, All rights reserved.
#  Creation Date:  2013/10/18
#  This module calculates amphiphilic alpha helix which has positively charged face.

#  Algorithm: calculate mu_N for highest scoring position of mu_H to avoid dubious match.
#  mu_N = mu_H/n - r*cosA*mu_C/n, where n is the length of window, r is the ratio parameter between mu_H and mu_C,
#   and A is the degree between two vectors, namely vectors of hydrophobic moment and charge moment

package ncfHmoment;
BEGIN{
    use FindBin qw ($Bin);
}

use strict;
use warnings;
#use diagnostics;
use Math::Trig;

sub new{
    my $pkg = shift;

    bless {
	windowS => 10,
	windowE => 20,
	moment_x => undef,
	moment_y => undef,
	moment => undef,
	theta => undef,
	maximum_index => undef,
	charge_coord => undef,
	windowSize => undef,
	best_window => undef,
    }, $pkg;

}

# Calculate hydrophobic moment by specified scales
## Arguments:
## seq a sequence to be calclated its hydrophobic moment
## window window size
## length region

sub calcHmoment{

    my $self = shift;

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

    my $Eis = {
        'F' => 0.61, 'M' => 0.26, 'I' => 0.73,
        'L' => 0.53, 'V' => 0.54, 'C' => 0.04,
        'W' => 0.37, 'A' => 0.25, 'T' => -0.18,
        'G' => 0.16, 'S' => -0.26, 'P' => -0.07,
        'Y' => 0.02, 'H' => -0.4, 'Q' => -0.69,
        'N' => -0.64, 'E' => -0.62, 'K' => -1.1,
        'D' => -0.72, 'R' => -1.8,
    };

    my $GES = {
        'F' => 3.7, 'M' => 3.4, 'I' => 3.1,
        'L' => 2.8, 'V' => 2.6, 'C' => 2.0,
        'W' => 1.9, 'A' => 1.6, 'T' => 1.2,
        'G' => 1.0, 'S' => 0.6, 'P' => -0.2,
        'Y' => -0.7, 'H' => -3.0, 'Q' => -4.1,
        'N' => -4.8, 'E' => -8.2, 'K' => -8.8,
        'D' => -9.2, 'R' => -12.3,
    };

    my $C1 = {
        # Aboderin. 1971.Int J Biochem 2:537-544.
        # J. CORNETTE, et al. 1987. J Mol Biol 195:659-685
        'S' => '0', 'F' => '0', 'T' => '0',
        'N' => '0', 'K' => '1', 'Y' => '0',
        'E' => '-1', 'V' => '0', 'Q' => '0',
        'M' => '0', 'C' => '0', 'L' => '0',
        'A' => '0', 'W' => '0', 'P' => '0',
        'H' => '1', 'D' => '-1', 'R' => '1',
        'I' => '0', 'G' => '0',
    };

    my ($seq, $degree, $ratio, $length) = @_;
    my @moment; #2Dmat
    my @moment_deg;
    my $windowS = $self->{windowS};
    my $windowE = $self->{windowE};

    if($seq =~ /[XxUu]/){
	$seq =~ s/[XxUu]//g;
    }

    if( !defined($length) ){
        $length = length($seq);
    }

    if( !defined($ratio) ){
        $ratio = 1;
    }

    if( !defined($degree) ){
	$degree = 100;
    }

    if($seq && $length){
        my @seq = split(//, $seq);

	#
	# First, Calculate normal hydrophobicity moment and scan highest match.
	#

	for my $window($windowS..$windowE){
	    my @temp;
	    my @t; #to contain theta
	    my @d;
	    for(my $i=0; $i<$length-$window+1;$i++){
		my $y = 0;
		my $x = 0;
		my $substr; #ウィンドウ内の部分文字列

		for my $j(1..$window){
		    #die "$seq[$i+$j-1]" if(!defined $ABODR->{$seq[$i+$j-1]});
		    $y += $ABODR->{$seq[$i+$j-1]}*sin(deg2rad(-$degree*$j+90+$degree));
		    $x += $ABODR->{$seq[$i+$j-1]}*cos(deg2rad(-$degree*$j+90+$degree));
		    $substr .= $seq[$i+$j-1];
		}

		my $y2 = $y * $y;
		my $x2 = $x * $x;
		my $mu_H = $x2 + $y2; $mu_H = $mu_H ** 0.5;

		push @temp, $mu_H/$window;
	    }

	    $moment[$window-$windowS] = \@temp;
	}

	($self->{windowSize}, $self->{maximum_index}) = _maxIndex($windowS, \@moment);

	#
        # Then, Calculate mu_N
	#
	my $x   =0;
	my $y   = 0;
	my $x_c = 0;
	my $y_c = 0;
	my $dot = 0;   # dot product of mu_H and mu_C
	my $theta = 0; # degree between given mu_H and mu_C

	for my $j (1..$self->{windowSize}){
	    $y += $ABODR->{$seq[$self->{maximum_index}+$j-1]}*sin(deg2rad(-$degree*$j+90+$degree));
	    $x += $ABODR->{$seq[$self->{maximum_index}+$j-1]}*cos(deg2rad(-$degree*$j+90+$degree));
	    $y_c += $C1->{$seq[$self->{maximum_index}+$j-1]}*sin(deg2rad(-$degree*$j+90+$degree));
	    $x_c += $C1->{$seq[$self->{maximum_index}+$j-1]}*cos(deg2rad(-$degree*$j+90+$degree));
	}

	my $deg = rad2deg(atan2($y, $x));

	$dot = $x*$x_c + $y*$y_c;
	my $mu_C = ($x_c*$x_c + $y_c*$y_c) ** 0.5; #used to calc theta between mu_H and mu_C
	my $mu_H = $moment[$self->{windowSize}-$self->{windowS}]->[$self->{maximum_index}] * $self->{windowSize};

	#print "DEBUG: ", $mu_H, "\n";

	if($mu_H && $mu_C){ #sometimes $mu_[HC] can be zero
	    $theta = $dot/($mu_H * $mu_C); #cos\theta
	    $theta = rad2deg(acos($theta));
	}
	
	my $rotDeg = deg2rad(90-$deg);
	my @coord_mu_C = ($x_c*cos($rotDeg)-$y_c*sin($rotDeg), $x_c*sin($rotDeg)+$y_c*cos($rotDeg) );

	my $new_moment;
	if($mu_H == 0){
	    $new_moment = 0;
	} else {
	    $dot /= $mu_H;
	    $new_moment = $mu_H - ($ratio*$dot);
	}

	$self->{moment} = $new_moment/$self->{windowSize};
	$self->{theta}  = $theta;
	$self->{charge_coord} = join(",", @coord_mu_C );
	$self->{best_window} = substr($seq, $self->{maximum_index}, $self->{windowSize});

	#debug
	#print STDERR $self->{moment_deg}, ", ", @seq[$self->{maximum_index}..$self->{maximum_index}+18], "\n";
    }
    else{
        print STDERR "==\tError: arguments for Hmoment() are not defined.\n";
        return -1;
    }
}

#getters
sub getMoment{
    my $self = shift;

    if(!defined($self->{moment})){
	print STDERR "Do calcHmoment()\n";
	return 0;
    };

    return $self->{moment};
}

sub getPos{
    my $self = shift;

    if(!defined($self->{maximum_index})){
	print STDERR "Do calcHmoment()\n";
	return 0;
    };

    return $self->{maximum_index}+1;
}

sub getWindowSize{
    my $self = shift;

    if(!defined($self->{windowSize})){
	print STDERR "Do calcHmoment()\n";
	return 0;
    };
    return $self->{windowSize};
}

sub getWindow{
    my $self = shift;

    if(!defined($self->{best_window})){
	print STDERR "Do calcHmoment()\n";
	return 0;
    };

    return $self->{best_window};
}

sub getTheta{
    my $self = shift;

    if(!defined($self->{theta})){
	print STDERR "Do calcHmoment()\n";
	return 0;
    };
    return $self->{theta};
}

sub getChargeCoord{
    my $self = shift;

    if(!defined($self->{charge_coord})){
	print STDERR "Do calcHmoment()\n";
	return 0;
    };
    return $self->{charge_coord};
}

sub _maxIndex{
    my ($windowS, $array) = @_;
    my ($max, $maxInd, $window) = (0,0,0);

    for(my $i = 0; $i<@{$array}; $i++){
	for (my $j = 0; $j<@{$array->[$i]}; $j++){
	    if($max < $array->[$i]->[$j]){
		$max = $array->[$i]->[$j];
		$maxInd = $j;
		$window = $i;
	    }
	}
    }

    #since window is just an index, adds $windowS
    return ($window+$windowS, $maxInd);
}

1;
