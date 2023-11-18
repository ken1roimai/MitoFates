#!/usr/bin/perl
#  Author: K Imai
#  Modified by: Yoshinori Fukasawa
#  Modification Date:  2012/7/21
#  Copyright (C) 2013, Kenichiro Imai and Yoshinori Fukasawa, All rights reserved.

package CleavageScorePWM;
use strict;
use warnings;
use Carp;
use Math::Cephes qw(gamma gdtr lgam :constants);

sub new{
    my $class = shift @_;
    my $self = shift @_;
    croak "Specify the path.\n" unless(defined($self->{pwmDatabasePath}));
    $self->{fungiFlag} = 1 unless(defined($self->{fungiFlag}));
    return bless $self, $class;
}

sub calcCleavageScore{
    my ( $self, $seq ) = @_;
    #Threshold of lenght for cleavage Site searching
    my $lengthThreshold = $self->{regionLength};
    #pwm PATH
    my $profile  = $self->{pwmDatabasePath};
    my $matrices = $self->_readMatrix($profile);
    my $flag     = $self->{fungiFlag};
    my @segments = $self->_cut10residuesSegment( $seq, $lengthThreshold );
    my $GammaParams = $self->_rtnGammaParams($flag);
    #Calc cleavage scores for each segments
    my @segmentScores;
    for(0..$#segments){
        my $score = $self->_calculate_pwm( $segments[$_], $matrices, 'MPP' );
        ##Weight Score By Gamma Distribution
        my $weightedScore = $self->_GammaWeighting( $_+3, \$score, $GammaParams );
        push( @segmentScores, $weightedScore);
    }

=forComment
    ##PickUp Best and 2nd Score
    my @sortScores = sort { $b <=> $a } @segmentScores;
    my $bestScore = $sortScores[0];
    my $SecondScore = $sortScores[1];
    return ($bestScore,$SecondScore);
=cut

    return @segmentScores;
}

sub _cut10residuesSegment{
    my ($self, $seq, $lengthThreshold) = @_;
    my @seq = split(//,$seq);
    my $end = $lengthThreshold-9;
    $end = $#seq -9 if($#seq < $lengthThreshold);
    my @segments;
    for(0..$end){
        my $segment = substr($seq, $_, 10);
        push(@segments, $segment);
    }
    return @segments;
}

sub _calculate_pwm{
    if(@_ != 4){
        print STDERR "==\tcalculate_pwm\n";
        print STDERR "==\tERROR: The number of arguments is incorrect.\n";
        exit(-1);
    }

    my %amino2index = (
        'A'=>0, 'C'=>1, 'D'=>2,
        'E'=>3, 'F'=>4, 'G'=>5,
        'H'=>6, 'I'=>7, 'K'=>8,
        'L'=>9, 'M'=>10,'N'=>11,
        'P'=>12,'Q'=>13,'R'=>14,
        'S'=>15,'T'=>16,'V'=>17,
        'W'=>18,'Y'=>19,
    );

    my ($self, $seq, $matrix, $identifier) = @_;
    my @seq = split(//, $seq);
    my $score=0;
    my @mat;
    my %checker;

    foreach my $key (keys %$matrix){
        if($key =~ /$identifier/ && !defined($checker{$identifier})){
            @mat = @{$matrix->{$key}};
            $checker{$identifier} = 1;
        } elsif(defined($checker{$identifier}) && $key =~ /$identifier/){
            print STDERR "Bad identifier $identifier was used.\n";
            print STDERR "There are more than one matrices which matched $identifier.\n";
            exit(-1);
        } else {
            next;
        }
    }

    die "There is no matrix identified by $identifier." if(!defined($checker{$identifier}));
    #die "Sequence length does not match with that of matrix. ",$#mat if(length($seq) != $#mat+1);

    for(my $i = 0; $i < @seq; $i++){
        $score += $mat[$i][$amino2index{$seq[$i]}] if($i != 7 && $i != 9); #!= 7 is masking test
    }
    return $score;
}

sub getPositionFromHighest{

    my ($self, $scoreArrayRef, $length) = @_;
    my @segmentScores = @{$scoreArrayRef};
    my %rank2pos;

    my @Rankings = ();
    foreach my $score (@segmentScores){
        my $iRes = grep {$_ > $score} @segmentScores;
        push(@Rankings, $iRes+1);
    }

    for(my $i=0;$i<@Rankings;$i++){
        $rank2pos{$Rankings[$i]} = $i;
    }

    my @tmp;
    for my $rank (sort {$a<=>$b} keys %rank2pos){
        last if($rank > $length);
        push (@tmp, $rank2pos{$rank});
    }

    return \@tmp;
}

sub _readMatrix{
    my $self = shift;
    my $path = shift;
    my $name='';
    my $length=0;
    my @array;
    my %hash;
    local $/ = "\n"; # change indent definition.

    open(my $MAT, "<", $path) || die "Can't open PWM file.\n $!\n";
    while(<$MAT>){
        chomp($_);
        if(/^MitoFates/){
            $name = '';
            $length = 0;
            @array = ();
        } elsif(/^NAME/){
            $_ =~ s/^NAME\t//;
            $name = $_;
        } elsif(/^LENGTH/){
            $_ =~ s/LENGTH\t//;
            $length = $_;
        } elsif(/^PWM/){
            for(my $i=0; $i<$length;$i++){
                my $line = <$MAT>;
                chomp($line);
                my @tmp = split(/\t/, $line);
                shift @tmp;
                $array[$i] = \@tmp;
            }
        } elsif(/^\/Matrix/){
            my @copyarray = @array; # To avoid overwriting of @array.
            $hash{$name} = \@copyarray;

        } elsif(/^#/ || /^PRIOR/ || /^DATE/ || $_ eq ''){
            next; # comment, supp info and empty lines are ignored.
        } else {
            print STDERR "==\tERROR: PWM file includes undefined line. $_\n";
            exit(-1);
        }
    }
    close($MAT);
    return \%hash;
}

## Calculate Gamma Distribution's CDF
sub _gammaCDF{
    my ($slef, $x, $shape, $scale) = @_;
    my $beta = 1/$scale; #a.k.a rate parameter
    my $y = gdtr($x, $shape, $beta);
    return $y;
}

## To weight logarithmic ratio of hmmer by gamma mixture's pdf
## Arguments:
## position same as the above
## scoreRef reference to hmmer score

## Define P($position) = \int_{$position-1}^{$position} f_{X}(x)dx
## => F_{X}($position) - F_{X}($position-1)
sub _GammaWeighting{
    my ($self, $position, $scoreRef, $paramsRef) = @_;
    my $score = $$scoreRef;
    my @params = @{$paramsRef};

    # magic number for minimum
    # if an error occurs, change this to an arbitarary small number.
    my $MINIMUM = 1e-308;

    # Calculate probability
    my $y=0;
    for(my $i=0; $i<@params; $i++){
        my $x = $position;
        my $a = $params[$i]->{"shape"};
        my $b = $params[$i]->{"scale"};
        my $r = 1/$params[$i]->{"scale"}; #rate

        #checkpoint for gdtr = 1-igamc(complemented incomplete gamma integral)
        if($a*$x > 1 && $a*$x > $r){  
            my $ax = $r*log($a*$x)-($a*$x)-lgam($r);
            # Check for underflow. MAXLOG is a constant of Math::Cephes
            if($ax >= -1*$MAXLOG){
                $y += $params[$i]->{"lambda"}*(
                    $self->_gammaCDF($x, $a, $b) - $self->_gammaCDF($x-1, $a, $b)
                );
            }
             else {
                 $y += $MINIMUM; #gammacdf(x)-$gammacdf(x-1) = 1-1
             }
        }
         else{
             $y += $params[$i]->{"lambda"}*(
                 $self->_gammaCDF($x, $a, $b) - $self->_gammaCDF($x-1, $a, $b)
             );
         }
    }
    $$scoreRef += $self->_Log2($y);
}

sub _rtnGammaParams{
    die "The number of arguments is incorrect." if(@_ != 2);
    my $self = shift;
    my $flag = shift;
    my $params;

    if($flag){

        $params = [
            {lambda => 1,
             shape  => 3.368054,
             scale  => 7.702288,
         },
        ];

    }
    else{

        # from Huang et al., 2009,Plant physiology
        $params = [
            {lambda => 0.6262201,
             shape  => 23.110278,
             scale  => 1.1410110,},

            {lambda => 0.2433623,
             shape  => 9.7935300,
             scale  => 4.6511930,},

            {lambda => 0.1304175,
             shape  => 65.305173,
             scale  => 1.4530050,},
        ];

    }

    return $params;
}


sub _Log2{
  my ($self, $x) = @_;
  return log($x)/log(2);
}

1;
