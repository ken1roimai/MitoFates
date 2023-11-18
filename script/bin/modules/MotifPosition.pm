#  Author:  Yoshinori Fukasawa
#  Organizations:  Department of Computational Biology, University of Tokyo
#  Copyright (C) 2014, Yoshinori Fukasawa, All rights reserved.
#  Creation Date:  2014/1/28

package MotifPosition;
use Exporter;
@ISA = (Exporter);
@EXPORT = qw(findMotifPosition);

use warnings;
use strict;

sub findMotifPosition{
    my %args = (
        Seq => "",
        Length => 100,
        Pattern => "spcbp2",
        @_
    );

    chomp($args{Seq});

    my $substr_seq = substr( $args{Seq}, 0, $args{Length} );

    my $pattern = _patternGenerator($args{Pattern});

    my @position_array;
    while($substr_seq =~ /(?=($pattern))/ig){
        my $StartPosition = pos($substr_seq) + 1;
        my $EndPosition = $StartPosition + length($1) - 1;
        push @position_array, "$StartPosition-$EndPosition";
    }

    return @position_array;
}

sub _patternGenerator{

    ### Definition for amino acids in terms of physico-chemical property
    my $_phi    = "FILVWYMCA";
    my $_beta   = "HKR";
    my $_sigma  = "IKPRSV";
    my $_chi    = "ACDEFGHIKLMNPQRSTVWY"; #any letter
    ####################################################################

    my $patternSeq = shift @_;
    my @motif = split(//, $patternSeq);
    my $temporary = "";
    foreach(@motif){
        if($_ =~ /[Ss]/){
            $temporary .= "[$_sigma]";
        } elsif( $_ =~ /[Pp]/ ){
            $temporary .= "[$_phi]";
        } elsif( $_ =~ /[Bb]/ ){
            $temporary .= "[$_beta]";
        } elsif( $_ =~ /[Cc]/ ){
            $temporary .= "[$_chi]";
        } elsif( $_ =~ /[1-9]/ ){
            $temporary .= "{$_}";
        } else {
            die "$_";
        }
    }
    return $temporary;
}
