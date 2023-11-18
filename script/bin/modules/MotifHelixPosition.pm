#  Author:  Yoshinori Fukasawa
#  Organizations:  Department of Computational Biology, University of Tokyo
#  Copyright (C) 2014, Yoshinori Fukasawa, All rights reserved.
#  Creation Date:  2014/1/28

#  This module returns positions of significant MTS motifs and helix for Tom20.
#  This class is a sub-class of ncfHmoment class 

package MotifHelixPosition;
BEGIN{
    use FindBin qw ($Bin);
}

use strict;
use warnings;
use Data::Dumper;
use base qw(ncfHmoment);

sub new{
    my $pkg = shift;

    my %hash;
    my $self = $pkg->SUPER::new();

    #add some member vars
    $self->{motifs} = undef;
    $self->{matchPos} = \%hash;

    bless $self, $pkg;

    my @motifs = (
	'HHPBHH', 'HHBPHH', 'HHHBPH', 'PHHBPH', 'HHBPHB', 'HHBHHB', 'HBHHBb',
	'BHHPPP', 'HHHBBH', 'HPBHHP', 'HBHHbB', 'HHHHBB', 'HHBPHP', 'BPHBHH'
	);

    $self->setMotifs(@motifs);
    return $self;
}

# easy method
sub searchAndGetPositions {
    my ($self, $seq, $length) = @_;

    if( !defined($seq) ){
        return 0;
    }

    if( !defined($length) ){
        $length = 30;
    }

    $self->searchMotifPositions($seq, $length);
    return $self->getAllMotifPositions;
}

sub searchMotifPositions {
    my $self   = shift;
    my $seq    = shift;
    my $length = shift;

    if( !defined($seq) ){
        return 0;
    }

    if( !defined($length) ){
        $length = 30;
    }

    if($seq && $length){
        my $patternedSeq = _changeAminoAcidToPhysicalPattern(substr($seq,0,$length));
	for my $motif (@{$self->{motifs}}){

	    my (@seq_array, @pos_array);
	    while($patternedSeq =~ /(?=($motif))/g){
		push @seq_array, $1;
		push @pos_array, pos($patternedSeq);
	    }

	    my @array;
	    push @array, join("-", ($_+1,$_+length($motif)-1+1)) foreach(@pos_array);
	    $self->{matchPos}->{$motif} = \@array;
	}

    } else {
	print STDERR "Sequence or length are not defined.\n";
	return 0;
    }

}

sub getAllMotifPositions {
    my $self = shift;

    if(keys %{$self->{matchPos}} == 0){
	return 0;
    }
    
    return wantarray ? %{$self->{matchPos}} : $self->{matchPos};
}

sub getMotifPositions {
    my $self = shift;
    my $motif = shift;

    if(!defined $self->{matchPos}->{$motif}){
	return 0;
    }

    return wantarray ? @{$self->{matchPos}->{$motif}} : $self->{matchPos}->{$motif};
}

# Overwrite method default 14 motif patterns.
sub setMotifs {
    my $self   = shift @_;
    my @motifs = @_;

    $self->{motifs}= \@motifs;
}

# Returns current list of motifs
sub getMotifs {
    my $self = shift;
    return wantarray ? @{$self->{motifs}} : $self->{motifs};
}

sub _changeAminoAcidToPhysicalPattern {
    my $seg = shift;
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
