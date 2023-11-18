#!/usr/bin/perl -w
#  Author: Paul Horton
#  Modified by: Kenichiro Imai
#  Organization: Computational Biology Research Center, AIST, Japan
#  Copyright (C) 2009, Paul Horton, All rights reserved.

package AminoAcidPropensity;
use Carp qw(confess);
use strict;
use warnings;



sub name{  $_[0]{name}  }
sub id  {  $_[0]{id}    }


# return propensity of amino acid $residue
sub value{
    ### assert: @_ == 2
    my ($self, $residue) = @_;
    ### assert: length( $residue ) == 1

    exists( $self->{$residue} )  ||  confess( "propensity not defined for \"$residue\"" );

    return $self->{$residue};
}



use overload(

    q{""}=>sub{
	$_[0]->name();
    }

    );


1;

=end

