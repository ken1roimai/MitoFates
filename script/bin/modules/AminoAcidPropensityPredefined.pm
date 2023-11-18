#!/usr/bin/perl -w
#  Author: Paul Horton
#  Modified by: Kenichiro Imai
#  Organization: Computational Biology Research Center, AIST, Japan
#  Copyright (C), Paul Horton and Kenichiro Imai, All rights reserved.

package AminoAcidPropensityPredefined;
use AminoAcidPropensity;
use Carp qw(confess);
use strict;
use warnings;



my %relativeAbundance =
    ( name => 'relative abundance',
      id   => 'NCBI_aminoAcidExplorer_relativeAbundance',
      L	=> 0.0994, A => 0.0884, G => 0.0703, V => 0.0677, S => 0.0672,
      E	=> 0.0624, I => 0.0595, R => 0.0570, T => 0.0543, D => 0.0539,
      K	=> 0.0527, P => 0.0471, N => 0.0417, F => 0.0400, Q => 0.0382,
      Y => 0.0300, M => 0.0237, H => 0.0220, C => 0.0124, W => 0.0121 );

my %EngelmanGES1986 = 
    ( name => 'Engelman GES hydrophobicity',
      id   => 'GES1986',
      A =>  -6.7, R =>  51.5, N => 20.1, D =>  38.5, C =>  -8.4, Q => 17.2,
      E =>  34.3, G =>  -4.2, H => 12.6, I => -13.0, L => -11.7, K => 36.8,
      M => -14.2, F => -15.5, P =>  0.8, S =>  -2.5, T =>  -5.0, W => -7.9,
      Y =>   2.9, V => -10.9, B => 30.5, Z =>  27.8, X =>   5.8 );


my %AcHydro =
    ( name => 'Aboderin Cornette',
      id   => 'AcHydro',
      A =>  0.83, R => -1.74, N => -2.90, D => -2.82, C =>  0.00, Q => -2.24,
      E => -1.91, G => 0.00, H => -2.07, I =>  4.31, L => 4.89, K => -2.32,
      M =>  3.81, F =>  4.56, P => 0.66, S => -0.83, T => -0.50, W => 4.23,
      Y => 3.23, V =>  3.65, X => 0.00 );

my %KYTJ820101 =
    ( name => 'Kyte Doolittle',
      id   => 'KYTJ820101',
      A =>  1.8, R => -4.5, N => -3.5, D => -3.5, C =>  2.5, Q => -3.5,
      E => -3.5, G => -0.4, H => -3.2, I =>  4.5, L =>  3.8, K => -3.9,
      M =>  1.9, F =>  2.8, P => -1.6, S => -0.8, T => -0.7, W => -0.9,
      Y => -1.3, V =>  4.2, B => -3.5, Z => -3.5, X => -0.19 );

my %KLEP840101 =
    ( name => 'Klein net charge',
      id   => 'KLEP840101',
      A =>  0, R =>  1, N => 0,     D => -1,     C => 0, Q => 0,
      E => -1, G =>  0, H => 0,     I =>  0,     L => 0, K => 1,
      M =>  0, F =>  0, P => 0,     S =>  0,     T => 0, W => 0,
      Y =>  0, V =>  0, B => -0.56, Z =>  -0.62, X => -0.01 );

my %Pcharge =
    ( name => 'Positivecharge',
      id   => 'Pcharge',
      A =>  0, R =>  1, N => 0,     D =>  0,     C => 0, Q => 0,
      E =>  0, G =>  0, H => 1,     I =>  0,     L => 0, K => 1,
      M =>  0, F =>  0, P => 0,     S =>  0,     T => 0, W => 0,
      Y =>  0, V =>  0, X => 0 );

my %Ncharge =
    ( name => 'Negativecharge',
      id   => 'Ncharge',
      A =>  0, R =>  0, N => 0,     D => -1,     C => 0, Q => 0,
      E => -1, G =>  0, H => 0,     I =>  0,     L => 0, K => 0,
      M =>  0, F =>  0, P => 0,     S =>  0,     T => 0, W => 0,
      Y =>  0, V =>  0, X => 0 );

my %Ncharge2 =
    ( name => 'Negativecharge2',
      id   => 'Ncharge2',
      A =>  0, R =>  0, N => 0,     D =>  1,     C => 0, Q => 0,
      E =>  1, G =>  0, H => 0,     I =>  0,     L => 0, K => 0,
      M =>  0, F =>  0, P => 0,     S =>  0,     T => 0, W => 0,
      Y =>  0, V =>  0, X => 0 );

my %Tcharge =
    ( name => 'Totalcharge',
      id   => 'Tcharge',
      A =>  0, R =>  1, N => 0,     D => -1,     C => 0, Q => 0,
      E => -1, G =>  0, H => 1,     I =>  0,     L => 0, K => 1,
      M =>  0, F =>  0, P => 0,     S =>  0,     T => 0, W => 0,
      Y =>  0, V =>  0, X => 0 );

my %Pro =
    ( name => 'Prolines',
      id   => 'Proline',
      A =>  0, R =>  0, N => 0,     D =>  0,     C => 0, Q => 0,
      E =>  0, G =>  0, H => 0,     I =>  0,     L => 0, K => 0,
      M =>  0, F =>  0, P => 1,     S =>  0,     T => 0, W => 0,
      Y =>  0, V =>  0, X => 0 );
      
my %Gly =
    ( name => 'Glycines',
      id   => 'Glycine',
      A =>  0, R =>  0, N => 0,     D =>  0,     C => 0, Q => 0,
      E =>  0, G =>  1, H => 0,     I =>  0,     L => 0, K => 0,
      M =>  0, F =>  0, P => 0,     S =>  0,     T => 0, W => 0,
      Y =>  0, V =>  0, X => 0 );

my %Ser =
    ( name => 'Serines',
      id   => 'Serine',
      A =>  0, R =>  0, N => 0,     D =>  0,     C => 0, Q => 0,
      E =>  0, G =>  0, H => 0,     I =>  0,     L => 0, K => 0,
      M =>  0, F =>  0, P => 0,     S =>  1,     T => 0, W => 0,
      Y =>  0, V =>  0, X => 0 );

my %Thr =
    ( name => 'Threonines',
      id   => 'Threonine',
      A =>  0, R =>  0, N => 0,     D =>  0,     C => 0, Q => 0,
      E =>  0, G =>  0, H => 0,     I =>  0,     L => 0, K => 0,
      M =>  0, F =>  0, P => 0,     S =>  0,     T => 1, W => 0,
      Y =>  0, V =>  0, X => 0 );

my %Aroma =
    ( name => 'AromaticResidues',
      id   => 'AromaResi',
      A =>  0, R =>  0, N => 0,     D =>  0,     C => 0, Q => 0,
      E =>  0, G =>  0, H => 0,     I =>  0,     L => 0, K => 0,
      M =>  0, F =>  1, P => 0,     S =>  0,     T => 0, W => 1,
      Y =>  1, V =>  0, X => 0 );
      
my %SMAmphiphilicity =
    ( name => 'Amphiphilicity',
      id   => 'Amphi',
      A =>  0,    R =>  2.45, N => 0,        D =>  0,     C => 0, Q => 1.25,
      E =>  1.27, G =>  0,    H => 1.45,     I =>  0,     L => 0, K => 3.67,
      M =>  0,    F =>  1,    P => 0,        S =>  0,     T => 0, W => 1,
      Y =>  1,    V =>  0,    X => 0 );

my %CSPreference =
    ( name => 'Cleavage site Preference',
      id   => 'GAS',
      A =>  1, R =>  0, N => 0,     D =>  0,     C => 0, Q => 0,
      E =>  0, G =>  1, H => 0,     I =>  0,     L => 0, K => 0,
      M =>  0, F =>  0, P => 0,     S =>  1,     T => 0, W => 0,
      Y =>  0, V =>  0, X => 0 );

my %nameToReference = 
    ( relativeabundance         => \%relativeAbundance,
      engelmangeshydrophobicity => \%EngelmanGES1986,
      kytedoolittle             => \%KYTJ820101,
      AboderinCornette          => \%AcHydro,
      kleinnetcharge            => \%KLEP840101,
      Totalcharge               => \%Tcharge,
      Positivecharge            => \%Pcharge,
      Negativecharge            => \%Ncharge,
      Negativecharge2           => \%Ncharge2,
      CleavageSitePreference    => \%CSPreference,
      Amphiphilicity            => \%SMAmphiphilicity,
      AromaticResidues          => \%Aroma,
      Threonines                => \%Thr,
      Serines                   => \%Ser,
      Glycines                  => \%Gly,
      Prolines                  => \%Pro
    );


my %idToReference = 
    ( NCBI_aminoAcidExplorer_relativeAbundance => \%relativeAbundance,
      GES1986                                  => \%EngelmanGES1986,
      KYTJ820101      => \%KYTJ820101,
      AcHydro         => \%AcHydro,
      KLEP840101      => \%KLEP840101,
      Tcharge         => \%Tcharge,
      Pcharge         => \%Pcharge,
      Ncharge         => \%Ncharge,
      Ncharge2        => \%Ncharge2,
      GAS             => \%CSPreference,
      Amphi           => \%SMAmphiphilicity,
      AromaResi       => \%Aroma,
      Threonine       => \%Thr,
      Serine          => \%Ser,
      Glycine         => \%Gly,
      Proline         => \%Pro
    );



sub get{
    ### assert: @_ == 1
    my $nameOrId = shift;
    ### assert: length( $nameOrId ) > 1

    my $objectRef
	= exists( $nameToReference{$nameOrId} )  ?  $nameToReference{$nameOrId}
    :     exists( $idToReference{$nameOrId}   )  ?  $idToReference{$nameOrId}
    :  confess( "cannot find propensity named \"$nameOrId\"" );

    return bless $objectRef, 'AminoAcidPropensity';
			  
}





1;

=end
