#  Author:  Yoshinori Fukasawa
#  Organizations:  Department of Computational Biology, University of Tokyo
#  Copyright (C) 2011, Yoshinori Fukasawa, All rights reserved.
#  Creation Date:  2011/10/11

package MitoFatesPWM;
require Exporter;
our @ISA = qw (Exporter);
our @EXPORT = qw (readMatrix calculate_pwm);

use strict;
use warnings;
use diagnostics;

sub readMatrix{
    my $path = shift;
    my $name='';
    my $length=0;
    my @array;
    my %hash;
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

# Calculate PWM score
# Argumets:
# $seq         Sequence to be calculated
# $matrix      Hash reference containing matrices
# $identifier  Identifier to pick up matrix from $matrix

sub calculate_pwm{
    if(@_ != 3){
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

    my ($seq, $matrix, $identifier) = @_;
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

1;
