#!/usr/bin/perl -w
#  Author: Kenichiro Imai and Paul Horton
#  Organization: Computational Biology Research Center, AIST, Japan
#  Copyright (C),  Kenichiro Imai and Paul Horton, All rights reserved.

#  Input:  Amino acid sequence Fasta File
#  Output: SVM feature File
#
#  Bug fix: fixed a bug related hash access by enforced accessing with defined key order (Feb/23/2015)

BEGIN{  
    use FindBin;
    push @INC, "$FindBin::Bin/modules";
}
use Perl6::Slurp;
use Getopt::Long qw(:config posix_default bundling permute );
use Pod::Usage;
use strict;
use warnings;
use AminoAcidPropensityPredefined;
use SignalPepFeature;
use MotifMatch6RScoreTop14;
use PhysicoChemicalFeatures;
use CleavageScorePWM;
use ncfHmoment;
use KeyList;

# options with default values
my $optMan    = 0;
my $optHelp   = 0;
my $optUsage  = 0;
my $optFeatureListFilename  = undef;
my $optOutputFeatureNames   = 0;

# mandatory arguments
my $argClassNumber;

my $getOptionsRetval =
    GetOptions( 
	'man'              => \$optMan,
	'help'             => \$optHelp,
	'usage'            => \$optUsage,
	'c|class-number=i'       =>  \$argClassNumber,
	'f|feature-list-file=s'  =>  \$optFeatureListFilename,
	'output-feature-names'   =>  \$optOutputFeatureNames,
    );

$getOptionsRetval || pod2usage( -verbose => 0 );

if( $optMan ){ 
    if( $ENV{TERM} eq "dumb" ){  pod2usage( -verbose => 2, -output => \*STDERR );  }
    else{                        pod2usage( -verbose => 2 ); }
}
$optHelp    && pod2usage( -verbose => 1 );
$optUsage   && pod2usage( -verbose => 1 );

# ***** Setup some global constants *****
my @FEATURE_NAME  =  sort qw(
A C D E F G H I K L M N P Q R S T V W Y
AA AC AD AE AF AG AH AI AK AL AM AN AP AQ AR AS AT AV AW AY
CA CC CD CE CF CG CH CI CK CL CM CN CP CQ CR CS CT CV CW CY
DA DC DD DE DF DG DH DI DK DL DM DN DP DQ DR DS DT DV DW DY
EA EC ED EE EF EG EH EI EK EL EM EN EP EQ ER ES ET EV EW EY
FA FC FD FE FF FG FH FI FK FL FM FN FP FQ FR FS FT FV FW FY
GA GC GD GE GF GG GH GI GK GL GM GN GP GQ GR GS GT GV GW GY
HA HC HD HE HF HG HH HI HK HL HM HN HP HQ HR HS HT HV HW HY
IA IC ID IE IF IG IH II IK IL IM IN IP IQ IR IS IT IV IW IY
KA KC KD KE KF KG KH KI KK KL KM KN KP KQ KR KS KT KV KW KY
LA LC LD LE LF LG LH LI LK LL LM LN LP LQ LR LS LT LV LW LY
MA MC MD ME MF MG MH MI MK ML MM MN MP MQ MR MS MT MV MW MY
NA NC ND NE NF NG NH NI NK NL NM NN NP NQ NR NS NT NV NW NY
PA PC PD PE PF PG PH PI PK PL PM PN PP PQ PR PS PT PV PW PY
QA QC QD QE QF QG QH QI QK QL QM QN QP QQ QR QS QT QV QW QY
RA RC RD RE RF RG RH RI RK RL RM RN RP RQ RR RS RT RV RW RY
SA SC SD SE SF SG SH SI SK SL SM SN SP SQ SR SS ST SV SW SY
TA TC TD TE TF TG TH TI TK TL TM TN TP TQ TR TS TT TV TW TY
VA VC VD VE VF VG VH VI VK VL VM VN VP VQ VR VS VT VV VW VY
WA WC WD WE WF WG WH WI WK WL WM WN WP WQ WR WS WT WV WW WY
YA YC YD YE YF YG YH YI YK YL YM YN YP YQ YR YS YT YV YW YY
AxxA AxxR AxxN AxxD AxxC AxxQ AxxE AxxG AxxH AxxI AxxL AxxK AxxM AxxF AxxP AxxS AxxT AxxW AxxY AxxV
CxxA CxxR CxxN CxxD CxxC CxxQ CxxE CxxG CxxH CxxI CxxL CxxK CxxM CxxF CxxP CxxS CxxT CxxW CxxY CxxV
DxxA DxxR DxxN DxxD DxxC DxxQ DxxE DxxG DxxH DxxI DxxL DxxK DxxM DxxF DxxP DxxS DxxT DxxW DxxY DxxV
ExxA ExxR ExxN ExxD ExxC ExxQ ExxE ExxG ExxH ExxI ExxL ExxK ExxM ExxF ExxP ExxS ExxT ExxW ExxY ExxV
FxxA FxxR FxxN FxxD FxxC FxxQ FxxE FxxG FxxH FxxI FxxL FxxK FxxM FxxF FxxP FxxS FxxT FxxW FxxY FxxV
GxxA GxxR GxxN GxxD GxxC GxxQ GxxE GxxG GxxH GxxI GxxL GxxK GxxM GxxF GxxP GxxS GxxT GxxW GxxY GxxV
HxxA HxxR HxxN HxxD HxxC HxxQ HxxE HxxG HxxH HxxI HxxL HxxK HxxM HxxF HxxP HxxS HxxT HxxW HxxY HxxV
IxxA IxxR IxxN IxxD IxxC IxxQ IxxE IxxG IxxH IxxI IxxL IxxK IxxM IxxF IxxP IxxS IxxT IxxW IxxY IxxV
KxxA KxxR KxxN KxxD KxxC KxxQ KxxE KxxG KxxH KxxI KxxL KxxK KxxM KxxF KxxP KxxS KxxT KxxW KxxY KxxV
LxxA LxxR LxxN LxxD LxxC LxxQ LxxE LxxG LxxH LxxI LxxL LxxK LxxM LxxF LxxP LxxS LxxT LxxW LxxY LxxV
MxxA MxxR MxxN MxxD MxxC MxxQ MxxE MxxG MxxH MxxI MxxL MxxK MxxM MxxF MxxP MxxS MxxT MxxW MxxY MxxV
NxxA NxxR NxxN NxxD NxxC NxxQ NxxE NxxG NxxH NxxI NxxL NxxK NxxM NxxF NxxP NxxS NxxT NxxW NxxY NxxV
PxxA PxxR PxxN PxxD PxxC PxxQ PxxE PxxG PxxH PxxI PxxL PxxK PxxM PxxF PxxP PxxS PxxT PxxW PxxY PxxV
QxxA QxxR QxxN QxxD QxxC QxxQ QxxE QxxG QxxH QxxI QxxL QxxK QxxM QxxF QxxP QxxS QxxT QxxW QxxY QxxV
RxxA RxxR RxxN RxxD RxxC RxxQ RxxE RxxG RxxH RxxI RxxL RxxK RxxM RxxF RxxP RxxS RxxT RxxW RxxY RxxV
SxxA SxxR SxxN SxxD SxxC SxxQ SxxE SxxG SxxH SxxI SxxL SxxK SxxM SxxF SxxP SxxS SxxT SxxW SxxY SxxV
TxxA TxxR TxxN TxxD TxxC TxxQ TxxE TxxG TxxH TxxI TxxL TxxK TxxM TxxF TxxP TxxS TxxT TxxW TxxY TxxV
VxxA VxxR VxxN VxxD VxxC VxxQ VxxE VxxG VxxH VxxI VxxL VxxK VxxM VxxF VxxP VxxS VxxT VxxW VxxY VxxV
WxxA WxxR WxxN WxxD WxxC WxxQ WxxE WxxG WxxH WxxI WxxL WxxK WxxM WxxF WxxP WxxS WxxT WxxW WxxY WxxV
YxxA YxxR YxxN YxxD YxxC YxxQ YxxE YxxG YxxH YxxI YxxL YxxK YxxM YxxF YxxP YxxS YxxT YxxW YxxY YxxV
NSigPchargeDensity NSigNchargeDensity NSigNchargeDensity2 NSigHydrophobicity NSigCSPreference
HHPBHH HHBPHH HHHBPH PHHBPH HHBPHB HHBHHB HBHHBb BHHPPP HHHBBH HPBHHP HBHHbB HHHHBB HHBPHP BPHBHH TotalMotifScore
hMomentSeg1 hMomentSeg2 hMomentSeg3 hMomentSeg4 hMomentSeg5 hMomentSeg6
sMomentSeg1 sMomentSeg2 sMomentSeg3 sMomentSeg4 sMomentSeg5 sMomentSeg6
hydroSeg1 hydroSeg2 hydroSeg3 hydroSeg4 hydroSeg5 hydroSeg6
pChargeSeg1 pChargeSeg2 pChargeSeg3 pChargeSeg4 pChargeSeg5 pChargeSeg6
nChargeSeg1 nChargeSeg2 nChargeSeg3 nChargeSeg4 nChargeSeg5 nChargeSeg6
amphiSeg1 amphiSeg2 amphiSeg3 amphiSeg4 amphiSeg5 amphiSeg6
aromaSeg1 aromaSeg2 aromaSeg3 aromaSeg4 aromaSeg5 aromaSeg6
glySeg1 glySeg2 glySeg3 glySeg4 glySeg5 glySeg6
serSeg1 serSeg2 serSeg3 serSeg4 serSeg5 serSeg6
thrSeg1 thrSeg2 thrSeg3 thrSeg4 thrSeg5 thrSeg6
proSeg1 proSeg2 proSeg3 proSeg4 proSeg5 proSeg6
hMomentWhole
sMomentWhole
hydroWhole
pchargeWhole
nchargeWhole
amphiWhole
aromaWhole
glyWhole
serWhole
thrWhole
proWhole
cleavageScore1 cleavageScore2
Hmoment
);

my %FEATURE_NAME  =  map {$_ =>1} @FEATURE_NAME;


my $minSeqLength  =  5;
my $aminoAcidAlphabet  =  'ACDEFGHIKLMNPQRSTVWY';
my @mer1Set  =  split( //, $aminoAcidAlphabet );


# ***** Set up global variables *****
my %featureValue;  # modify with setFeatureValue.
    

# ** Read and parse selected features list **
my %selectedFeature;

if( $optFeatureListFilename ){
    my @lines  =  grep m/^[^#]/,  slurp $optFeatureListFilename, {chomp => 1};
    my $listString  =  join( ' ', @lines );
    %selectedFeature = map {$_ => 1}  split( "[ \t]+", $listString );

    for (keys %selectedFeature){
	exists $FEATURE_NAME{$_} || die "unknown feature name \"$_\"";
    }
}else{
    %selectedFeature = %FEATURE_NAME;
}


# ** make sure mandatory args are set with reasonable values
defined( $argClassNumber )
    ||  pod2usage( -verbose => 0, -msg => 'class number not given\n');

( $argClassNumber =~ m{ [-+]? \d+ }x )
    ||  pod2usage( -verbose => 0, -msg => 'class number must be an integer' );



# ****** Set up global variables *****

# initialize set of all possible dimers {AA, AC, ... YY}
my @mer2Set;
for my $r1 (@mer1Set){
    for my $r2 (@mer1Set){  push( @mer2Set, "$r1$r2" );  }
}
my @mer1xx1Set;
for my $r1 (@mer1Set){
    for my $r2 (@mer1Set){  push( @mer1xx1Set, $r1."xx".$r2 );  }
}
my @mer1xxx1Set;
for my $r1 (@mer1Set){
    for my $r2 (@mer1Set){  push( @mer1xxx1Set, $r1."xxx".$r2 );  }
}


# ***** Open input file *****
my $usage = "Usage: $0 fastaFile";

(@ARGV < 2) || die pod2usage( -verbose => 0 );

my $fastaFile;

(@ARGV == 0) || ($ARGV[0] eq '-')
    ?  $fastaFile = \*STDIN 
    :  open( $fastaFile, $ARGV[0] ) || die "could not open file \"$ARGV[0]\"\n";



local $/ = "\n>"; # to read fasta records as lines.

# ***** Loop over fasta file, one record at a time.
while( <$fastaFile> ){
    m{  ^ [#]  }x  &&  next;  # skip comment line at head of file.

	# parse sequence.
	m{ \A [^\n]* \n (.*) \z }xs;  # remove header line
	( my $seq = $1 ) =~ s/[^$aminoAcidAlphabet]//og;
	
	my $seqLength = length( $seq );
	# !! Reset feature values !!
	%featureValue = ();
	
	##***** get N terminal phyicochemical feature
	my $propensity = AminoAcidPropensityPredefined::get( 'AcHydro' );
	my $pc_propensity = AminoAcidPropensityPredefined::get( 'Pcharge' );
	my $nc_propensity = AminoAcidPropensityPredefined::get( 'Ncharge' );
	my $cs_propensity = AminoAcidPropensityPredefined::get( 'GAS' );
	my $getSigPepFeature
    = SignalPepFeature->new( { h_propensity  => $propensity,
				     pc_propensity   => $pc_propensity,
				     nc_propensity   => $nc_propensity,
				     cs_propensity   => $cs_propensity,
				     halfWindowSize   => 7,
				     regionLength   => 70,
				     hydroThresHold => 0.5 }
    );
	my @spFeature = $getSigPepFeature->getSpFeature( $seq );
    setFeatureValue( 'NSigPchargeDensity', $spFeature[0] );
    setFeatureValue( 'NSigNchargeDensity', $spFeature[1] );
    setFeatureValue( 'NSigNchargeDensity2', $spFeature[2] );
    setFeatureValue( 'NSigHydrophobicity', $spFeature[3] );
    setFeatureValue( 'NSigCSPreference', $spFeature[4] );
    
    ##get physicochemical features
	my $am_propensity = AminoAcidPropensityPredefined::get( 'Amphi' );
	my $ar_propensity = AminoAcidPropensityPredefined::get( 'AromaResi' );
	my $t_propensity = AminoAcidPropensityPredefined::get( 'Threonine' );
	my $s_propensity = AminoAcidPropensityPredefined::get( 'Serine' );
	my $g_propensity = AminoAcidPropensityPredefined::get( 'Glycine' );
	my $p_propensity = AminoAcidPropensityPredefined::get( 'Proline' );

	my $getPhysicoChemicalFeature
	= PhysicoChemicalFeatures->new( { 
	                 h_propensity    => $propensity,
				     pc_propensity   => $pc_propensity,
				     nc_propensity   => $nc_propensity,
				     am_propensity   => $am_propensity,
				     ar_propensity   => $ar_propensity,
				     t_propensity    => $t_propensity,
				     s_propensity    => $s_propensity,
				     g_propensity    => $g_propensity,
				     p_propensity    => $p_propensity,
				     divisionNumber  => 6,
				     regionLength    => 90 }
    );
    my @physicoChemicalFeatures = $getPhysicoChemicalFeature->getPhyscoChemicalFeature( $seq );
	setFeatureValue( 'hMomentSeg1', $physicoChemicalFeatures[0] );
	setFeatureValue( 'hMomentSeg2', $physicoChemicalFeatures[1] );
	setFeatureValue( 'hMomentSeg3', $physicoChemicalFeatures[2] );
	setFeatureValue( 'hMomentSeg4', $physicoChemicalFeatures[3] );
	setFeatureValue( 'hMomentSeg5', $physicoChemicalFeatures[4] );
	setFeatureValue( 'hMomentSeg6', $physicoChemicalFeatures[5] );
	setFeatureValue( 'sMomentSeg1', $physicoChemicalFeatures[6] );
	setFeatureValue( 'sMomentSeg2', $physicoChemicalFeatures[7] );
	setFeatureValue( 'sMomentSeg3', $physicoChemicalFeatures[8] );
	setFeatureValue( 'sMomentSeg4', $physicoChemicalFeatures[9] );
	setFeatureValue( 'sMomentSeg5', $physicoChemicalFeatures[10] );
	setFeatureValue( 'sMomentSeg6', $physicoChemicalFeatures[11] );
	setFeatureValue( 'hydroSeg1', $physicoChemicalFeatures[12] );
	setFeatureValue( 'hydroSeg2', $physicoChemicalFeatures[13] );
	setFeatureValue( 'hydroSeg3', $physicoChemicalFeatures[14] );
	setFeatureValue( 'hydroSeg4', $physicoChemicalFeatures[15] );
	setFeatureValue( 'hydroSeg5', $physicoChemicalFeatures[16] );
	setFeatureValue( 'hydroSeg6', $physicoChemicalFeatures[17] );
	setFeatureValue( 'pChargeSeg1', $physicoChemicalFeatures[18] );
	setFeatureValue( 'pChargeSeg2', $physicoChemicalFeatures[19] );
	setFeatureValue( 'pChargeSeg3', $physicoChemicalFeatures[20] );
	setFeatureValue( 'pChargeSeg4', $physicoChemicalFeatures[21] );
	setFeatureValue( 'pChargeSeg5', $physicoChemicalFeatures[22] );
	setFeatureValue( 'pChargeSeg6', $physicoChemicalFeatures[23] );
	setFeatureValue( 'nChargeSeg1', $physicoChemicalFeatures[24] );
	setFeatureValue( 'nChargeSeg2', $physicoChemicalFeatures[25] );
	setFeatureValue( 'nChargeSeg3', $physicoChemicalFeatures[26] );
	setFeatureValue( 'nChargeSeg4', $physicoChemicalFeatures[27] );
	setFeatureValue( 'nChargeSeg5', $physicoChemicalFeatures[28] );
	setFeatureValue( 'nChargeSeg6', $physicoChemicalFeatures[29] );
	setFeatureValue( 'amphiSeg1', $physicoChemicalFeatures[30] );
	setFeatureValue( 'amphiSeg2', $physicoChemicalFeatures[31] );
	setFeatureValue( 'amphiSeg3', $physicoChemicalFeatures[32] );
	setFeatureValue( 'amphiSeg4', $physicoChemicalFeatures[33] );
	setFeatureValue( 'amphiSeg5', $physicoChemicalFeatures[34] );
	setFeatureValue( 'amphiSeg6', $physicoChemicalFeatures[35] );
	setFeatureValue( 'aromaSeg1', $physicoChemicalFeatures[36] );
	setFeatureValue( 'aromaSeg2', $physicoChemicalFeatures[37] );
	setFeatureValue( 'aromaSeg3', $physicoChemicalFeatures[38] );
	setFeatureValue( 'aromaSeg4', $physicoChemicalFeatures[39] );
	setFeatureValue( 'aromaSeg5', $physicoChemicalFeatures[40] );
	setFeatureValue( 'aromaSeg6', $physicoChemicalFeatures[41] );
	setFeatureValue( 'glySeg1', $physicoChemicalFeatures[42] );
	setFeatureValue( 'glySeg2', $physicoChemicalFeatures[43] );
	setFeatureValue( 'glySeg3', $physicoChemicalFeatures[44] );
	setFeatureValue( 'glySeg4', $physicoChemicalFeatures[45] );
	setFeatureValue( 'glySeg5', $physicoChemicalFeatures[46] );
	setFeatureValue( 'glySeg6', $physicoChemicalFeatures[47] );
	setFeatureValue( 'serSeg1', $physicoChemicalFeatures[48] );
	setFeatureValue( 'serSeg2', $physicoChemicalFeatures[49] );
	setFeatureValue( 'serSeg3', $physicoChemicalFeatures[50] );
	setFeatureValue( 'serSeg4', $physicoChemicalFeatures[51] );
	setFeatureValue( 'serSeg5', $physicoChemicalFeatures[52] );
	setFeatureValue( 'serSeg6', $physicoChemicalFeatures[53] );
	setFeatureValue( 'thrSeg1', $physicoChemicalFeatures[54] );
	setFeatureValue( 'thrSeg2', $physicoChemicalFeatures[55] );
	setFeatureValue( 'thrSeg3', $physicoChemicalFeatures[56] );
	setFeatureValue( 'thrSeg4', $physicoChemicalFeatures[57] );
	setFeatureValue( 'thrSeg5', $physicoChemicalFeatures[58] );
	setFeatureValue( 'thrSeg6', $physicoChemicalFeatures[59] );
	setFeatureValue( 'proSeg1', $physicoChemicalFeatures[60] );
	setFeatureValue( 'proSeg2', $physicoChemicalFeatures[61] );
	setFeatureValue( 'proSeg3', $physicoChemicalFeatures[62] );
	setFeatureValue( 'proSeg4', $physicoChemicalFeatures[63] );
	setFeatureValue( 'proSeg5', $physicoChemicalFeatures[64] );
	setFeatureValue( 'proSeg6', $physicoChemicalFeatures[65] );
	setFeatureValue( 'hMomentWhole', $physicoChemicalFeatures[66] );
	setFeatureValue( 'sMomentWhole', $physicoChemicalFeatures[67] );
	setFeatureValue( 'hydroWhole', $physicoChemicalFeatures[68] );
	setFeatureValue( 'pchargeWhole', $physicoChemicalFeatures[69] );
	setFeatureValue( 'nchargeWhole', $physicoChemicalFeatures[70] );
	setFeatureValue( 'amphiWhole', $physicoChemicalFeatures[71] );
	setFeatureValue( 'aromaWhole', $physicoChemicalFeatures[72] );
	setFeatureValue( 'glyWhole', $physicoChemicalFeatures[73] );
	setFeatureValue( 'serWhole', $physicoChemicalFeatures[74] );
	setFeatureValue( 'thrWhole', $physicoChemicalFeatures[75] );
	setFeatureValue( 'proWhole', $physicoChemicalFeatures[76] );
	
	#get cleavage Score
	my $getCleavageScore
	= CleavageScorePWM->new( {
	    regionLength => 100,
	    pwmDatabasePath => "$FindBin::Bin/pwm/MPP.pwm",
	    fungiFlag => 1,
				 }
	);
	my @cleavageScores = $getCleavageScore->calcCleavageScore($seq);
	my @sortCleavageScores = sort { $b <=> $a } @cleavageScores;
	setFeatureValue( 'cleavageScore1', $sortCleavageScores[0] );
	setFeatureValue( 'cleavageScore2', $sortCleavageScores[1] );
	
	##get Motif score
	my $getMotifFeature
	= MotifMatch6RScoreTop14->new( {
					     regionLength   => 30,
					        }
	);
	my @candidateRegionFeature = $getMotifFeature->getFeature( $seq );
#	my $candidateSegment = $candidateRegionFeature[0];
	setFeatureValue( 'HHPBHH', $candidateRegionFeature[1] );
	setFeatureValue( 'HHBPHH', $candidateRegionFeature[2] );
	setFeatureValue( 'HHHBPH', $candidateRegionFeature[3] );
	setFeatureValue( 'PHHBPH', $candidateRegionFeature[4] );
	setFeatureValue( 'HHBPHB', $candidateRegionFeature[5] );
	setFeatureValue( 'HHBHHB', $candidateRegionFeature[6] );
	setFeatureValue( 'HBHHBb', $candidateRegionFeature[7] );
	setFeatureValue( 'BHHPPP', $candidateRegionFeature[8] );
	setFeatureValue( 'HHHBBH', $candidateRegionFeature[9] );
	setFeatureValue( 'HPBHHP', $candidateRegionFeature[10] );
	setFeatureValue( 'HBHHbB', $candidateRegionFeature[11] );
	setFeatureValue( 'HHHHBB', $candidateRegionFeature[12] );
	setFeatureValue( 'HHBPHP', $candidateRegionFeature[13] );
	setFeatureValue( 'BPHBHH', $candidateRegionFeature[14] );
	setFeatureValue( 'TotalMotifScore', $candidateRegionFeature[15] );
	
	my $compositionSearchLength = 30;
	my $candidateSegment = substr($seq,0,$compositionSearchLength);
	
	##get Hmoment
	my $hmom = ncfHmoment->new();
	$hmom->calcHmoment($candidateSegment, 96, 8.5);
	my $HmomentScore = $hmom->getMoment;
	setFeatureValue( 'Hmoment', $HmomentScore);
	
	my $segLength = length($candidateSegment);
	##aa composition
	(  $segLength > $minSeqLength  )
	||  die( "expected seq length at least $minSeqLength, but read seq: \"$candidateSegment\"" );
	
	# compute list of 1-mers
	my @mer1 = split( //, $candidateSegment);
	# ****** Compute list of 2-mers ******
	my @mer2 = ();
	{
	#  "abcd" --> (a,b,c) + (b,c,d) --> (ab,bc,cd)
	my @mer1CopyForResidue1 = @mer1;  pop   @mer1CopyForResidue1;
	my @mer1CopyForResidue2 = @mer1;  shift @mer1CopyForResidue2;
		@mer2 = map { my $r1 = $_;
			my $r2 = shift(@mer1CopyForResidue2);  # LOOP PROGRESS (2)
			"$r1$r2"
		} @mer1CopyForResidue1;                              # LOOP PROGRESS (1)
	}
	
	# ****** Compute list of AxxA-mers ******
	my @mer1xx1 = ();
	{
	#  "abcdef" --> (a,b,c,d) + (d,e,f,g) --> (axxd,bxxe,cxxf)
		my @mer1CopyForResidue1x = @mer1;  pop   @mer1CopyForResidue1x;  pop   @mer1CopyForResidue1x; pop   @mer1CopyForResidue1x;
		my @mer1CopyForResidue4x = @mer1;  shift @mer1CopyForResidue4x;  shift @mer1CopyForResidue4x; shift @mer1CopyForResidue4x;
		@mer1xx1 = map { my $x1 = $_;
			my $x4 = shift(@mer1CopyForResidue4x);  # LOOP PROGRESS (2)
			$x1."xx".$x4
		} @mer1CopyForResidue1x;                              # LOOP PROGRESS (1)
	}
	
	# ****** Compute 1-mer, 2-mer and 1x1-mers frequencies ******
	# tally 1-mer, 2-mer AxxA, AxxxA counts
	my %mer1Freq = ();  ++$mer1Freq{$_}  for( @mer1 );
	my %mer2Freq = ();  ++$mer2Freq{$_}  for( @mer2 );
	my %mer1xx1Freq = ();  ++$mer1xx1Freq{$_}  for( @mer1xx1 );
	
    # normalize
	for my $r (@mer1Set){
	if( exists($mer1Freq{$r}) ) {  $mer1Freq{$r} /=  $segLength;     }
	else                        {  $mer1Freq{$r} =  0;               }
	}
	for my $r (@mer2Set){
	if( exists($mer2Freq{$r}) ) {  $mer2Freq{$r} /= ($segLength-1);  }
	else                        {  $mer2Freq{$r} =  0;               }
	}
	for my $r (@mer1xx1Set){
	if( exists($mer1xx1Freq{$r}) ) {  $mer1xx1Freq{$r} /= ($segLength-3);  }
	else                          {  $mer1xx1Freq{$r} =  0;               }
	}
	# ***** copy %mer1Freq, %mer2Freq values into %featureValue
	{
	my ( $s, $freq );
		setFeatureValue( $s, 100*$freq )  while(  ($s, $freq) = each %mer1Freq  );
		setFeatureValue( $s, 100*$freq )  while(  ($s, $freq) = each %mer2Freq  );
		setFeatureValue( $s, 100*$freq )  while(  ($s, $freq) = each %mer1xx1Freq );
	}
	
	# ***** Print Output *****
	printf "%+2d ", $argClassNumber;
	
	my $featureCount = 0;
	#for my $feature (keys %selectedFeature){
	foreach my $feature (@{$KeyList::key_list}){
		if( $optOutputFeatureNames ){
			printf  "%s:%-5.3f ",  $feature,        $featureValue{$feature};
		}else{
			printf  "%d:%-5.3f ",  ++$featureCount, $featureValue{$feature};
		}
	}
	
	print "\n";
} # go read next record.



# Open featureListFilename and return feature name list
sub readFeatureListFile{
    (@_ == 1) || die;
    my $featureListFilename = shift;

    my $featureListFile;

    if( $optFeatureListFilename ){
	open( $featureListFile, $optFeatureListFilename )
	    || pod2usage( -verbose => 0,
			  -msg => "could not open feature list file \"$optFeatureListFilename\"" );
    }

}

sub setFeatureValue{
    (@_ == 2) || die;
    my ( $feature, $value ) = @_;
    exists $FEATURE_NAME{$feature} || die( "unknown feature \"$feature\"" );
    $featureValue{$feature} = $value;
}




=pod

=head1 NAME
computeSVMFeatureForPresequence.pl

=head1 SYNOPSIS

B<computeSVMFeatureForPresequence.pl> [B<-f> I<featureListFile>] B<-c> I<classNumber> [I<fastaInputFile>]

B<computeSVMFeatureForPresequence.pl> (B<--usage>|B<--help>|B<--man>)

=head1 OPTIONS

=over 8

=item B<--output-feature-names>

Output feature names with values, instead of the integer feature indices.


=item B<-f>,<--feature-list-file> I<feature-list-file>

Load selected features from I<features-list-file>.  The file format
is feature names separated by white space or new lines. Spaces are
not allowed in feature names. Lines starting with '#' mark comments.

=back

=head1 DEPENDENCIES

Depends on modules Perl6::Slurp.

=cut
