#  Author:  Yoshinori Fukasawa
#  Organizations:  Department of Computational Biology, University of Tokyo
#  Copyright (C) 2015, Yoshinori Fukasawa, All rights reserved.
#  Creation Date:  2015/2/22
#  This module contains an array variable to output features correctly. 
#   Its dimension is very long, so set aside this var to keep the source clean.

package KeyList;

{
    our $key_list = [
	HxxV,
	KxxD,
	AT,
	CxxY,
	WM,
	YT,
	HHBPHH,
	GxxK,
	KK,
	RxxC,
	FH,
	HHBPHB,
	VW,
	SH,
	IxxK,
	hydroSeg6,
	TxxQ,
	YxxH,
	MP,
	WD,
	nChargeSeg2,
	IxxG,
	QxxG,
	QxxS,
	amphiSeg3,
	WxxE,
	EN,
	PW,
	RQ,
	HxxT,
	FxxW,
	FC,
	TK,
	HxxM,
	EE,
	FxxP,
	MxxE,
	YF,
	CK,
	sMomentSeg5,
	AA,
	EK,
	HHBPHP,
	M,
	AxxF,
	EP,
	LG,
	WxxN,
	NSigPchargeDensity,
	AH,
	VxxF,
	MG,
	SR,
	KxxS,
	CxxC,
	FT,
	MxxM,
	WxxY,
	SxxT,
	YY,
	YL,
	RY,
	thrSeg1,
	LxxQ,
	CC,
	RR,
	NC,
	IxxR,
	proSeg5,
	ExxL,
	FxxM,
	HQ,
	LL,
	FQ,
	LxxA,
	MxxW,
	sMomentWhole,
	AF,
	VxxT,
	KE,
	TxxF,
	CxxH,
	WY,
	HD,
	F,
	ID,
	IxxY,
	FxxE,
	QT,
	GM,
	AxxQ,
	GK,
	HL,
	thrSeg6,
	HxxE,
	RxxK,
	FI,
	ExxT,
	SxxN,
	VE,
	IL,
	R,
	RC,
	NK,
	VxxL,
	YxxI,
	VxxG,
	TS,
	CxxP,
	TG,
	RxxV,
	KxxR,
	KC,
	FA,
	nChargeSeg3,
	HV,
	WxxF,
	PA,
	amphiSeg2,
	CxxD,
	QA,
	CV,
	GD,
	amphiSeg6,
	CxxI,
	TxxL,
	MF,
	GxxE,
	HPBHHP,
	NSigCSPreference,
	HBHHBb,
	IxxF,
	QxxF,
	QxxP,
	PxxE,
	KxxC,
	QxxA,
	IR,
	RxxD,
	PxxL,
	SxxQ,
	KxxK,
	WxxM,
	sMomentSeg6,
	AxxW,
	MK,
	EQ,
	TxxR,
	AG,
	HG,
	KT,
	aromaSeg5,
	hMomentSeg2,
	AI,
	pchargeWhole,
	AxxE,
	YM,
	AN,
	VV,
	WL,
	ExxQ,
	NxxR,
	TxxE,
	TxxY,
	HY,
	EH,
	FxxR,
	IG,
	IxxS,
	KN,
	TxxD,
	aromaSeg4,
	SC,
	DxxM,
	IxxL,
	T,
	GxxN,
	MS,
	NL,
	EC,
	ExxF,
	MY,
	NxxM,
	sMomentSeg3,
	FxxH,
	ExxK,
	KxxY,
	cleavageScore2,
	EV,
	DxxA,
	RS,
	ExxS,
	YxxY,
	CN,
	thrWhole,
	SxxI,
	serSeg1,
	LxxD,
	PxxV,
	NxxP,
	QxxI,
	AxxP,
	KD,
	LxxT,
	G,
	proSeg2,
	FxxN,
	FS,
	HxxW,
	AxxH,
	LxxI,
	DxxW,
	VxxY,
	ND,
	TP,
	VxxK,
	DT,
	GT,
	DF,
	QI,
	DxxY,
	HxxN,
	HC,
	HF,
	WG,
	PxxS,
	DW,
	PI,
	SK,
	SY,
	MxxF,
	VxxS,
	S,
	DxxV,
	AxxI,
	GL,
	RxxE,
	KI,
	CD,
	CxxT,
	VC,
	WxxW,
	NM,
	LC,
	TxxW,
	RxxW,
	MN,
	EY,
	PxxC,
	YxxK,
	WS,
	VxxA,
	D,
	NV,
	FxxF,
	IxxE,
	WxxV,
	PxxD,
	EL,
	GxxC,
	nchargeWhole,
	KL,
	YR,
	AxxV,
	FV,
	CT,
	C,
	LQ,
	FE,
	GxxI,
	CxxA,
	GxxM,
	NxxQ,
	CM,
	nChargeSeg1,
	HS,
	HxxG,
	serSeg4,
	GA,
	ExxN,
	TM,
	hydroSeg2,
	YxxR,
	IY,
	IxxC,
	LA,
	PH,
	SQ,
	YxxF,
	TE,
	TF,
	LxxS,
	QxxD,
	TQ,
	thrSeg4,
	EW,
	KxxH,
	VxxR,
	FxxQ,
	NxxY,
	LxxC,
	DM,
	NE,
	DxxH,
	AD,
	RA,
	WT,
	LxxN,
	HxxC,
	GS,
	WN,
	ExxE,
	ME,
	CxxW,
	WxxR,
	VN,
	SxxG,
	YE,
	IP,
	HHHHBB,
	DxxL,
	MxxC,
	KxxM,
	MxxP,
	PxxW,
	aromaSeg6,
	DV,
	FK,
	DE,
	MV,
	WxxD,
	NxxN,
	MI,
	GI,
	WxxL,
	CS,
	HxxQ,
	RT,
	YK,
	WxxP,
	TxxH,
	PS,
	hMomentSeg4,
	sMomentSeg1,
	WE,
	YxxL,
	K,
	FxxV,
	E,
	Y,
	RxxI,
	VP,
	SxxL,
	CE,
	KM,
	QH,
	pChargeSeg4,
	PC,
	GxxL,
	MQ,
	serSeg3,
	SW,
	CxxK,
	CxxN,
	VY,
	LF,
	WA,
	VxxQ,
	glyWhole,
	EM,
	AL,
	DxxN,
	NxxC,
	HW,
	VD,
	GxxH,
	IxxD,
	MD,
	GN,
	L,
	NSigNchargeDensity2,
	amphiSeg4,
	nChargeSeg6,
	HxxA,
	QQ,
	ExxM,
	AP,
	NW,
	DxxK,
	QY,
	pChargeSeg3,
	AQ,
	YxxT,
	sMomentSeg4,
	RxxR,
	YxxQ,
	MxxV,
	CL,
	HxxF,
	GxxT,
	ET,
	IS,
	DL,
	amphiSeg5,
	ExxH,
	NxxK,
	HI,
	TI,
	IF,
	PxxY,
	SxxA,
	GxxS,
	Hmoment,
	IxxM,
	YD,
	QG,
	TxxG,
	MH,
	DxxG,
	KxxN,
	HHBHHB,
	HA,
	TL,
	FxxT,
	DY,
	ExxA,
	LxxG,
	PxxG,
	SP,
	WxxQ,
	YS,
	AxxG,
	PxxM,
	KxxT,
	WxxK,
	IM,
	MxxD,
	RI,
	hMomentSeg5,
	SxxF,
	IxxT,
	DR,
	RxxQ,
	LN,
	FL,
	KxxI,
	MW,
	NN,
	glySeg4,
	GY,
	IA,
	PK,
	FxxC,
	QxxR,
	WK,
	NxxT,
	HxxY,
	NF,
	TD,
	CxxV,
	HP,
	TxxP,
	GV,
	MxxT,
	YxxS,
	AxxK,
	CQ,
	ExxW,
	HHHBBH,
	YxxD,
	AY,
	TxxI,
	PR,
	proSeg1,
	HxxP,
	WV,
	ExxP,
	DI,
	QxxW,
	DP,
	sMomentSeg2,
	MA,
	pChargeSeg2,
	SV,
	MT,
	aromaSeg3,
	LxxV,
	hMomentSeg6,
	PxxI,
	thrSeg2,
	RxxS,
	MxxY,
	hMomentSeg1,
	VxxP,
	LR,
	TxxT,
	TN,
	MxxQ,
	PxxN,
	PxxF,
	SxxK,
	Q,
	IC,
	KR,
	WxxT,
	TotalMotifScore,
	YH,
	GC,
	AM,
	CxxR,
	KxxA,
	NT,
	RH,
	DxxT,
	QxxK,
	LxxM,
	GxxR,
	DD,
	YxxM,
	ExxG,
	WxxI,
	CF,
	LE,
	AE,
	proWhole,
	TH,
	FxxI,
	VH,
	NSigNchargeDensity,
	EA,
	amphiWhole,
	RxxA,
	aromaSeg1,
	hydroSeg1,
	IxxA,
	NxxL,
	KF,
	SxxC,
	VL,
	WP,
	HxxI,
	EF,
	WQ,
	QF,
	NxxD,
	SE,
	IV,
	NQ,
	FxxS,
	NG,
	QP,
	YC,
	H,
	LK,
	RV,
	VxxI,
	WH,
	KxxP,
	LxxH,
	FM,
	DxxF,
	AxxA,
	LV,
	HHPBHH,
	IH,
	KQ,
	FxxD,
	NxxS,
	LH,
	YP,
	VxxC,
	DS,
	TV,
	VI,
	QN,
	RxxH,
	SN,
	HH,
	GG,
	SxxS,
	YxxW,
	WR,
	NxxV,
	KH,
	AxxY,
	MxxS,
	RxxP,
	HxxR,
	QV,
	DK,
	LxxP,
	DC,
	QxxT,
	TA,
	YxxV,
	ExxV,
	FxxA,
	NH,
	LxxW,
	TxxV,
	FW,
	DQ,
	LxxE,
	LP,
	HR,
	GxxG,
	DxxE,
	KxxV,
	GxxQ,
	PxxH,
	VQ,
	hMomentWhole,
	SxxV,
	BPHBHH,
	SxxM,
	WxxA,
	WxxH,
	ST,
	KS,
	GxxY,
	hydroSeg3,
	TY,
	IxxN,
	HxxS,
	EI,
	CG,
	VR,
	GW,
	WC,
	YxxE,
	GR,
	VxxV,
	FxxL,
	P,
	FG,
	DxxS,
	RxxT,
	hydroWhole,
	PP,
	AxxL,
	TxxA,
	RG,
	serSeg5,
	AxxT,
	SxxD,
	HxxH,
	TxxN,
	HK,
	VM,
	CR,
	HxxK,
	VA,
	QS,
	LD,
	glySeg3,
	SD,
	NxxF,
	CxxS,
	VxxW,
	thrSeg5,
	PY,
	WxxC,
	YQ,
	WxxS,
	LI,
	CH,
	KG,
	GH,
	DxxI,
	QE,
	NxxA,
	VS,
	YV,
	I,
	PE,
	RxxY,
	PxxP,
	SL,
	NSigHydrophobicity,
	KxxG,
	SxxE,
	GxxF,
	glySeg1,
	QM,
	MxxA,
	ER,
	VxxD,
	IQ,
	VxxH,
	AS,
	CxxL,
	AV,
	RM,
	LY,
	RN,
	RxxG,
	GxxA,
	QxxE,
	IxxH,
	SxxR,
	pChargeSeg1,
	serWhole,
	LT,
	FN,
	PM,
	SG,
	FD,
	YxxA,
	QxxN,
	TT,
	IK,
	SxxY,
	AC,
	MxxG,
	PG,
	NY,
	YxxG,
	CxxE,
	ExxI,
	YA,
	NA,
	CxxG,
	FxxK,
	KV,
	NxxG,
	IxxI,
	IT,
	RE,
	LxxF,
	AxxM,
	ExxD,
	TxxC,
	DxxR,
	PV,
	RF,
	QW,
	CW,
	EG,
	FR,
	GxxD,
	nChargeSeg4,
	RW,
	LW,
	VF,
	MxxR,
	SxxP,
	LS,
	SA,
	proSeg4,
	DG,
	MxxI,
	CxxM,
	PT,
	CA,
	AxxD,
	PQ,
	NP,
	KxxL,
	IxxW,
	hMomentSeg3,
	NxxI,
	PL,
	ExxR,
	LM,
	KxxE,
	pChargeSeg5,
	KY,
	QxxQ,
	serSeg2,
	HHHBPH,
	proSeg3,
	NI,
	N,
	QL,
	glySeg6,
	QxxH,
	MM,
	YxxN,
	NR,
	LxxK,
	QxxM,
	ExxY,
	MR,
	SxxH,
	TxxS,
	NxxE,
	glySeg2,
	VxxN,
	GxxV,
	amphiSeg1,
	RxxM,
	VG,
	CI,
	RL,
	RxxF,
	WxxG,
	QR,
	GE,
	FxxY,
	QD,
	CxxF,
	V,
	YW,
	DxxP,
	PxxQ,
	FxxG,
	hydroSeg4,
	PN,
	PHHBPH,
	cleavageScore1,
	WF,
	IN,
	AR,
	GP,
	AW,
	DxxD,
	MC,
	HN,
	AxxS,
	MxxL,
	DA,
	TxxK,
	YN,
	YG,
	LxxR,
	pChargeSeg6,
	QxxV,
	HT,
	RxxN,
	ED,
	LxxY,
	RP,
	KxxF,
	VxxE,
	IxxV,
	SS,
	GF,
	PxxK,
	KW,
	TC,
	CP,
	NxxH,
	SF,
	HxxL,
	nChargeSeg5,
	ES,
	PxxT,
	AxxN,
	AxxR,
	YI,
	WW,
	HE,
	NS,
	thrSeg3,
	AxxC,
	FY,
	HBHHbB,
	DH,
	TR,
	YxxC,
	DxxQ,
	TxxM,
	PxxA,
	serSeg6,
	SI,
	glySeg5,
	VT,
	aromaSeg2,
	SM,
	VK,
	AK,
	ML,
	YxxP,
	MxxN,
	LxxL,
	NxxW,
	PD,
	WI,
	QxxL,
	RK,
	W,
	HxxD,
	ExxC,
	CxxQ,
	FP,
	CY,
	hydroSeg5,
	QxxY,
	IE,
	SxxW,
	TW,
	KP,
	KxxW,
	DxxC,
	PF,
	QC,
	RxxL,
	VxxM,
	DN,
	IxxP,
	GxxP,
	QxxC,
	IW,
	A,
	proSeg6,
	FF,
	GQ,
	MxxH,
	aromaWhole,
	IxxQ,
	BHHPPP,
	PxxR,
	II,
	KA,
	MxxK,
	QK,
	RD,
	GxxW,
	KxxQ,
	HM
	];
}

1;
