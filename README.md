# MitoFates (ver. 1.2) 
![](image/MitoFatesTitle.png)


MitoFates predicts presequences (cleavable mitochondrial localization signal 
in N-termini, and their cleaved position.

Web server of MitoFates is available at
[http://mitf.cbrc.jp/MitoFates/](https://mitf.cbrc.pj.aist.go.jp/MitoFates)

The source code is available at 
[https://github.com/ken1roimai/MitoFates](https://github.com/ken1roimai/MitoFates)

Here are some very brief notes on using the MitoFates software.

## Installation
MitoFates is perl script, thus instration of MitoFates is not necessary.
but following tools and Perl modules is required.

### Required tools and Perl modules
Please make sure installation of libsvm and required perl modules.
1. LIBSVM ver 3.0 (or later)
2. Math::Cephes (perl library)
3. Perl6::Slurp (perl library)
4. Inline::C    (perl library)

To install libsvm 3.0, please follow their installation guide.
Using Other versions of libsvm, Mitofates seems to work. But we cannot gurantee correct prediction.
In addition, make sure their binary files are located in a directory and make path to the directory.
Perl modules are available at CPAN (http://www.cpan.org/).

Inline::C is used in /bin/module/DirichletRegulator_fast.pm.
When you encounter a problem installing this module, you have two options:

1. Consult with other web sites or perl users to install them correctly.
2. Replace this module with the previous pure perl module named DirichletRegulator.pm.
    All you need is delete "_fast" in line 31 of /bin/cleavage.pl and line 22 of MitoFatesPWMUtil.pm.
    This makes program slow, however you can skip annoying errors related to installation above.

## Usage

[MitoFates.pl] [FastaFile] [Organism Flag: fungi, metazoa or plant]
Only fasta file format is currently supported.

### Example,
    > perl MitoFates.pl example.fasta fungi

## Output format
MitoFates outputs 21 Columns with the header (see below).

Column 1    :Sequence ID

Column 2    :Probability of presequence (default cutoff: 0.385).

Column 3    :Prediction result (e.g. Possessing mitochondrial presequence or No mitochondrial presequence)

Column 4    :Predicted cleavage sites for processing enzymes (MPP, Oct1 and Icp55)(e.g. 23(MPP), 31(Oct1)).

Column 5    :Net charge of predcited presequence

Column 6    :Sequence region matching to TOM20 recognition motif (e.g. 15-19).
             In case of multiple hits, each region is deliminated by comma.
             
Column 7    :Region having the maximum score of positively charged amphiphilicity (PA) score in 
             N-terminal 30 residues. High and low means that the score is higher and lower than 
             the sensitivity and specificity versus cutoff value, respectively (e.g. 2-11-high or 2-11-low).
             
Column 8-21 :Regions matching to any of 14 types of statistically significant 6-mers in N-terminal 
             30 residues of presequence (e.g. 13-18). In case of multiple hits, each region is 
             deliminated by comma; 14 motifs is described below. 
             
             BHHPPP
             BPHBHH
             HBHHBb
             HBHHbB
             HHBHHB
             HHBPHB
             HHBPHH
             HHBPHP
             HHHBBH
             HHHBPH
             HHHHBB
             HHPBHH
             HPBHHP
             PHHBPH
The letters of the 14 motifs are defined as follows; H, B, P and g indicate
hydrophobic (L, F, I, V, W, Y, M, C, A), basic (R, K, H), polar (S, T, N, Q),
and secondary structure breaker (P, G), respectively. 

## Additional Information
If you use MitoFates for your research, please cite our paper.
https://doi.org/10.1074%2Fmcp.M114.043083

