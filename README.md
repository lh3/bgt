## Getting started
```sh
# Installation
git clone https://github.com/lh3/bgt.git
cd bgt; make
# Download demo BCF (1st 1Mbp of chr11 from 1000g), and convert to BGT
wget -O- http://bit.ly/BGTdemo | tar xf -
./bgt import 1kg11-1M.bgt 1kg11-1M.raw.bcf
gzip -dc 1kg11-1M.raw.samples.gz > 1kg11-1M.bgt.spl  # sample meta data
# Get all sample genotypes
./bgt view -C 1kg11-1M.bgt | less -S
# Get genotypes of HG00171 and HG00173 in region 11:100,000-200,000
./bgt view -s,HG00171,HG00173 -f'AC>0' -r 11:100000-200000 1kg11-1M.bgt
# Get alleles high-frequency in CEU but absent from YRI
./bgt view -s'population=="CEU"' -s'population=="YRI"' -f'AC1/AN1>=0.1&&AC2==0' -G 1kg11-1M.bgt
# Select high-impact sites (var annotation provided with -d)
./bgt view -d anno11-1M.fmf.gz -a'impact=="HIGH"' -CG 1kg11-1M.bgt
# Server and client
go build bgt-server.go
./bgt-server -d anno11-1M.fmf.gz 1kg11-1M.bgt 2> server.log &
curl -s '0.0.0.0:8000' | less -S  # help
curl -s '0.0.0.0:8000/?a=(impact=="HIGH")&s=(population=="FIN")&f=(AC>0)'
```

## Introduction

BGT is a compact file format for efficiently storing and querying whole-genome
genotypes of tens to hundreds of samples. It can be considered as an alternative
to BCFv2 but it keeps genotypes only, *without per-genotype metadata*. In
comparison to BCFv2, BGT is more compact, more efficient for sample subsetting,
and supports sample metadata and allele annotations.

BGT comes with a command line tool and a web application which largely mirrors
the command line uses. It supports expressive and powerful query syntax. The
"Getting started" section shows a few examples.

## Users' guide

### 1. Import genotypes

#### 1.1 Import VCF/BCF

```sh
# Import BCFv2
bgt import prefix.bgt in.bcf
# Import VCF with "##contig" header lines
bgt import -S prefix.bgt in.vcf.gz
# Import VCF without "##contig" header lines
bgt import -St ref.fa.fai prefix.bgt in.vcf.gz
```

#### 1.2 Import sample pheonotypes

After importing VCF/BCF, BGT generates `prefix.bgt.spl` text file. You can add
pheotype data to this file in a format like (fields TAB-delimited):
```
sample1   gender:Z:M    height:f:1.73     region:Z:WestEurasia     foo:i:10
sample2   gender:Z:F    height:f:1.64     region:Z:WestEurasia     foo:i:20
```
where each meta annotation takes a format `key:type:value` with `type` being
`Z` for a string, `f` for a real number and `i` for an integer. We call this
format as Flat Metadata Format or FMF in brief. You can get samples matching
certain conditions with:
```sh
bgt fmf prefix.bgt.spl 'height>1.7&&region=="WestEurasia"'
bgt fmf prefix.bgt.spl 'mass/height**2>25&&region=="WestEurasia"'
```
You can most common arithmetic and logical operators in the condition.

#### 1.3 Import variant annotations

Each VCF keeps variant annotations in its own INFO field. BGT is different. It
recommends to keep annotations of all important variants in a separate FMF file
and share it across BGT files. This FMF looks like:
```
11:209621:1:T  effect:Z:missense_variant   gene:Z:RIC8A  CCDS:Z:CCDS7690.1  CDSpos:i:347
11:209706:1:T  effect:Z:synonymous_variant gene:Z:RIC8A  CCDS:Z:CCDS7690.1  CDSpos:i:432
```
We provide a script `misc/vep2fmf.pl` to convert the VEP output with the
`--pick` option to FMF.

Note that due to an implementation limitation, we recommend to use a subset of
"important" variants with BGT, for example:
```sh
gzip -dc vep-all.fmf.gz | grep -v "effect:Z:MODIFIER" | gzip > vep-important.fmf.gz
```
Using the full set of variants is fine, but is much slower with the current
implementation.

### 2. Query genotypes or genotype summary statistics

#### 2.1 General query patterns


## Discussions

### Unary VCF

During import, BGT converts a multi-allelic VCF to unary VCF on the fly. In
this process, BGT first decompose a complex variant consisting of multiple SNPs
and INDELs into atomic events and then describes each atomic event with one and
only one VCF line. If a sample has two overlapping non-reference alleles, BGT
will use a `<M>` symbolic allele. For example, BGT internally converts this VCF:
```
#CHROM POS ID REF ALT     QUAL FILTER INFO FORMAT S1   S2   S3   S4
11  101  .  GCGT  G,GCGA,GTGA,CCGT 199  .  .  GT  0/1  1/2  2/3  2/4
```
to
```
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2 S3 S4
11  101  .  G     C,<M>  0  .  .  GT  0/2  2/0  0/0  0/1
11  101  .  GCGT  G,<M>  0  .  .  GT  0/1  1/2  2/2  2/2
11  102  .  C     T,<M>  0  .  .  GT  0/2  2/0  0/1  0/0
11  104  .  T     A,<M>  0  .  .  GT  0/2  2/1  1/1  1/0
```
While [vt][vt] aims to achieve a similar goal, it is unable to properly set the
genotype fields as it [cannot][vtissue] perform decomposition and "unarization"
at the same time.

### Theoretical advnatages of BGT

Internally, BGT uses BCF to keep allele positions and sequences, and uses PBWT
to store sample genotypes. The time complexity to extract *k* samples across
*n* sites is *O*(*n*(*mH*<sub>0</sub>+*k*)), where *m* is the total number of
samples in BGT and *H*<sub>0</sub> the average 0-order empirical entroy of PBWT. When *k* is
larger than a few hundred, the time complexity approaches *O*(*nk*),
independent of *m*. This is a significant improvement over the *O*(*nm*) time
complexity if we use BCF. Although BGT comes with a larger constant - probably
due to more cache misses - and is thus slower on extracting all samples, it
can be parallelized by extracting a subset of samples and combining samples
later, which can't be achieved with BCF.

In addition to genotype extraction, BGT may in theory support advanced PBWT
operations such as imputation and phasing. Implementing these functionalities
is not on the plan for the time being.

### Performance on 1000g VCF

The experiment is done on the 1000g autosomal VCF version v5a. There are
81.7 million alleles in VCF. The size of BGT is 3.46GB (1GB=1024\*1024\*1024
bytes), or 45.4 bytes per site in average. Of the 45.4 bytes, 6.6 bytes
are used to keep ALT sequences and positions, and 38.8 bytes for genotypes
in the indexed PBWT. Converting this BGT file(s) to BCF results in a 13.1GB
file. BGT is much smaller for the same information.

The following shows the time on a few operations (the following is done with
BGT-r98; newer versions should have a similar performance):
```sh
# extract two samples across all autosomes
1m10s    bgt view -bl0 -s:HG00100,HG00121 1000g.bgt > /dev/null
5m45s    htsbox vcfview -bl0 -s:HG00100,HG00121 1000g.bcf > /dev/null
# extract the genotypes of the 99 CEU samples on chr11
14s      bgt view -bl0 -s CEU.txt -r 11 1000g.bgt > /dev/null
22s      htsbox vcfview -bl0 -s CEU.txt 1000g.bcf 11 > /dev/null
# list all genotypes on chr11
2m09s    bgt view -bl0 -r 11 1000g.bgt > /dev/null
1m12s    htsbox vcfview -bl0 1000g.bcf 11 > /dev/null
# variants common in GBR but rare in YRI
28s      bgt view -s"population=='GBR'" -s"population=='YRI'" \
             -Gf'AC1/AN1>.2&&AC2/AN2<.05' -r11 1000g.bgt > /dev/null
```

### Performance on ExAC VCF

ExAC release 3 consists of 60 thousand samples and 9.3 million filtered alleles.
The BGT file size is 7.4GB, or 856 bytes per site. The genotype compression
ratio is much worse than 1000g because: 1) ExAC contains many missing data; 2)
ExAC are unphased and 3) exons separated by long introns or intergenes are
unlinked. Nonetheless, the BGT file is less than half of the equivalent BCF
file (17GB in size) converted from BGT.

The following shows the time on a few operations (there are 565k alleles on chr11):
```sh
# extract the genotypes of 2 public samples, from chr11
1.6s     bgt view -bl0 -s:NA10851,NA12044 -r11 exac3.bgt > /dev/null
40.3s    htsbox vcfview -bl0 -s:NA10851,NA12044 exac3.bcf 11 > /dev/null
# extract the genotypes of the 74 1000g CEU samples on chr11
3.2s     bgt view -bl0 -s CEU.txt -r 11 exac3.bgt > /dev/null
41.4s    htsbox vcfview -bl0 -s CEU.txt exac3.bcf 11 > /dev/null
# get chr11 alleles polymorphic in the 74 CEU samples
2.9s     bgt view -GC1 -s CEU.txt -r 11 exac3.bgt > /dev/null
# extract the genotypes of 1000 random samples on chr11
19.2s    bgt view -bl0 -s random-1000.txt -r 11 exac3.bgt > /dev/null
47.2s    htsbox vcfview -bl0 -s random-1000.txt exac3.bcf 11 > /dev/null
# extract the genotypes of 10000 random samples on chr11
3m15s    bgt view -bl0 -s random-1000.txt -r 11 exac3.bgt > /dev/null
1m42s    htsbox vcfview -bl0 -s random-10000.txt exac3.bcf 11 > /dev/null
```

[vt]: https://github.com/atks/vt
[vtissue]: https://github.com/atks/vt/issues/26
