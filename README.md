## <a name="started"></a>Getting Started
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
# Server and client (requiring the Go compiler)
go build bgt-server.go
GOMAXPROCS=4 ./bgt-server -d anno11-1M.fmf.gz 1kg11-1M.bgt 2> server.log &
curl -s '0.0.0.0:8000' | less -S  # help
curl -s '0.0.0.0:8000/?a=(impact=="HIGH")&s=(population=="FIN")&f=(AC>0)'
```

## Table of Contents

- [Getting Started](#started)
- [Users' Guide](#guide)
  - [Data model overview](#model)
  - [Import](#import)
    - [Import genotypes](#igenotype)
    - [Import sample phenotypes](#iphenotype)
    - [Import site annotations](#isite)
  - [Query](#query)
    - [Genotype-independent site selection](#givs)
    - [Genotype-independent sample selection](#giss)
    - [Genotype-dependent site selection](#gdvs)
    - [Tabular output](#tabout)
    - [Miscellaneous output](#miscout)
  - [BGT server](#server)
    - [Privacy](#privacy)
- [Further Notes](#notes)
  - [Other genotype formats](#others)
  - [Performance evaluation](#perf)

## <a name="guide"></a>Users' Guide

BGT is a compact file format for efficiently storing and querying whole-genome
genotypes of tens to hundreds of thousands of samples. It can be considered as
an alternative to genotype-only BCFv2. BGT is more compact in size, more
efficient to process, and more flexible on query.

BGT comes with a command line tool and a web application which largely mirrors
the command line uses. The tool supports expressive and powerful query syntax.
The "Getting Started" section shows a few examples.

### <a name="model"></a>1. Data model overview

BGT models a genotype data set as a matrix of genotypes with rows indexed by
site and columns by sample. Each BGT database keeps a genetype matrix and a
sample annotation file. Site annotations are kept in a separate file which is
intended to be shared across multiple BGT databases. This model is different
from VCF in that VCF 1) keeps sample information in the header and 2) stores
site annotations in INFO along with genotypes which are not meant to be shared
across VCFs.

### <a name="import"></a>2. Import

A BGT database always has a genotype matrix and sample names, which are
acquired from VCF/BCF. Site annotations and sample phenotypes are optional but 
are recommended. Flexible meta data query is a distinguishing feature of BGT.

#### <a name="igenotype"></a>2.1 Import genotypes

```sh
# Import BCFv2
bgt import prefix.bgt in.bcf
# Import VCF with "##contig" header lines
bgt import -S prefix.bgt in.vcf.gz
# Import VCF without "##contig" header lines
bgt import -St ref.fa.fai prefix.bgt in.vcf.gz
```
During import, BGT separates multiple alleles on one VCF line. It discards all
INFO fields and FORMAT fields except GT. See section 2.3 about how to use
variant annotations with BGT.

#### <a name="iphenotype"></a>2.2 Import sample phenotypes

After importing VCF/BCF, BGT generates `prefix.bgt.spl` text file, which for
now only has one column of sample names. You can add pheotype data to this file
in a format like (fields TAB-delimited):
```
sample1   gender:Z:M    height:f:1.73     region:Z:WestEurasia     foo:i:10
sample2   gender:Z:F    height:f:1.64     region:Z:WestEurasia     bar:i:20
```
where each meta annotation takes a format `key:type:value` with `type` being
`Z` for a string, `f` for a real number and `i` for an integer. We call this
format Flat Metadata Format or FMF in brief. You can get samples matching
certain conditions with:
```sh
bgt fmf prefix.bgt.spl 'height>1.7&&region=="WestEurasia"'
bgt fmf prefix.bgt.spl 'mass/height**2>25&&region=="WestEurasia"'
```
You can most common arithmetic and logical operators in the condition.

#### <a name="isite"></a>2.3 Import site annotations

Site annotations are also kept in a FMF file like:
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

### <a name="query"></a>3. Query

A BGT query is composed of output and conditions. The output is VCF by default
or can be a TAB-delimited table if requsted. Conditions include
genotype-independent site selection with option `-r` and `-a` (e.g. variants in
a region), genotype-independent sample selection with option `-s` (e.g. a list
of samples), and genotype-dependent site selection with option `-f` (e.g.
allele frequency among selected samples above a threshold). BGT has limited
support of genotype-dependent sample selection (e.g. samples having an allele).

BGT has an important concept "sample group". On the command line, each option
`-s` creates a sample group. The #-th option `-s` populates a pair of `AC#` and
`AN#` aggregate variables. These variables can be used in output or
genotype-dependent site selection.

#### <a name="givs"></a>3.1 Genotype-independent site selection

```sh
# Select by a region
bgt view -r 11:100,000-200,000 1kg11-1M.bgt > out.vcf
# Select by regions in a BED (BGT will read through the entire BGT)
bgt view -B regions.bed 1kg11-1M.bgt > out.vcf
# Select a list of alleles (if on same chr, use random access)
bgt view -a,11:151344:1:G,11:110992:AACTT:A,11:160513::G 1kg11-1M.bgt
# Select by annotations (-d specifies the site annotation database)
bgt view -d anno11-1M.fmf.gz -a'impact=="HIGH"' -CG 1kg11-1M.bgt
```
It should be noted that in the last command line, BGT will read through the
entire annotation file to find the list of matching alleles. It may take
several minutes if the site annotation files contains 100 million lines.
That is why we recommend to use a subset of important alleles (section 2.3).

#### <a name="giss"></a>3.2 Genotype-independent sample selection

```sh
# Select a list of samples
bgt view -s,HG00171,HG00173 1kg11-1M.bgt
# Select by phenotypes (see also section 2.2)
bgt view -s'population=="CEU"' 1kg11-1M.bgt
# Create sample groups (there will be AC1/AN1 and AC2/AN2 in VCF INFO)
bgt view -s'population=="CEU"' -s'population=="YRI"' -G 1kg11-1M.bgt
```

#### <a name="gdvs"></a>3.3 Genotype-dependent site selection

```sh
# Select by allele frequency
bgt view -f'AN>0&&AC/AN>.05' 1kg11-1M.bgt
# Select by group frequnecy
bgt view -s'population=="CEU"' -s'population=="YRI"' -f'AC1>10&&AC2==0' -G 1kg11-1M.bgt
```
Of course, we can mix all the three types of conditions in one command line:
```sh
bgt view -G -s'population=="CEU"' -s'population=="YRI"' -f'AC1/AN1>.1&&AC2==0' \
         -r 11:100,000-500,000 -d anno11-1M.fmf.gz -a'CDSpos>0' 1kg11-1M.bgt
```

#### <a name="tabout"></a>3.4 Tabular output

```sh
# Output position, sequence and allele counts
bgt view -t CHROM,POS,REF,ALT,AC1,AC2 -s'population=="CEU"' -s'population=="YRI"' 1kg11-1M.bgt
```

#### <a name="miscout"></a>3.5 Miscellaneous output

```sh
# Get samples having a set of alleles (option -S)
bgt view -S -a,11:151344:1:G,11:110992:AACTT:A,11:160513::G -s'population=="CEU"' 1kg11-1M.bgt
# Count haplotypes
bgt view -Hd anno11-1M.fmf.gz -a'gene=="SIRT3"' -f 'AC/AN>.01' 1kg11-1M.bgt
# Count haplotypes in multiple populations
bgt view -Hd anno11-1M.fmf.gz -a'gene=="SIRT3"' -f 'AC/AN>.01' \
         -s'region=="Africa"' -s'region=="EastAsia"' 1kg11-1M.bgt
```

### <a name="server"></a>4. BGT server

In addition to a command line tool, we also provide a prototype web application
for genotype query. The query syntax is similar to `bgt view` as is shown in
"Getting Started", but with some notable differences:

1. The server uses `.and.` for the logical AND operator `&&` (as `&` is a special character to HTML).
2. The server can't load a list of samples from a local file (for security).
3. The server doesn't support BCF output for now (can be implemented on request).
4. The server doesn't output genotypes by default (option `g` required for server).
5. The server loads site annotations into RAM (for real-time response but requiring more memory).
6. By default (tunable), the server processes up to 10 million genotypes and then truncates the result.
7. The server may forbid the output of genotypes of some samples (see below).

#### <a name="privacy"></a>4.1 Privacy

The BGT server implements a simple mechanism to keep the privacy of samples or
a subset of samples. It is controlled by a single parameter: minimal sample
group size or MGS.  The server refuses to create a sample group if the size of
the group is smaller than the MGS of one of the samples in the group. In
particular, if MGS is above one, the server doesn't report sample name or
sample genotypes.  Each sample may have a different MGS as is marked by the
`_mgs` integer tag in `prefix.bgt.spl`. For samples without this tag, a
default MGS is applied.

## <a name="notes"></a>Further Notes

### <a name="others"></a>Other genotype formats

* BGT vs [PBWT][pbwt]. BGT uses the same data structure as PBWT and is inspired
  by PBWT. PBWT supports advanced query such as haplotype matching, phasing
  and imputation, while BGT puts more emphasis on fast random access and data
  retrieval.

* BGT vs [BCF2][vcf]. BCF is more versatile. It is able to keep per-genotype
  meta information (e.g. per-genotype read depth and genotype likelihood). BGT
  is generally more efficient and times smaller. It scales better to many
  samples. BGT also supports more flexible queries, although technically,
  nothing prevents us from implementing similar functionalities on top of BCF.

* BGT vs [GQT][gqt]. GQT should be much faster on traversing sites across whole
  chromosomes without considering LD. It is however inefficient to retrieve
  data in small regions or to get haplotype information due to its design.
  For this reason, GQT is regarded as a complement to BCF or BGT, not a
  replacement. On file size, GQT is usually larger than genotype-only BCF and
  is thus larger than BGT.

### <a name="perf"></a>Performance evaluation

The test is run on the first release of [Haplotype Reference Consortium][hrc]
(HRC) data. There are ~39 million phased SNPs in 32,488 samples. We have
generated the BGT for the entire dataset, but We are only running tools in
region chr11:10,000,001-20,000,000. The following table shows the time and
command line. Note that the table omits option `-r 11:10,000,001-20,000,000`
which has been applied to all command lines below.

|Time   |Command line|
|------:|:------------|
|11s    |bgt view -G HRC-r1.bgt|
|13s    |bcftools view -Gu HRC-r1.bcf|
|30s    |bgt view -GC HRC-r1.bgt|
|4s     |bgt view -GC -s'source=="1000G"'|
|19s    |bcftools view -Gu -S 1000G.txt HRC-r1.bcf|
|8s     |bgt view -G -s 'source=="UK10K"' -s 'source=="1000G"&&population!="GBK"'|

On file sizes, the BGT database for HRC-r1 is 7.4GB (this excludes the `.pb1`
file which is generated but not used for now; 1GB=1024\*1024\*1024 bytes). In comparison,
BCFv2 for the same data takes 65GB, GQT 93GB and PBWT 4.4GB. BGT and PBWT,
which are based on the same data structure, are much more compact. BGT is
larger than PBWT primarily because BGT keeps an extra bit per haplotype to
distinguish reference and multi allele, and stores markers to enable fast
random access.

[hrc]: http://www.haplotype-reference-consortium.org
[gqt]: https://github.com/ryanlayer/gqt
[pbwt]: https://github.com/richarddurbin/pbwt
[vcf]: https://samtools.github.io/hts-specs/
