## Getting started
```sh
# installation
git clone https://github.com/lh3/bgt.git
cd bgt; make
# import
./bgt import -St ref.fa.fai in.vcf.gz prefix
# export
./bgt view -f"AC>0" -s:sample1,sample2 -r 1:100-200 prefix > sub.vcf
./bgt view -f"AN>=100&&AC/AN>.05" -s"gender=='M'" prefix1 prefix2 > sub.vcf
```

## Introduction

BGT is an alternative to BCF2 for storing multi-sample genotypes *without
per-genotype metadata*. It is much smaller than BCF2 and often faster for
accessing a subset of samples. It also supports more types of genotype queries.
Here are a few examples demonstrating how to query with BGT. In these examples,
we assume users have attached per-sample annotation to `pre.spl` in a format
like:
```
HG00096  gender:Z:M  population:Z:GBR
HG00315  gender:Z:F  population:Z:FIN
NA19144  gender:Z:M  population:Z:YRI
```
```sh
# GBR samples with ALT allele frequency over 5%
bgt view -G -f"AN>=20&&AC/AN>.05" -s"population=='GBR'" pre
# alleles present in both GBR and FIN; -g defines sample groups and populates AC1 etc
bgt view -s"population=='GBR'||population=='FIN'" -g"population=='GBR'" \
    -g"population=='FIN'" -Gf'AC1>0&&AC2>0' pre
# a shorter CMD for the above
bgt view -s"population=='GBR'||population=='FIN'" -g"population=='GBR'" \
    -Gf'AC0>0&&AC1>0' pre
# alleles >20% in GBR but <5% in YRI
bgt view -s"population=='GBR'||population=='YRI'" -g"population=='GBR'" \
    -g"population=='YRI'" -Gf'AC1/AN1>.2&&AC2/AN2<.05' pre
```

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
28s      bgt view -s"population=='GBR'||population=='YRI'" -g"population=='GBR'" \
         -g"population=='YRI'" -Gf'AC1/AN1>.2&&AC2/AN2<.05' -r11 1000g-p3v5a.bgt > /dev/null
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
