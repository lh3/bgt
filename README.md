## Getting started
```sh
# installation
git clone https://github.com/lh3/bgt.git
cd bgt; make
# import
./bgt import -St ref.fa.fai in.vcf.gz prefix
# export
./bgt view -s:sample1,sample2 -r 1:100-200 prefix > sub.vcf
./bgt mview -C1 -s:sample1,sample3 prefix1 prefix2 > sub.vcf
```

## Introduction

BGT is an alternative binary representation of unary VCF where multi-allelic
records are split into multiple lines with additional alleles represented by
the `<M>` symbolic allele. BGT imports from and exports to VCF/BCF *without
per-genotype metadata*. It is smaller and more efficient to process than BCF.

Internally, BGT uses BCF to keep allele positions and sequences, and uses PBWT to store
sample genotypes. BGT is able to extract a small subset of samples in time
proportional to the compressed size of genotype matrix, which is in practice
much smaller than the original matrix. This is a significant improvement over
the quadratic time complexity with BCF. BGT potentially supports advanced PBWT
functionailities (e.g. haplotype matching, phasing and imputation), but these
are not on the plan for now.

## Random comments

*This section needs to be cleaned later...*

### Performance on 1000g VCF

The experiment is done on the 1000g autosomal VCF version v5a. There are
81.7 million alleles in VCF. The size of BGT is 3.46GB (1GB=1024\*1024\*1024
bytes), or 45.4 bytes per site in average. Of the 45.4 bytes, 6.6 bytes
are used to keep ALT sequences and positions, and 38.8 bytes for genotypes
in the indexed PBWT. Converting this BGT file(s) to BCF results in a 13.1GB
file. BGT is much smaller for the same information.

The following shows the time on a few operations:
```sh
# list all genotypes on chr11
2m09s    bgt view -bl0 -r 11 1000g.bgt > /dev/null
1m12s    htsbox vcfview -bl0 1000g.bcf 11 > /dev/null
# extract the genotypes of the 99 CEU samples on chr11
14s      bgt mview -bl0 -s CEU.txt -r 11 1000g.bgt > /dev/null
22s      htsbox vcfview -bl0 -s CEU.txt 1000g.bcf 11 > /dev/null
# extract two samples across all autosomes
1m10s    bgt view -bl0 -s:HG00100,HG00121 1000g.bgt > /dev/null
5m45s    htsbox vcfview -bl0 -s:HG00100,HG00121 1000g.bcf > /dev/null
# extract chr11 alleles polymorphic in CEU; no sample genotypes
12s      bgt mview -aGC1 -bl0 -s CEU.txt -r 11 1000g.bgt > /dev/null
```

Generally, on decoding genotypes of all samples, BGT has the same time
complexity as BCF but is associated with a larger constant, probably because
BGT incurs many cache misses. Nonetheless, BGT is times smaller than BCF and
faster to retrieve a subset of samples. The difference should be more prominent
given more samples because BGT has better space complexity and time complexity
on sample subsetting.

### Performance on ExAC VCF

ExAC release 3 consists of 60 thousand samples and 9.3 million filtered alleles.
The BGT file size is 7.4GB, or 856 bytes per site. The genotype compression
ratio is much worse than 1000g because: 1) ExAC contains many missing data; 2)
ExAC are unphased and 3) exons separated by long introns or intergenes are
unlinked. Nonetheless, the BGT file is still a thrid of the equivalent BCF file
converted from BGT (accurate numbers to appear).

The following shows the time on a few operations (there are 565k alleles on chr11):
```sh
# extract the genotypes of the 74 1000g CEU samples on chr11
3.2s     bgt mview -bl0 -s CEU.txt -r 11 exac3.bgt > /dev/null
# get chr11 alleles polymorphic in the 74 CEU samples
2.9s     bgt mview -GC1 -s CEU.txt -r 11 exac3.bgt > /dev/null
# extract the genotypes of 1000 random samples on chr11
19.2s    bgt mview -bl0 -s random-1000.txt -r 11 exac3.bgt > /dev/null
# extract the genotypes of 10000 random samples on chr11
3m15s    bgt mview -bl0 -s random-1000.txt -r 11 exac3.bgt > /dev/null
```
