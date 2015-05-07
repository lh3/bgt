***WARNING:*** *This is a tool in the alpha phase. Use with caution!*

## Getting started
```sh
# installation
git clone https://github.com/lh3/bgt.git
cd bgt; make
# download precompiled k8 javascript shell (on Mac, replace k8-linux with k8-darwin)
wget -O- http://bit.ly/k8-021 | tar -jxf - k8-linux
# import; vcf2ucf.js turns a multi-allelic VCF to unary VCF
./k8-linux vcf2ucf.js /path/to/srt.vcf.gz | ./bgt import -St ref.fa.fai - prefix
# export
./bgt view -s:sample1,sample2 -r 1:100-200 prefix > sub.vcf
```

## Introduction

BGT is an alternative binary representation of unary VCF where multi-allelic
records are split into multiple lines with additional alleles represented by
the `<M>` symbolic allele. BGT imports from and exports to VCF/BCF *without
per-genotype metadata*. It is smaller and more efficient to process than BCF.

Internally, BGT uses BCF to keep per-allele annotations and PBWT to store
sample genotypes. BGT is able to extract a small subset of samples in time
proportional to the compressed size of genotype matrix, which is in practice
much smaller than the original matrix. This is a significant improvement over
the quadratic time complexity with BCF. BGT potentially supports advanced PBWT
functionailities (e.g. haplotype matching, phasing and imputation), but these
are not on the plan for now.

## Random comments

*This section needs to be cleaned later...*

### Performance

The experiment is done on the 1000g autosomal VCF version v5a. There are
81,708,305 alleles in VCF. The size of the original genotype-only VCFs is
14.4GB (1GB=1024\*1024\*1024 bytes) in size. The corresponding BCF is sized
12.3GB. Running `vcf2ucf.js` on the raw VCF results in a unary BCF 12.4GB in
size, or 163.4 bytes per allele. For this dataset, BGT uses 65.9 bytes/allele,
about 2.5-fold as small as BCF. Note that BGT keeps site-only information in
BCF which uses 27.3 bytes/allele. It uses only 38.6 (=65.9-27.3) bytes/allele
for the genotype matrix.

The following shows the time on a few operations:
```sh
# list all genotypes on chr11
2m10s    bgt view -bl0 -r 11 1000g.bgt > /dev/null
1m19s    htsbox vcfview -bl0 1000g.bcf 11 > /dev/null
# list two samples across whole autosomes
1m51s    bgt view -bl0 -s:HG00100,HG00121 1000g.bgt > /dev/null
8m14s    htsbox vcfview -bl0 -s:HG00100,HG00121 1000g.bcf > /dev/null
```
Note again that for a few samples, BGT spends significant time on inflating
site annotations.  If we only retrieve genotypes (with `pbfview -c 6 -c 7 -c 46
-c 47`), it only takes 46 seconds.

Generally, on decoding genotypes of all samples, BGT has the same time
complexity as BCF but is associated with a larger constant, probably because
BGT incurs many cache misses. Nonetheless, BGT is times smaller than BCF and
faster to retrieve a subset of samples. The difference should be more prominent
given more samples because BGT has better space complexity and time complexity
on sample subsetting.

