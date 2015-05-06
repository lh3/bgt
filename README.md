**WARNING:** This is a tool in the alpha phase. Use with caution!

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
per-sample metadata*. It is smaller and more efficient to process than BCF.

Internally, BGT uses BCF to keep per-allele annotations and PBWT to store
sample genotypes. BGT is able to extract a small subset of samples in time
proportional to the compressed size of genotype matrix, which is in practice
much smaller than the original matrix. This is a significant improvement over
the quadratic time complexity with BCF. BGT potentially supports advanced PBWT
functionailities (e.g. haplotype matching, phasing and imputation), but these
have not been implemented yet.
