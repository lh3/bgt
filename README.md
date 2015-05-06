**WARNING:** This is a tool in the alpha phase. Use with caution!

## Getting started
```sh
git clone https://github.com/lh3/bgt.git
cd bgt; make
wget -O- http://bit.ly/k8-021 | tar -jxf  -
./k8-linux vcf2ucf.js srt.vcf.gz | ./bgt ucf2bgt -S - prefix
./bgt view -s:sample1,sample2 -r 1:100-200 prefix > sub.vcf
```

## What is BGT?

BGT is a PBWT-based binary representation of unary VCF without per-sample
metadata. It imports from and exports to VCF/BCF. In comparison to BCF, BGT is
smaller, much faster to extract a subset of samples and potentially supports
more advanced haplotype-based queries.

(more explanations coming soon...)
