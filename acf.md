# The Atomic VCF Format (ACF)

The Atomic VCF format (ACF) is a strict subset of VCF to describe genetic
sequence variations in a population. Unlike VCF which may encode multiple
variants on one line (aka record), ACF only encodes one atomic allele per
record. The representation for a set of variants is unique in ACF. BGT always
converts VCF to ACF.

## Definitions

Given a reference genome, an *allele* is a 4-tuple (`ctg`,`pos`,`len`,`seq`),
indicating sequence `seq` replaces a `len`-long reference subsequence at
`ctg`:`pos`.  An allele is a *reference allele* (or REF) if `seq` is the same
as the reference subsequence; otherwise the allele is an *alternate allele* (or
ALT). An ALT is *atomic* if it is a single substitution, a single insertion or
a single deletion.

In ACF, each record encodes one atomic ALT. For a particular allele *X*, if
there are other allele(s) overlapping with *X*, these other allele(s) will be
indicated by a symbolic allele `<M>`. `<M>` plays the same role as the `*`
allele in VCF for deletions, but it is more general. Because we use `<M>` for
all types of other overlapping alleles, in ACF, the allele number in the GT
field can only be `.`, `0`, `1` or `2`. In BGT, these four allele numbers are
encoded with 2 bits.

## Examples

### A simple example

```txt
Pos: 12345678901234567890123 4567890123
Ref: XXXXXXXXXCATATGCAAGTCGT-TATTAGAGCTXXXXX
H1:  XXXXXXXXXCATGTGC--GTCGTATATT----CTXXXXX
H2:  XXXXXXXXXCATATGCAAGTCGTATATT--AGCTXXXXX
```
The corresponding ACF is
```txt
chr1  13  A      G      ..  1|0
chr1  16  CAA    C      ..  1|0
chr1  23  T      TA     ..  1|1
chr1  27  TAGAG  T,<M>  ..  1|2
chr1  27  TAG    T,<M>  ..  2|1
```
On the first three records, there are no overlapping variants. ACF is identical
VCF in this case. There are two overlapping deletions at pos 27. In ACF, we
describe each separately and use `<M>` as a placeholder for alleles not
described on the line. When we parse ACF at pos 27, we have to read both
records to reconstruct the underlying alignment.

### A contrived example

```txt
Pos: 123456789012345678 90
Ref: XXXXXXXXXGTATATAGC-GAXXXXX
H1:  XXXXXXXXXGTATA-------XXXXX
H2:  XXXXXXXXXG------GCTGAXXXXX
```
Can be encoded in ACF as
```txt
chr1  10  GTATATA  G,<M>   ..  2|1
chr1  14  ATAGCGA  A,<M>   ..  1|2
chr1  18  C        CT,<M>  ..  2|1
```
Similarly, an ACF parser has to memorize all three records to correctly
reconstruct the haplotypes.

### A multi-sample example

The following VCF
```txt
11  101  GCGT  G,GCGA,GTGA,CCGT ..  0|1 1|2 2|3 2|4
```
can be converted with `bgt atomize -MS` to ACF
```txt
11  101  G     C,<M>  ..  0|2   2|0   0|0  0|1
11  101  GCGT  G,<M>  ..  0|1   1|2   2|2  2|2
11  102  C     T,<M>  ..  0|2   2|0   0|1  0|0
11  104  T     A,<M>  ..  0|2   2|1   1|1  1|0
```

## An Extension to ACF

A major issue with ACF is its inconenience: when there are overlapping
variants, we don't know what `<M>` or allele number `2` refers to. We have to
combine multiple records to reconstruct genotypes. A potential workaround is to
introduce VCF INFO and FORMAT key-value pairs to spell out `<M>`. For the 2nd
example, we can
```txt
chr1  10  GTATATA  G,<M>   A2=chr1_14_7_A               GT:A2  2|1:3|1
chr1  14  ATAGCGA  A,<M>   A2=chr1_10_7_G,chr1_18_1_CT  GT:A2  1|2:1|3,4
chr1  18  C        CT,<M>  A2=chr1_14_7_A               GT:A2  2|1:3|1
```
Here the `A2` INFO tag uses a succinct way to encode the list of overlapping
alleles; the `A2` FORMAT tag indexes into the combine ALT and A2 array and
gives the complete genotype. However, this extension leads to duplicated
information which makes it more likely to introduce errors. We can't go far
without breaking the compatibility with VCF.
