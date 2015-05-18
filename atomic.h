#ifndef ATOMIC_H
#define ATOMIC_H

#include "vcf.h"

typedef struct {
	kstring_t ref;
	char *alt;
	int rid, pos, rlen, anum;
	uint8_t *gt;
} bcf_atom_t;

typedef struct {
	int n, m;
	bcf_atom_t *a;
} bcf_atom_v;

typedef struct {
	bcf_atom_v a, tmp;
	bcf1_t *b;
} bcf_atom_reader;

void bcf_atomize(const bcf_hdr_t *h, bcf1_t *b, bcf_atom_v *a);

static inline int bcf_atom_cmp(const bcf_atom_t *a, const bcf_atom_t *b)
{
	if (a->rid != b->rid) return a->rid - b->rid;
	if (a->pos != b->pos) return a->pos - b->pos;
	if (a->rlen!=b->rlen) return a->rlen-b->rlen;
	return strcmp(a->alt, b->alt);
}

#endif
