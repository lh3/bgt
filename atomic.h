#ifndef ATOMIC_H
#define ATOMIC_H

#include "vcf.h"

typedef struct {
	kstring_t ref;
	char *alt;
	int rid, pos, rlen, n_gt;
	uint32_t anum:30, has_multi:1, from_new:1;
	uint8_t *gt;
} bcf_atom_t;

typedef struct {
	int n, m;
	bcf_atom_t *a;
} bcf_atom_v;

typedef struct {
	htsFile *in;
	bcf_atom_v a;
	bcf1_t *b;
	int start, no_vcf;
	bcf_hdr_t *h;
} bcf_atombuf_t;

void bcf_atomize(const bcf_hdr_t *h, bcf1_t *b, bcf_atom_v *a);
bcf_atombuf_t *bcf_atombuf_init(htsFile *in);
void bcf_atombuf_destroy(bcf_atombuf_t *buf);
const bcf_atom_t *bcf_atom_read(bcf_atombuf_t *buf);

static inline int bcf_atom_cmp(const bcf_atom_t *a, const bcf_atom_t *b)
{
	int ret;
	if (a->rid != b->rid) return a->rid - b->rid;
	if (a->pos != b->pos) return a->pos - b->pos;
	if (a->rlen!=b->rlen) return a->rlen-b->rlen;
	ret = strcmp(a->alt, b->alt);
	if (ret != 0) return ret;
	else return (int)a->from_new - (int)b->from_new;
}

#endif
