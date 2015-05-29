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
	int start;
	uint32_t no_vcf:16, keep_flt:16;
	bcf_hdr_t *h;
} bcf_atombuf_t;

void bcf_atomize_first(const bcf1_t *b, bcf_atom_t *a);
void bcf_atomize(const bcf_hdr_t *h, bcf1_t *b, bcf_atom_v *a);
bcf_atombuf_t *bcf_atombuf_init(htsFile *in, int keep_flt);
void bcf_atombuf_destroy(bcf_atombuf_t *buf);
const bcf_atom_t *bcf_atom_read(bcf_atombuf_t *buf);
void bcf_atom2bcf(const bcf_atom_t *a, bcf1_t *b, int write_M, int id_GT);
void bcf_atom_print(const bcf_hdr_t *h, int n, const bcf_atom_t *aa);

static inline int bcf_atom_cmp(const bcf_atom_t *a, const bcf_atom_t *b)
{
	if (a->rid != b->rid) return a->rid - b->rid;
	if (a->pos != b->pos) return a->pos - b->pos;
	if (a->rlen!=b->rlen) return a->rlen-b->rlen;
	return strcmp(a->alt, b->alt);
}

#endif
