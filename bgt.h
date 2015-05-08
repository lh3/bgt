#ifndef BGT_H
#define BGT_H

#include "vcf.h"
#include "pbwt.h"

typedef struct {
	int n_samples;
	char **samples;
	void *h_samples; // hash table
	pbf_t *pb;
	htsFile *bcf;
	bcf_hdr_t *h0; // site-only BCF header
	bcf1_t *b0; // site-only BCF record
	hts_idx_t *idx; // BCF index
	hts_itr_t *itr;

	int n_sub, *sub;
	bcf_hdr_t *h_sub;
} bgt_t;

typedef struct {
	bcf1_t *b0;
	uint8_t *a[2];
} bgt_rec_t;

typedef struct {
	int rid, pos;
	int m_b, n_b, row, finished;
	bgt_rec_t *b;
} bgt_pos_t;

typedef struct {
	int n_bgt;
	bgt_t **bgt;
	bgt_pos_t *p;
	bcf_hdr_t *h;
} bgtm_t;

bgt_t *bgt_open(const char *prefix);
void bgt_close(bgt_t *bgt);
void bgt_set_samples(bgt_t *bgt, int n, char *const* samples);
int bgt_set_region(bgt_t *bgt, const char *reg);
int bgt_read(bgt_t *bgt, bcf1_t *b);

bgtm_t *bgtm_open(int n_files, char *const*fns);
void bgtm_close(bgtm_t *bm);
void bgtm_set_samples(bgtm_t *bm, int n, char *const* samples);
int bgtm_set_region(bgtm_t *bm, const char *reg);

/******************************
 * The missing BCF operations *
 ******************************/

#include <string.h>

static inline void bcf_copy(bcf1_t *dst, const bcf1_t *src)
{
	kstring_t ts = dst->shared, ti = dst->indiv;
	free(dst->d.id); free(dst->d.allele); free(dst->d.flt); free(dst->d.info); free(dst->d.fmt);
	*dst = *src;
	memset(&dst->d, 0, sizeof(bcf_dec_t));
	dst->unpacked = 0;
	dst->unpack_ptr = 0;
	ts.l = ti.l = 0;
	dst->shared = ts; dst->indiv = ti;
	kputsn(src->shared.s, src->shared.l, &dst->shared);
	kputsn(src->indiv.s, src->indiv.l, &dst->indiv);
}

static inline int bcf_cmp(const bcf1_t *a, const bcf1_t *b)
{
	int i, l[2];
	uint8_t *ptr[2];
	if (a->rid != b->rid) return a->rid - b->rid;
	if (a->pos != b->pos) return a->pos - b->pos;
	if (a->rlen!=b->rlen) return a->rlen-b->rlen;
	for (i = 0; i < 2; ++i) {
		int x, type;
		ptr[i] = (uint8_t*)a->shared.s;
		x = bcf_dec_size(ptr[i], &ptr[i], &type); // size of ID
		ptr[i] += x << bcf_type_shift[type]; // skip ID
		x = bcf_dec_size(ptr[i], &ptr[i], &type); // size of REF
		ptr[i] += x << bcf_type_shift[type]; // skip REF
		l[i] = bcf_dec_size(ptr[i], &ptr[i], &type); // size of ALT1
	}
	if (l[0] != l[1]) return l[0] - l[1];
	return strncmp((char*)ptr[0], (char*)ptr[1], l[0]);
}

#endif
