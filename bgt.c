#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "bgt.h"
#include "kstring.h"

#include "khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)

#include "ksort.h"
#define generic_key(x) (x)
KRADIX_SORT_INIT(i, int, generic_key, 4)

bgt_t *bgt_open(const char *prefix)
{
	char *fn;
	bgt_t *bgt;
	int i, absent;
	khash_t(s2i) *h;

	bgt = (bgt_t*)calloc(1, sizeof(bgt_t));
	fn = (char*)malloc(strlen(prefix) + 9);

	// read samples
	sprintf(fn, "%s.spl", prefix);
	bgt->samples = hts_readlines(fn, &bgt->n_samples);
	h = kh_init(s2i);
	for (i = 0; i < bgt->n_samples; ++i) {
		khint_t k;
		k = kh_put(s2i, h, bgt->samples[i], &absent);
		assert(absent); // otherwise there are duplicated sample names
		kh_val(h, k) = i;
	}
	bgt->h_samples = h;

	// open PBF
	sprintf(fn, "%s.pbf", prefix);
	bgt->pb = pbf_open_r(fn);

	// open BCF and load index
	sprintf(fn, "%s.bcf", prefix);
	bgt->bcf = vcf_open(fn, "rb", 0);
	bgt->h0 = vcf_hdr_read(bgt->bcf);
	bgt->idx = bcf_index_load(fn);
	bgt->b0 = bcf_init1();

	free(fn);
	bgt_set_samples(bgt, bgt->n_samples, bgt->samples);
	return bgt;
}

void bgt_close(bgt_t *bgt)
{
	int i;
	if (bgt->b0) bcf_destroy1(bgt->b0);
	free(bgt->sub);
	if (bgt->h_sub) bcf_hdr_destroy(bgt->h_sub);
	hts_itr_destroy(bgt->itr);
	hts_idx_destroy(bgt->idx);
	bcf_hdr_destroy(bgt->h0);
	vcf_close(bgt->bcf);
	pbf_close(bgt->pb);
	kh_destroy(s2i, (kh_s2i_t*)bgt->h_samples);
	for (i = 0; i < bgt->n_samples; ++i)
		free(bgt->samples[i]);
	free(bgt->samples);
	free(bgt);
}

void bgt_set_samples(bgt_t *bgt, int n, char *const* samples)
{
	int i, last, *t;
	const khash_t(s2i) *h = (khash_t(s2i)*)bgt->h_samples;
	kstring_t str = {0,0,0};

	bgt->sub = (int*)realloc(bgt->sub, n * sizeof(int));
	for (i = 0; i < n; ++i) {
		khint_t k;
		k = kh_get(s2i, h, samples[i]);
		bgt->sub[i] = k != kh_end(h)? i : bgt->n_samples;
	}
	radix_sort_i(bgt->sub, bgt->sub + n);
	for (i = bgt->n_sub = 1, last = bgt->sub[0]; i < n; ++i) // remove unidentified or duplicated samples
		if (bgt->sub[i] != bgt->n_samples && bgt->sub[i] != last)
			bgt->sub[bgt->n_sub++] = bgt->sub[i], last = bgt->sub[i];

	if (bgt->h_sub) bcf_hdr_destroy(bgt->h_sub);
	bgt->h_sub = bcf_hdr_init();
	kputsn(bgt->h0->text, bgt->h0->l_text, &str);
	if (str.s[str.l-1] == 0) --str.l;
	kputs("\tFORMAT", &str);
	for (i = 0; i < bgt->n_sub; ++i) {
		kputc('\t', &str);
		kputs(bgt->samples[bgt->sub[i]], &str);
	}
	bgt->h_sub->text = str.s;
	bgt->h_sub->l_text = str.l + 1; // including the last NULL
	bcf_hdr_parse(bgt->h_sub);

	t = (int*)malloc(bgt->n_sub * 2 * sizeof(int));
	for (i = 0; i < bgt->n_sub; ++i)
		t[i<<1|0] = bgt->sub[i]<<1|0, t[i<<1|1] = bgt->sub[i]<<1|1;
	pbf_subset(bgt->pb, bgt->n_sub<<1, t);
	free(t);
}

int bgt_set_region(bgt_t *bgt, const char *reg)
{
	bgt->itr = bcf_itr_querys(bgt->idx, bgt->h0, reg);
	return bgt->itr? 0 : -1;
}

void bcf_copy(bcf1_t *dst, const bcf1_t *src)
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

int bgt_bits2gt[4] = { (0+1)<<1, (1+1)<<1, 0<<1, (2+1)<<1 };

int bgt_read_core(bgt_t *bgt)
{
	int i, id, row;
	row = bgt->itr? bcf_itr_next((BGZF*)bgt->bcf->fp, bgt->itr, bgt->b0) : vcf_read1(bgt->bcf, bgt->h0, bgt->b0);
	if (row < 0) return row;
	assert(bgt->b0->n_sample == 0); // there shouldn't be any sample fields
	row = -1;
	id = bcf_id2int(bgt->h0, BCF_DT_ID, "_row");
	assert(id > 0);
	bcf_unpack(bgt->b0, BCF_UN_INFO);
	for (i = 0; i < bgt->b0->n_info; ++i) {
		bcf_info_t *p = &bgt->b0->d.info[i];
		if (p->key == id) row = p->v1.i;
	}
	assert(row >= 0);
	return row;
}

void bgt_gen_gt(const bcf_hdr_t *h, bcf1_t *b, int m, const uint8_t **a)
{
	int id, i;
	id = bcf_id2int(h, BCF_DT_ID, "GT");
	b->n_fmt = 1; b->n_sample = m;
	bcf_enc_int1(&b->indiv, id);
	bcf_enc_size(&b->indiv, 2, BCF_BT_INT8);
	ks_resize(&b->indiv, b->indiv.l + b->n_sample*2 + 1);
	for (i = 0; i < b->n_sample<<1; ++i)
		b->indiv.s[b->indiv.l++] = bgt_bits2gt[a[1][i]<<1 | a[0][i]];
	b->indiv.s[b->indiv.l] = 0;
}

int bgt_read(bgt_t *bgt, bcf1_t *b)
{
	int ret;
	const uint8_t **a;
	ret = bgt_read_core(bgt);
	if (ret < 0) return ret;
	pbf_seek(bgt->pb, ret);
	a = pbf_read(bgt->pb);
	bcf_copy(b, bgt->b0);
	bgt_gen_gt(bgt->h_sub, b, bgt->n_sub, a);
	return ret;
}
/*
static inline void bgt_copy_raw(int m, bgt_raw1_t *dst, const bgt_raw1_t *src)
{
	if (!dst->copied) {
		dst->bit[0] = (uint8_t*)malloc(m);
		dst->bit[1] = (uint8_t*)malloc(m);
		dst->copied = 1;
	}
	memcpy(dst->bit[0], src->bit[0], m);
	memcpy(dst->bit[1], src->bit[1], m);
	dst->row = src->row;
	bcf_copy(dst->b, src->b);
}

int bgt_rawpos_read(bgt_t *bgt, bgt_rawpos_t *p)
{
	p->n_b = 0;
	if (p->next.row < 0) return -1; // end-of-file
	if (p->next.b == 0 && bgt_read_raw(bgt, &p->next, 0) < 0) // this is the first call
		return -2;
	p->rid = p->next.b->rid, p->pos = p->next.b->pos;
	do {
		if (p->n_b == p->m_b) {
			int i, oldm = p->m_b;
			p->m_b = p->m_b? p->m_b<<1 : 4;
			p->b = (bgt_raw1_t*)realloc(p->b, p->m_b * sizeof(bgt_raw1_t));
			memset(&p->b[oldm], 0, (p->m_b - oldm) * sizeof(bgt_raw1_t));
			for (i = oldm; i < p->m_b; ++i) p->b[i].b = bcf_init1();
		}
		bgt_copy_raw(bgt->n_sub, &p->b[p->n_b++], &p->next);
		if (bgt_read_raw(bgt, &p->next, 0) < 0)
			p->next.row = -1;
	} while (p->next.row >= 0 && p->next.b->rid == p->rid && p->next.b->pos == p->pos);
	return p->n_b;
}

bgtm_t *bgtm_open(int n_files, char **fns)
{
	bgtm_t *bm;
	int i;
	bm = (bgtm_t*)calloc(1, sizeof(bgtm_t));
	bm->n_bgt = n_files;
	bm->bgt = (bgt_t**)calloc(bm->n_bgt, sizeof(void*));
	for (i = 0; i < bm->n_bgt; ++i)
		bm->bgt[i] = bgt_open(fns[i]);
	bm->rp = (bgt_rawpos_t*)calloc(bm->n_bgt, sizeof(bgt_rawpos_t));
	return bm;
}

void bgtm_close(bgtm_t *bm)
{
	int i;
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_close(bm->bgt[i]);
	free(bm->bgt);
	free(bm);
}

int bgtm_read_raw(bgtm_t *bm)
{
	int i, n_rest = 0;
	for (i = 0; i < bm->n_bgt; ++i) {
		if (bm->rp[i].n_b == 0)
			bgt_rawpos_read(bm->bgt[i], &bm->rp[i]);
		n_rest += bm->rp[i].n_b;
	}
	if (n_rest == 0) return -1;
	// search for the smallest allele
	return 0;
}
*/
