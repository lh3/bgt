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
	free(bgt->sub);
	if (bgt->h_sub) bcf_hdr_destroy(bgt->h_sub);
	hts_itr_destroy(bgt->itr);
	hts_idx_destroy(bgt->idx);
	bcf_hdr_destroy(bgt->h0);
	vcf_close(bgt->bcf);
	pbf_close(bgt->pb);
	kh_destroy(s2i, bgt->h_samples);
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
	kputs("\tFORMAT", &str);
	for (i = 0; i < bgt->n_sub; ++i) {
		kputc('\t', &str);
		kputs(bgt->samples[bgt->sub[i]], &str);
	}
	bgt->h_sub->text = str.s;
	bgt->h_sub->l_text = str.l;
	bcf_hdr_parse(bgt->h_sub);

	t = malloc(bgt->n_sub * 2 * sizeof(int));
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

int bgt_read(bgt_t *bgt, bcf1_t *b)
{
	int ret, i, id, row = -1;
	const uint8_t **a;
	int32_t *gt;

	ret = bgt->itr? bcf_itr_next((BGZF*)bgt->bcf->fp, bgt->itr, bgt->b0) : vcf_read1(bgt->bcf, bgt->h0, bgt->b0);
	if (ret < 0) return ret;
	assert(bgt->b0->n_sample == 0); // there shouldn't be any sample fields

	id = bcf_id2int(bgt->h0, BCF_DT_ID, "_row");
	assert(id > 0);
	bcf_unpack(bgt->b0, BCF_UN_INFO);
	for (i = 0; i < bgt->b0->n_info; ++i) {
		bcf_info_t *p = &bgt->b0->d.info[i];
		if (p->key == id) row = p->v1.i;
	}
	assert(row >= 0);
	bcf_copy(b, bgt->b0);

	id = bcf_id2int(bgt->h_sub, BCF_DT_ID, "GT");
	b->n_fmt = 1; b->n_sample = bgt->n_sub;
	bcf_enc_int1(&b->indiv, id);
	pbf_seek(bgt->pb, row);
	a = pbf_read(bgt->pb);
	gt = calloc(b->n_sample * 2, 4);
	for (i = 0; i < b->n_sample<<1; ++i) {
		int g = a[1][i]<<1 | a[0][i];
		if (g == 2) gt[i] = 0;
		else if (g == 3) gt[i] = (2+1)<<1;
		else gt[i] = (g+1)<<1;
	}
	bcf_enc_vint(&b->indiv, 2 * b->n_sample, gt, 2);
	free(gt);
	return 0;
}
