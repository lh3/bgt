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

	free(fn);
	return bgt;
}

void bgt_close(bgt_t *bgt)
{
	int i;
	free(bgt->sub);
	bcf_hdr_destroy(bgt->h_sub);
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

void bgt_set_sub(bgt_t *bgt, int n, char *const* samples)
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
	kputsn(bgt->h0->text, bgt->h0->l_text, &str);

	if (bgt->h_sub) bcf_hdr_destroy(bgt->h_sub);
	bgt->h_sub = bcf_hdr_init();
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
