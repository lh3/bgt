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

/**********************
 * Single BGT reading *
 **********************/

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
	free(bgt->out);
	if (bgt->h_out) bcf_hdr_destroy(bgt->h_out);
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

	bgt->out = (int*)realloc(bgt->out, n * sizeof(int));
	for (i = 0; i < n; ++i) {
		khint_t k;
		k = kh_get(s2i, h, samples[i]);
		bgt->out[i] = k != kh_end(h)? kh_val(h, k) : bgt->n_samples;
	}
	radix_sort_i(bgt->out, bgt->out + n);
	for (i = bgt->n_out = 1, last = bgt->out[0]; i < n; ++i) // remove unidentified or duplicated samples
		if (bgt->out[i] != bgt->n_samples && bgt->out[i] != last)
			bgt->out[bgt->n_out++] = bgt->out[i], last = bgt->out[i];
	if (bgt->n_out == 0) return;

	if (bgt->h_out) bcf_hdr_destroy(bgt->h_out);
	bgt->h_out = bcf_hdr_init();
	kputsn(bgt->h0->text, bgt->h0->l_text, &str);
	if (str.s[str.l-1] == 0) --str.l;
	kputs("\tFORMAT", &str);
	for (i = 0; i < bgt->n_out; ++i) {
		kputc('\t', &str);
		kputs(bgt->samples[bgt->out[i]], &str);
	}
	bgt->h_out->text = str.s;
	bgt->h_out->l_text = str.l + 1; // including the last NULL
	bcf_hdr_parse(bgt->h_out);

	t = (int*)malloc(bgt->n_out * 2 * sizeof(int));
	for (i = 0; i < bgt->n_out; ++i)
		t[i<<1|0] = bgt->out[i]<<1|0, t[i<<1|1] = bgt->out[i]<<1|1;
	pbf_subset(bgt->pb, bgt->n_out<<1, t);
	free(t);

	bgt->b0->shared.l = 0; // mark b0 unread
}

int bgt_set_region(bgt_t *bgt, const char *reg)
{
	bgt->itr = bcf_itr_querys(bgt->idx, bgt->h0, reg);
	bgt->b0->shared.l = 0; // mark b0 unread
	return bgt->itr? 0 : -1;
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
	bcfcpy(b, bgt->b0);
	bgt_gen_gt(bgt->h_out, b, bgt->n_out, a);
	return ret;
}

/*************************************
 * reading records with the same pos *
 *************************************/

static void append_to_pos(bgt_pos_t *p, const bcf1_t *b0, int m, const uint8_t **a)
{
	bgt_rec_t *r;
	if (p->n_b == p->m_b) {
		int oldm = p->m_b;
		p->m_b = p->m_b? p->m_b<<1 : 4;
		p->b = (bgt_rec_t*)realloc(p->b, p->m_b * sizeof(bgt_rec_t));
		memset(&p->b[oldm], 0, (p->m_b - oldm) * sizeof(bgt_rec_t));
	}
	r = &p->b[p->n_b++];
	if (r->b0 == 0) r->b0 = bcf_init1();
	bcfcpy(r->b0, b0);
	r->a[0] = (uint8_t*)realloc(r->a[0], m);
	r->a[1] = (uint8_t*)realloc(r->a[1], m);
	memcpy(r->a[0], a[0], m);
	memcpy(r->a[1], a[1], m);
}

int bgt_read_pos(bgt_t *bgt, bgt_pos_t *p)
{
	p->n_b = 0;
	if (p->finished || bgt->n_out == 0) return -1; // end-of-file or nothing to read
	if (bgt->b0->shared.l == 0 && (p->row = bgt_read_core(bgt)) < 0) {
		p->finished = 1;
		return p->row;
	}
	p->rid = bgt->b0->rid, p->pos = bgt->b0->pos;
	do {
		const uint8_t **a;
		pbf_seek(bgt->pb, p->row);
		a = pbf_read(bgt->pb);
		append_to_pos(p, bgt->b0, bgt->n_out<<1, a);
		p->row = bgt_read_core(bgt); // read the next b0
	} while (p->row >= 0 && bgt->b0->rid == p->rid && bgt->b0->pos == p->pos);
	if (p->row < 0) p->finished = 1;
	return p->n_b;
}

/*********************
 * Multi BGT reading *
 *********************/

bgtm_t *bgtm_open(int n_files, char *const*fns)
{
	bgtm_t *bm;
	int i, j, k, n_samples = 0;
	char **samples;
	bm = (bgtm_t*)calloc(1, sizeof(bgtm_t));
	bm->n_bgt = n_files;
	bm->bgt = (bgt_t**)calloc(bm->n_bgt, sizeof(void*));
	for (i = 0; i < bm->n_bgt; ++i) {
		bm->bgt[i] = bgt_open(fns[i]);
		n_samples += bm->bgt[i]->n_samples;
	}
	samples = (char**)malloc(n_samples * sizeof(char*));
	for (i = k = 0; i < bm->n_bgt; ++i) {
		bgt_t *bgt = bm->bgt[i];
		for (j = 0; j < bgt->n_samples; ++j)
			samples[k++] = bgt->samples[j];
	}
	bgtm_set_samples(bm, n_samples, samples);
	free(samples);
	bm->p = (bgt_pos_t*)calloc(bm->n_bgt, sizeof(bgt_pos_t));
	return bm;
}

void bgtm_close(bgtm_t *bm)
{
	int i, j;
	bcf_hdr_destroy(bm->h_out);
	free(bm->a[0]); free(bm->a[1]);
	for (i = 0; i < bm->n_bgt; ++i) {
		bgt_pos_t *p = &bm->p[i];
		for (j = 0; j < p->m_b; ++j) {
			if (p->b[j].b0 == 0) continue;
			bcf_destroy1(p->b[j].b0);
			free(p->b[j].a[0]); free(p->b[j].a[1]);
		}
		free(p->b);
		bgt_close(bm->bgt[i]);
	}
	free(bm->p);
	free(bm->bgt);
	free(bm);
}

void bgtm_set_samples(bgtm_t *bm, int n, char *const* samples)
{
	int i, j;
	kstring_t h = {0,0,0};
	bcf_hdr_t *h0;
	if (bm->n_bgt == 0) return;
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_set_samples(bm->bgt[i], n, samples);

	h0 = bm->bgt[0]->h0; // FIXME: test if headers are consistent
	kputs("##fileformat=VCFv4.1\n", &h);
	kputs("##INFO=<ID=AC,Number=A,Type=String,Description=\"Count of alternate alleles\">\n", &h);
	kputs("##INFO=<ID=AN,Number=A,Type=String,Description=\"Count of total alleles\">\n", &h);
	kputs("##INFO=<ID=END,Number=1,Type=Integer,Description=\"Ending position\">\n", &h);
	kputs("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", &h);
	kputs("##ALT=<ID=M,Description=\"Multi-allele\">\n", &h);
	kputs("##ALT=<ID=DEL,Description=\"Deletion\">\n", &h);
	kputs("##ALT=<ID=DUP,Description=\"Duplication\">\n", &h);
	kputs("##ALT=<ID=INS,Description=\"Insertion\">\n", &h);
	kputs("##ALT=<ID=INV,Description=\"Inversion\">\n", &h);
	kputs("##ALT=<ID=DUP:TANDEM,Description=\"Tandem duplication\">\n", &h);
	kputs("##ALT=<ID=DEL:ME,Description=\"Deletion of mobile element\">\n", &h);
	kputs("##ALT=<ID=INS:ME,Description=\"Insertion of mobile element\">\n", &h);
	for (i = 0; i < h0->n[BCF_DT_CTG]; ++i)
		ksprintf(&h, "##contig=<ID=%s,length=%d>\n", h0->id[BCF_DT_CTG][i].key, h0->id[BCF_DT_CTG][i].val->info[0]);
	kputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", &h);
	for (i = bm->n_out = 0; i < bm->n_bgt; ++i) {
		bgt_t *bgt = bm->bgt[i];
		bm->n_out += bgt->n_out;
		for (j = 0; j < bgt->n_out; ++j) {
			kputc('\t', &h);
			kputs(bgt->samples[bgt->out[j]], &h);
		}
	}
	if (bm->h_out) bcf_hdr_destroy(bm->h_out);
	bm->h_out = bcf_hdr_init();
	bm->h_out->l_text = h.l + 1, bm->h_out->m_text = h.m, bm->h_out->text = h.s;
	bcf_hdr_parse(bm->h_out);
	bm->a[0] = (uint8_t*)realloc(bm->a[0], bm->n_out<<1);
	bm->a[1] = (uint8_t*)realloc(bm->a[1], bm->n_out<<1);
}

int bgtm_set_region(bgtm_t *bm, const char *reg)
{
	int i;
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_set_region(bm->bgt[i], reg);
	return 0;
}

int bgtm_read(bgtm_t *bm, bcf1_t *b)
{
	int i, j, off = 0, n_rest = 0, max_allele = 0, l_ref;
	bcf1_t *b0 = 0;
	// fill the buffer
	for (i = n_rest = 0; i < bm->n_bgt; ++i) {
		if (bm->p[i].n_b == 0)
			bgt_read_pos(bm->bgt[i], &bm->p[i]);
		n_rest += bm->p[i].n_b;
	}
	if (n_rest == 0) return -1;
	// search for the smallest allele
	for (i = 0; i < bm->n_bgt; ++i) {
		bgt_pos_t *p = &bm->p[i];
		for (j = 0; j < p->n_b; ++j) {
			if (b0) {
				int c;
				c = bcfcmp(b0, p->b[j].b0);
				if (c > 0) b0 = p->b[j].b0, max_allele = b0->n_allele;
				else if (c == 0)
					max_allele = p->b[j].b0->n_allele > max_allele? p->b[j].b0->n_allele : max_allele;
			} else b0 = p->b[j].b0, max_allele = b0->n_allele;
		}
	}
	assert(b0 && max_allele >= 2);
	l_ref = bcfcpy_min(b, b0, max_allele > 2? "<M>" : 0);
	if (l_ref != b->rlen)
		bcf_append_info_int(bm->h_out, b, "END", b->pos + b->rlen);
	// generate bm->a
	for (i = 0; i < bm->n_bgt; ++i) {
		bgt_pos_t *p = &bm->p[i];
		bgt_t *bgt = bm->bgt[i];
		bgt_rec_t q, swap;
		int n_min = 0, end = p->n_b - 1;
		if (bgt->n_out == 0) continue;
		memset(&q, 0, sizeof(bgt_rec_t)); // suppress a gcc warning
		for (j = 0; j <= end; ++j) {
			if (bcfcmp(b, p->b[j].b0) == 0) { // then swap to the end
				q = p->b[j];
				++n_min;
				if (j != end) swap = p->b[end], p->b[end] = p->b[j], p->b[j] = swap;
				--end, --j;
			}
		}
		p->n_b = end + 1;
		if (n_min) { // copy; when there are multiple, we only copy the last one
			memcpy(bm->a[0] + off, q.a[0], bgt->n_out<<1);
			memcpy(bm->a[1] + off, q.a[1], bgt->n_out<<1);
		} else { // add missing values
			memset(bm->a[0] + off, 0, bgt->n_out<<1);
			memset(bm->a[1] + off, 1, bgt->n_out<<1);
		}
		off += bgt->n_out<<1;
	}
	bgt_gen_gt(bm->h_out, b, bm->n_out, (const uint8_t**)bm->a);
	return 0;
}
