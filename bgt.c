#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "bgt.h"
#include "kstring.h"

#include "khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)

void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *_h);

/**********************
 * Single BGT reading *
 **********************/

bgt_t *bgt_open(const char *prefix)
{
	char *fn;
	bgt_t *bgt;

	bgt = (bgt_t*)calloc(1, sizeof(bgt_t));
	fn = (char*)malloc(strlen(prefix) + 9);

	// read samples
	sprintf(fn, "%s.spl", prefix);
	bgt->samples = hts_readlines(fn, &bgt->n_samples);

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
	for (i = 0; i < bgt->n_samples; ++i)
		free(bgt->samples[i]);
	free(bgt->samples);
	free(bgt);
}

void bgt_set_samples(bgt_t *bgt, int n, char *const* samples)
{
	int i, *t, tmp;
	khash_t(s2i) *h;
	kstring_t str = {0,0,0};

	// build hash table for samples
	h = kh_init(s2i);
	for (i = 0; i < n; ++i) {
		khint_t k;
		int absent;
		k = kh_put(s2i, h, samples[i], &absent);
	}
	tmp = kh_size(h) < bgt->n_samples? kh_size(h) : bgt->n_samples;

	// fill bgt->out
	bgt->n_out = 0;
	bgt->out = (int*)realloc(bgt->out, tmp * sizeof(int));
	for (i = 0; i < bgt->n_samples; ++i)
		if (kh_get(s2i, h, bgt->samples[i]) != kh_end(h))
			bgt->out[bgt->n_out++] = i;
	kh_destroy(s2i, h);

	// build ->h_out VCF header
	if (bgt->h_out) bcf_hdr_destroy(bgt->h_out);
	bgt->h_out = bcf_hdr_init();
	kputsn(bgt->h0->text, bgt->h0->l_text, &str);
	if (str.s[str.l-1] == 0) --str.l;
	if (bgt->n_out > 0) {
		kputs("\tFORMAT", &str);
		for (i = 0; i < bgt->n_out; ++i) {
			kputc('\t', &str);
			kputs(bgt->samples[bgt->out[i]], &str);
		}
	}
	bgt->h_out->text = str.s;
	bgt->h_out->l_text = str.l + 1; // including the last NULL
	bcf_hdr_parse(bgt->h_out);

	// subsetting pbf
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

int bgt_read_core0(bgt_t *bgt)
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
	b->indiv.l = 0;
	bcf_enc_int1(&b->indiv, id);
	bcf_enc_size(&b->indiv, 2, BCF_BT_INT8);
	ks_resize(&b->indiv, b->indiv.l + b->n_sample*2 + 1);
	for (i = 0; i < b->n_sample<<1; ++i)
		b->indiv.s[b->indiv.l++] = bgt_bits2gt[a[1][i]<<1 | a[0][i]];
	b->indiv.s[b->indiv.l] = 0;
}

int bgt_read_core(bgt_t *bgt)
{
	if (bgt->bed) {
		int ret;
		while ((ret = bgt_read_core0(bgt)) >= 0) {
			int r;
			r = bed_overlap(bgt->bed, bgt->h_out->id[BCF_DT_CTG][bgt->b0->rid].key, bgt->b0->pos, bgt->b0->pos + bgt->b0->rlen);
			if (bgt->bed_excl && r) continue;
			if (!bgt->bed_excl && !r) continue;
			break;
		}
		return ret;
	} else return bgt_read_core0(bgt);
}

int bgt_read_rec(bgt_t *bgt, bgt_rec_t *r)
{
	int row;
	const uint8_t **a;
	r->b0 = 0, r->a[0] = r->a[1] = 0;
	if (bgt->n_out == 0) return -1;
	if ((row = bgt_read_core(bgt)) < 0) return row;
	r->b0 = bgt->b0;
	pbf_seek(bgt->pb, row);
	a = pbf_read(bgt->pb);
	r->a[0] = (uint8_t*)a[0], r->a[1] = (uint8_t*)a[1];
	return row;
}

int bgt_read(bgt_t *bgt, bcf1_t *b)
{
	int ret;
	bgt_rec_t r;
	if ((ret = bgt_read_rec(bgt, &r)) < 0) return ret;
	bcfcpy(b, r.b0);
	bgt_gen_gt(bgt->h_out, b, bgt->n_out, r.a);
	return ret;
}

void bgt_set_bed(bgt_t *bgt, const void *bed, int excl) { bgt->bed = bed, bgt->bed_excl = excl; }

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
	bm->r = (bgt_rec_t*)calloc(bm->n_bgt, sizeof(bgt_rec_t));
	return bm;
}

void bgtm_close(bgtm_t *bm)
{
	int i;
	free(bm->group);
	free(bm->sample_idx);
	bcf_hdr_destroy(bm->h_out);
	free(bm->a[0]); free(bm->a[1]);
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_close(bm->bgt[i]);
	free(bm->r); free(bm->bgt); free(bm);
}

void bgtm_set_samples(bgtm_t *bm, int n, char *const* samples)
{
	int i, j, m;
	kstring_t h = {0,0,0};
	bcf_hdr_t *h0;

	if (bm->n_bgt == 0) return;
	for (i = bm->n_out = 0; i < bm->n_bgt; ++i) {
		bgt_set_samples(bm->bgt[i], n, samples);
		bm->n_out += bm->bgt[i]->n_out;
	}
	bm->group = (uint8_t*)realloc(bm->group, bm->n_out);
	memset(bm->group, 0, bm->n_out);
	bm->sample_idx = (uint64_t*)realloc(bm->sample_idx, bm->n_out * 8);
	for (i = m = 0; i < bm->n_bgt; ++i)
		for (j = 0; j < bm->bgt[i]->n_out; ++j)
			bm->sample_idx[m++] = (uint64_t)i<<32 | bm->bgt[i]->out[j];

	h0 = bm->bgt[0]->h0; // FIXME: test if headers are consistent
	kputs("##fileformat=VCFv4.1\n", &h);
	kputs("##INFO=<ID=AC,Number=A,Type=String,Description=\"Count of alternate alleles\">\n", &h);
	kputs("##INFO=<ID=AN,Number=A,Type=String,Description=\"Count of total alleles\">\n", &h);
	for (i = 1; i <= BGT_MAX_GROUPS; ++i) {
		ksprintf(&h, "##INFO=<ID=AC%d,Number=A,Type=String,Description=\"Count of alternate alleles for sample group %d\">\n", i, i);
		ksprintf(&h, "##INFO=<ID=AN%d,Number=A,Type=String,Description=\"Count of total alleles for sample group %d\">\n", i, i);
	}
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
	kputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", &h);
	if (!(bm->flag & BGT_F_NO_GT)) {
		kputs("\tFORMAT", &h);
		for (i = 0; i < bm->n_bgt; ++i) {
			bgt_t *bgt = bm->bgt[i];
			for (j = 0; j < bgt->n_out; ++j) {
				kputc('\t', &h);
				kputs(bgt->samples[bgt->out[j]], &h);
			}
		}
	}
	if (bm->h_out) bcf_hdr_destroy(bm->h_out);
	bm->h_out = bcf_hdr_init();
	bm->h_out->l_text = h.l + 1, bm->h_out->m_text = h.m, bm->h_out->text = h.s;
	bcf_hdr_parse(bm->h_out);
	bm->a[0] = (uint8_t*)realloc(bm->a[0], bm->n_out<<1);
	bm->a[1] = (uint8_t*)realloc(bm->a[1], bm->n_out<<1);
}

void bgtm_add_group(bgtm_t *bm, int n, char *const* samples)
{
	int i, j, k, absent;
	khash_t(s2i) *h;
	h = kh_init(s2i);
	for (i = 0; i < n; ++i)
		kh_put(s2i, h, samples[i], &absent);
	for (i = k = 0; i < bm->n_bgt; ++i) {
		bgt_t *bgt = bm->bgt[i];
		for (j = 0; j < bgt->n_out; ++j, ++k)
			if (kh_get(s2i, h, bgt->samples[bgt->out[j]]) != kh_end(h))
				bm->group[k] |= 1<<bm->n_groups;

	}
	kh_destroy(s2i, h);
	++bm->n_groups;
}

int bgtm_set_region(bgtm_t *bm, const char *reg)
{
	int i;
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_set_region(bm->bgt[i], reg);
	return 0;
}

int bgtm_read_core(bgtm_t *bm, bcf1_t *b)
{
	int i, j, off = 0, n_rest = 0, max_allele = 0, l_ref;
	const bcf1_t *b0 = 0;

	// fill the buffer
	for (i = n_rest = 0; i < bm->n_bgt; ++i) {
		if (bm->r[i].b0 == 0)
			bgt_read_rec(bm->bgt[i], &bm->r[i]);
		n_rest += (bm->r[i].b0 != 0);
	}
	if (n_rest == 0) return -1;
	// search for the smallest allele
	for (i = 0; i < bm->n_bgt; ++i) {
		bgt_rec_t *r = &bm->r[i];
		if (r->b0 == 0) continue;
		if (b0) {
			j = bcfcmp(b0, r->b0);
			if (j > 0) b0 = r->b0, max_allele = b0->n_allele;
			else if (j == 0)
				max_allele = r->b0->n_allele > max_allele? r->b0->n_allele : max_allele;
		} else b0 = r->b0, max_allele = b0->n_allele;
	}
	assert(b0 && max_allele >= 2);
	l_ref = bcfcpy_min(b, b0, max_allele > 2? "<M>" : 0);
	if (l_ref != b->rlen) {
		int32_t val = b->pos + b->rlen;
		bcf_append_info_ints(bm->h_out, b, "END", 1, &val);
	}
	// generate bm->a
	for (i = 0; i < bm->n_bgt; ++i) {
		bgt_rec_t *r = &bm->r[i];
		bgt_t *bgt = bm->bgt[i];
		if (bgt->n_out == 0) continue;
		if (r->b0 && bcfcmp(b, r->b0) == 0) { // copy
			r->b0 = 0;
			memcpy(bm->a[0] + off, r->a[0], bgt->n_out<<1);
			memcpy(bm->a[1] + off, r->a[1], bgt->n_out<<1);
		} else { // add missing values
			memset(bm->a[0] + off, 0, bgt->n_out<<1);
			memset(bm->a[1] + off, 1, bgt->n_out<<1);
		}
		off += bgt->n_out<<1;
	}
	assert(off == bm->n_out<<1);
	if ((bm->flag & BGT_F_SET_AC) || bm->filter_func) {
		int32_t an, ac[2], cnt[4], gac1[BGT_MAX_GROUPS+1], gan[BGT_MAX_GROUPS+1];
		memset(cnt, 0, 4 * 4);
		for (i = 0; i < bm->n_out<<1; ++i)
			++cnt[bm->a[1][i]<<1 | bm->a[0][i]];
		an = cnt[0] + cnt[1] + cnt[3];
		ac[0] = cnt[1], ac[1] = cnt[3];
		bcf_append_info_ints(bm->h_out, b, "AN", 1, &an);
		bcf_append_info_ints(bm->h_out, b, "AC", b->n_allele - 1, ac);
		if (bm->n_groups > 0) {
			int32_t gcnt[BGT_MAX_GROUPS+1][4];
			memset(gcnt, 0, 4 * (BGT_MAX_GROUPS+1) * 4);
			// the following two blocks achieve the same goal. The 1st is faster if there are not many samples
			if (bm->n_out<<1 < 1024) {
				int32_t j;
				for (i = 0; i < bm->n_out<<1; ++i) {
					int ht = bm->a[1][i]<<1 | bm->a[0][i];
					if (bm->group[i>>1])
						for (j = 0; j < bm->n_groups; ++j)
							if (bm->group[i>>1] & 1<<j) ++gcnt[j+1][ht];
				}
			} else {
				int32_t j, k, gcnt256[256][4];
				memset(gcnt256, 0, 256 * 4 * 4);
				for (i = 0; i < bm->n_out<<1; ++i)
					++gcnt256[bm->group[i>>1]][bm->a[1][i]<<1 | bm->a[0][i]];
				for (i = 0; i < 256; ++i)
					for (j = 0; j < bm->n_groups; ++j)
						if (i & 1<<j)
							for (k = 0; k < 4; ++k)
								gcnt[j+1][k] += gcnt256[i][k];
			}
			for (i = 1; i <= bm->n_groups; ++i) {
				char key[4];
				int32_t gac[2];
				key[3] = 0;
				gan[i] = gcnt[i][0] + gcnt[i][1] + gcnt[i][3];
				gac[0] = gac1[i] = gcnt[i][1];
				gac[1] = gcnt[i][3];
				key[0] = 'A', key[1] = 'N', key[2] = '0' + i; bcf_append_info_ints(bm->h_out, b, key, 1, &gan[i]);
				key[0] = 'A', key[1] = 'C', key[2] = '0' + i; bcf_append_info_ints(bm->h_out, b, key, b->n_allele - 1, gac);
			}
		}
		if (bm->filter_func && bm->filter_func(bm->h_out, b, an, ac[0], bm->n_groups, gan, gac1, bm->filter_data))
			return 1;
	}
	if ((bm->flag & BGT_F_NO_GT) == 0)
		bgt_gen_gt(bm->h_out, b, bm->n_out, (const uint8_t**)bm->a);
	return 0;
}

int bgtm_read(bgtm_t *bm, bcf1_t *b)
{
	int ret;
	while ((ret = bgtm_read_core(bm, b)) > 0);
	return ret;
}

void bgtm_set_bed(bgtm_t *bm, const void *bed, int excl)
{
	int i;
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_set_bed(bm->bgt[i], bed, excl);
}

void bgtm_set_flag(bgtm_t *bm, int flag) { bm->flag = flag; }
void bgtm_set_filter(bgtm_t *bm, bgt_filter_f func, void *data) { bm->filter_func = func; bm->filter_data = data; }

/******************
 * Allele parsing *
 ******************/

int bgt_al_parse(const char *al, bgt_allele_t *a)
{
	char *p = (char*)al, *ref = 0, *alt = 0;
	int sep = ':', off, tmp;
	a->chr.l = 0; a->alt = 0; a->pos = -1; a->rlen = -1;
	for (; *p && *p != sep; ++p);
	if (*p == 0) return -1;
	kputsn(al, p - al, &a->chr); kputc(0, &a->chr);
	++p; // skip the delimiter
	if (!isdigit(*p)) return -1;
	a->pos = strtol(p, &p, 10) - 1; // read position
	if (*p != sep) return -1;
	++p; // skip the delimiter
	if (isdigit(*p)) {
		a->rlen = strtol(p, &p, 10);
	} else if (isalpha(*p)) {
		ref = p;
		for (; isalpha(*p); ++p);
		a->rlen = p - ref;
	} else return -1;
	if (*p != sep) return -1;
	alt = ++p;
	for (off = 0; *p && isalpha(*p); ++p)
		if (ref && toupper(*p) == toupper(ref[off])) ++off;
		else break;
	a->pos += off; a->rlen -= off;
	tmp = a->chr.l;
	kputs(alt + off, &a->chr);
	a->alt = a->chr.s + tmp;
	if (ref) {
		int l_alt = a->chr.s + a->chr.l - a->alt;
		int min_l = l_alt < a->rlen? l_alt : a->rlen;
		ref += off;
		for (off = 0; off < min_l && isalpha(ref[a->rlen - 1 - off]) && toupper(ref[a->rlen - 1 - off]) == toupper(a->alt[l_alt - 1 - off]); ++off);
		a->rlen -= off;
		a->alt[l_alt - off] = 0;
		a->chr.l -= off;
	}
	return 0;
}
