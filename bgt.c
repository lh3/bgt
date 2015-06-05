#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <ctype.h>
#include "bgt.h"
#include "kstring.h"

#include "khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)

void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *_h);

/************
 * BGT file *
 ************/

bgt_file_t *bgt_open(const char *prefix)
{
	char *fn;
	bgt_file_t *bf;
	BGZF *fp;
	bf = (bgt_file_t*)calloc(1, sizeof(bgt_file_t));
	fn = (char*)malloc(strlen(prefix) + 9);
	bf->prefix = strdup(prefix);
	sprintf(fn, "%s.spl", prefix);
	bf->f = fmf_read(fn);
	sprintf(fn, "%s.bcf", prefix);
	fp = bgzf_open(fn, "rb");
	bf->h0 = bcf_hdr_read(fp);
	bf->idx = bcf_index_load(fn);
	bgzf_close(fp);
	free(fn);
	return bf;
}

void bgt_close(bgt_file_t *bf)
{
	if (bf == 0) return;
	hts_idx_destroy(bf->idx);
	bcf_hdr_destroy(bf->h0);
	fmf_destroy(bf->f);
	free(bf->prefix); free(bf);
}

/**********************
 * Single BGT reading *
 **********************/

/*** reader allocation/deallocation ***/

bgt_t *bgt_reader_init(const bgt_file_t *bf)
{
	char *fn;
	bgt_t *bgt;
	assert(BGT_MAX_GROUPS <= 32);
	bgt = (bgt_t*)calloc(1, sizeof(bgt_t));
	bgt->f = bf;
	fn = (char*)malloc(strlen(bf->prefix) + 9);
	sprintf(fn, "%s.pbf", bf->prefix);
	bgt->pb = pbf_open_r(fn);
	sprintf(fn, "%s.bcf", bf->prefix);
	bgt->bcf = bgzf_open(fn, "rb");
	bgt->b0 = bcf_init1();
	bcf_seekn(bgt->bcf, bgt->f->idx, 0);
	bgt->flag = (uint32_t*)calloc(bgt->f->f->n_rows, 4);
	free(fn);
	return bgt;
}

void bgt_reader_destroy(bgt_t *bgt)
{
	bcf_destroy1(bgt->b0);
	free(bgt->flag); free(bgt->group); free(bgt->out);
	if (bgt->h_out) bcf_hdr_destroy(bgt->h_out);
	hts_itr_destroy(bgt->itr);
	pbf_close(bgt->pb);
	bgzf_close(bgt->bcf);
	free(bgt);
}

/*** set samples, regions, etc. ***/

void bgt_add_group_core(bgt_t *bgt, int n, char *const* samples, const char *expr)
{
	int i;
	const fmf_t *f = bgt->f->f;

	if (n == BGT_SET_ALL_SAMPLES) {
		for (i = 0; i < f->n_rows; ++i)
			bgt->flag[i] |= 1<<bgt->n_groups;
	} else if (n > 0 || expr != 0) {
		int err, absent;
		khash_t(s2i) *h;
		khint_t k;
		kexpr_t *ke = 0;

		if (expr) {
			ke = ke_parse(expr, &err);
			if (err && ke) {
				ke_destroy(ke);
				ke = 0;
			}
		}
		h = kh_init(s2i);
		for (i = 0; i < n; ++i)
			k = kh_put(s2i, h, samples[i], &absent);
		for (i = 0; i < f->n_rows; ++i)
			if ((kh_get(s2i, h, f->rows[i].name) != kh_end(h)) || (ke && fmf_test(f, i, ke)))
				bgt->flag[i] |= 1<<bgt->n_groups;
		kh_destroy(s2i, h);
		ke_destroy(ke);
	}
	++bgt->n_groups;
}

void bgt_add_group(bgt_t *bgt, const char *expr)
{
	int is_file = 0;
	FILE *fp;
	if ((fp = fopen(expr, "r")) != 0) { // test if expr is a file
		is_file = 1;
		fclose(fp);
	}
	if (*expr == ':' || (*expr != '?' && is_file)) {
		int i, n;
		char **samples;
		samples = hts_readlines(expr, &n);
		bgt_add_group_core(bgt, n, samples, 0);
		for (i = 0; i < n; ++i) free(samples[i]);
		free(samples);
	} else bgt_add_group_core(bgt, 0, 0, expr);
}

int bgt_set_region(bgt_t *bgt, const char *reg)
{
	if (bgt->itr) bcf_itr_destroy(bgt->itr);
	bgt->itr = bcf_itr_querys(bgt->f->idx, bgt->f->h0, reg);
	bgt->b0->shared.l = 0; // mark b0 unread
	return bgt->itr? 0 : -1;
}

int bgt_set_start(bgt_t *bgt, int64_t i)
{
	return bcf_seekn(bgt->bcf, bgt->f->idx, i);
}

void bgt_set_bed(bgt_t *bgt, const void *bed, int excl) { bgt->bed = bed, bgt->bed_excl = excl; }

/*** prepare for the output ***/

void bgt_prepare(bgt_t *bgt)
{
	int i, *t;
	const fmf_t *f = bgt->f->f;
	kstring_t str = {0,0,0};

	if (bgt->n_groups == 0) bgt_add_group_core(bgt, BGT_SET_ALL_SAMPLES, 0, 0);
	for (i = 0, bgt->n_out = 0; i < f->n_rows; ++i)
		if (bgt->flag[i]) ++bgt->n_out;
	bgt->out = (int*)realloc(bgt->out, bgt->n_out * sizeof(int));
	bgt->group = (uint32_t*)realloc(bgt->group, bgt->n_out * 4);
	for (i = 0, bgt->n_out = 0; i < f->n_rows; ++i)
		if (bgt->flag[i]) bgt->out[bgt->n_out] = i, bgt->group[bgt->n_out++] = bgt->flag[i];

	// build ->h_out VCF header
	if (bgt->h_out) bcf_hdr_destroy(bgt->h_out);
	bgt->h_out = bcf_hdr_init();
	kputsn(bgt->f->h0->text, bgt->f->h0->l_text, &str);
	if (str.s[str.l-1] == 0) --str.l;
	if (bgt->n_out > 0) {
		kputs("\tFORMAT", &str);
		for (i = 0; i < bgt->n_out; ++i) {
			kputc('\t', &str);
			kputs(f->rows[bgt->out[i]].name, &str);
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

/*** read into BCF ***/

int bgt_bits2gt[4] = { (0+1)<<1, (1+1)<<1, 0<<1, (2+1)<<1 };

int bgt_read_core0(bgt_t *bgt)
{
	int i, id, row;
	row = bgt->itr? bcf_itr_next(bgt->bcf, bgt->itr, bgt->b0) : bcf_read1(bgt->bcf, bgt->b0);
	if (row < 0) return row;
	assert(bgt->b0->n_sample == 0); // there shouldn't be any sample fields
	row = -1;
	id = bcf_id2int(bgt->f->h0, BCF_DT_ID, "_row");
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
	if (bgt->h_out == 0) bgt_prepare(bgt);
	if ((ret = bgt_read_rec(bgt, &r)) < 0) return ret;
	bcfcpy(b, r.b0);
	bgt_gen_gt(bgt->h_out, b, bgt->n_out, r.a);
	return ret;
}

/*********************
 * Multi BGT reading *
 *********************/

/*** reader allocation/deallocation ***/

bgtm_t *bgtm_reader_init(int n_files, bgt_file_t *const* bf)
{
	bgtm_t *bm;
	int i;
	bm = (bgtm_t*)calloc(1, sizeof(bgtm_t));
	bm->n_bgt = n_files;
	bm->bgt = (bgt_t**)calloc(bm->n_bgt, sizeof(void*));
	for (i = 0; i < bm->n_bgt; ++i)
		bm->bgt[i] = bgt_reader_init(bf[i]);
	bm->r = (bgt_rec_t*)calloc(bm->n_bgt, sizeof(bgt_rec_t));
	return bm;
}

void bgtm_reader_destroy(bgtm_t *bm)
{
	int i;
	free(bm->hap);
	free(bm->alcnt);
	if (bm->site_flt) ke_destroy(bm->site_flt);
	free(bm->group);
	free(bm->sample_idx);
	if (bm->h_out) bcf_hdr_destroy(bm->h_out);
	free(bm->a[0]); free(bm->a[1]);
	for (i = 0; i < bm->n_al; ++i) {
		if (bm->aal) free(bm->aal[i].chr.s);
		free(bm->al[i].chr.s);
	}
	free(bm->aal);
	free(bm->al);
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_reader_destroy(bm->bgt[i]);
	free(bm->r); free(bm->bgt); free(bm);
}

/*** set samples, regions, etc. ***/

void bgtm_add_group_core(bgtm_t *bm, int n, char *const* samples, const char *expr)
{
	int i;
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_add_group_core(bm->bgt[i], n, samples, expr);
	++bm->n_groups;
}

void bgtm_add_group(bgtm_t *bm, const char *expr)
{
	int i;
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_add_group(bm->bgt[i], expr);
	++bm->n_groups;
}

int bgtm_set_region(bgtm_t *bm, const char *reg)
{
	int i;
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_set_region(bm->bgt[i], reg);
	return 0;
}

int bgtm_set_start(bgtm_t *bm, int64_t n)
{
	int i;
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_set_start(bm->bgt[i], n);
	return 0;
}

void bgtm_set_bed(bgtm_t *bm, const void *bed, int excl)
{
	int i;
	for (i = 0; i < bm->n_bgt; ++i)
		bgt_set_bed(bm->bgt[i], bed, excl);
}

void bgtm_set_flag(bgtm_t *bm, int flag) { bm->flag = flag; }

int bgtm_set_flt_site(bgtm_t *bm, const char *expr)
{
	int err;
	if (bm->site_flt) ke_destroy(bm->site_flt);
	bm->site_flt = ke_parse(expr, &err);
	if (err != 0) {
		if (bm->site_flt) ke_destroy(bm->site_flt);
		bm->site_flt = 0;
		return err;
	}
	return 0;
}

int bgtm_add_allele(bgtm_t *bm, const char *al)
{
	int ret;
	if (bm->n_al == bm->m_al) {
		int oldm = bm->m_al;
		bm->m_al = bm->m_al? bm->m_al<<1 : 4;
		bm->al = (bgt_allele_t*)realloc(bm->al, bm->m_al * sizeof(bgt_allele_t));
		memset(bm->al + oldm, 0, (bm->m_al - oldm) * sizeof(bgt_allele_t));
	}
	ret = bgt_al_parse(al, &bm->al[bm->n_al]);
	if (ret == 0) ++bm->n_al;
	return ret;
}

/*** prepare for the output ***/

void bgtm_prepare(bgtm_t *bm)
{
	int i, j, m;
	kstring_t h = {0,0,0};
	bcf_hdr_t *h0;

	if (bm->n_bgt == 0) return;

	// prepare group and sample_idx
	for (i = bm->n_out = 0; i < bm->n_bgt; ++i) {
		bgt_prepare(bm->bgt[i]);
		bm->n_out += bm->bgt[i]->n_out;
	}
	bm->group = (uint32_t*)realloc(bm->group, bm->n_out * 4);
	bm->sample_idx = (uint64_t*)realloc(bm->sample_idx, bm->n_out * 8);
	for (i = m = 0; i < bm->n_bgt; ++i) {
		for (j = 0; j < bm->bgt[i]->n_out; ++j) {
			bm->sample_idx[m] = (uint64_t)i<<32 | bm->bgt[i]->out[j];
			bm->group[m++] = bm->bgt[i]->group[j];
		}
	}

	// prepare header
	h0 = bm->bgt[0]->f->h0; // FIXME: test if headers are consistent
	kputs("##fileformat=VCFv4.1\n", &h);
	kputs("##INFO=<ID=AC,Number=A,Type=String,Description=\"Count of alternate alleles\">\n", &h);
	kputs("##INFO=<ID=AN,Number=A,Type=String,Description=\"Count of total alleles\">\n", &h);
	for (i = 1; i <= bm->n_groups; ++i) {
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
				kputs(bgt->f->f->rows[bgt->out[j]].name, &h);
			}
		}
	}
	if (bm->h_out) bcf_hdr_destroy(bm->h_out);
	bm->h_out = bcf_hdr_init();
	bm->h_out->l_text = h.l + 1, bm->h_out->m_text = h.m, bm->h_out->text = h.s;
	bcf_hdr_parse(bm->h_out);

	// prepare the haplotype arrays
	bm->a[0] = (uint8_t*)realloc(bm->a[0], bm->n_out<<1);
	bm->a[1] = (uint8_t*)realloc(bm->a[1], bm->n_out<<1);

	// set al[].rid
	for (i = 0; i < bm->n_al; ++i)
		bm->al[i].rid = bcf_name2id(bm->h_out, bm->al[i].chr.s);

	// set region if necessary
	if (bm->n_al && bm->bgt[0]->itr == 0) {
		char *reg;
		int min_pos = INT_MAX, max_pos = INT_MIN;
		for (i = 0; i < bm->n_al; ++i) {
			bgt_allele_t *a = &bm->al[i];
			min_pos = min_pos < a->pos? min_pos : a->pos;
			max_pos = max_pos > a->pos? max_pos : a->pos;
			if (a->rid < 0 || a->rid != bm->al[i].rid) break;
		}
		if (i == bm->n_al) {
			reg = (char*)alloca(strlen(bm->al[0].chr.s) + 23);
			sprintf(reg, "%s:%d-%d", bm->al[0].chr.s, min_pos+1, max_pos+1);
			bgtm_set_region(bm, reg);
		}
	}

	if (bm->n_al > 0) {
		if (bm->flag&BGT_F_CNT_AL)
			bm->alcnt = (int*)calloc(bm->n_out, sizeof(int));
		if (bm->flag&BGT_F_CNT_HAP) {
			bm->hap = (uint64_t*)calloc(bm->n_out<<1, 8);
			bm->aal = (bgt_allele_t*)calloc(bm->n_al, sizeof(bgt_allele_t));
		}
	}
}

/*** read into BCF ***/

static inline char *gen_group_key(char key[5], char nc, int g)
{
	key[0] = 'A'; key[1] = nc;
	if (g < 9) key[2] = '0' + (g+1), key[3] = 0;
	else key[2] = '0' + (g+1)/10, key[3] = '0' + (g+1)%10, key[4] = 0;
	return key;
}

int bgtm_pass_site_flt(const bgt_info_t *ss, kexpr_t *flt)
{
	int is_true, err, i;
	char key[5];
	if (flt == 0) return 1;
	ke_set_int(flt, "AN", ss->an);
	ke_set_int(flt, "AC", ss->ac[0]);
	for (i = 0; i < ss->n_groups; ++i) {
		ke_set_int(flt, gen_group_key(key, 'N', i), ss->gan[i]);
		ke_set_int(flt, gen_group_key(key, 'C', i), ss->gac[i][0]);
	}
	is_true = !!ke_eval_int(flt, &err);
	return err? 0 : is_true;
}

void bgtm_fill_info(const bcf_hdr_t *h, const bgt_info_t *ss, bcf1_t *b)
{
	bcf_append_info_ints(h, b, "AN", 1, &ss->an);
	bcf_append_info_ints(h, b, "AC", b->n_allele - 1, ss->ac);
	if (ss->n_groups > 1) {
		int i;
		char key[5];
		for (i = 0; i < ss->n_groups; ++i) {
			bcf_append_info_ints(h, b, gen_group_key(key, 'N', i), 1, &ss->gan[i]);
			bcf_append_info_ints(h, b, gen_group_key(key, 'C', i), b->n_allele - 1, ss->gac[i]);
		}
	}
}

void bgtm_cal_info(const bgtm_t *bm, bgt_info_t *ss)
{
	int32_t cnt[4], i;
	memset(cnt, 0, 4 * 4);
	ss->n_groups = bm->n_groups;
	for (i = 0; i < bm->n_out<<1; ++i)
		++cnt[bm->a[1][i]<<1 | bm->a[0][i]];
	ss->an = cnt[0] + cnt[1] + cnt[3];
	ss->ac[0] = cnt[1], ss->ac[1] = cnt[3];
	if (bm->n_groups > 1) {
		int32_t gcnt[BGT_MAX_GROUPS][4];
		memset(gcnt, 0, 4 * BGT_MAX_GROUPS * 4);
		// the following two blocks achieve the same goal. The 1st is faster if there are not many samples
		if (bm->n_out<<1 < 1024) {
			int32_t j;
			for (i = 0; i < bm->n_out<<1; ++i) {
				int ht = bm->a[1][i]<<1 | bm->a[0][i];
				if (bm->group[i>>1])
					for (j = 0; j < bm->n_groups; ++j)
						if (bm->group[i>>1] & 1<<j) ++gcnt[j][ht];
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
							gcnt[j][k] += gcnt256[i][k];
		}
		for (i = 0; i < bm->n_groups; ++i) {
			ss->gan[i] = gcnt[i][0] + gcnt[i][1] + gcnt[i][3];
			ss->gac[i][0] = gcnt[i][1];
			ss->gac[i][1] = gcnt[i][3];
		}
	}
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
	// fill bcf1_t up to INFO, excluding AC/AN/etc
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
	// find samples having a set of alleles, or do haplotype counting
	if (bm->n_al > 0) {
		int ret;
		// test if the current record matches an allele
		for (i = 0; i < bm->n_al; ++i) // NOTE: this is a quadratic algorithm; could be better
			if ((ret = bgt_al_test(b, &bm->al[i])) != 0) break;
		if (i == bm->n_al) return 1; // not matching any requested alleles
		// +1 to samples having the allele
		if ((bm->flag&BGT_F_CNT_AL) && bm->alcnt) {
			int is_ref = (ret == 2);
			for (i = 0; i < bm->n_out; ++i) {
				int g1 = bm->a[0][i<<1|0] | bm->a[1][i<<1|0]<<1;
				int g2 = bm->a[0][i<<1|1] | bm->a[1][i<<1|1]<<1;
				if (is_ref) bm->alcnt[i] += (g1 == 0 || g2 == 0);
				else bm->alcnt[i] += (g1 == 1 || g2 == 1);
			}
		}
		// generate haplotype
		if ((bm->flag&BGT_F_CNT_HAP) && bm->hap) {
			for (i = 0; i < bm->n_out<<1; ++i) {
				int h = bm->a[0][i] | bm->a[1][i]<<1;
				if (h == 1) bm->hap[i] |= 1ULL<<bm->n_aal;
			}
			bgt_al_from_bcf(bm->h_out, b, &bm->aal[bm->n_aal++]);
		}
	}
	// fill AC/AN/etc and test site_flt
	if ((bm->flag & BGT_F_SET_AC) || bm->site_flt || bm->n_groups > 1) {
		bgt_info_t ss;
		bgtm_cal_info(bm, &ss);
		bgtm_fill_info(bm->h_out, &ss, b);
		if (!bgtm_pass_site_flt(&ss, bm->site_flt))
			return 1;
	}
	return 0;
}

int bgtm_read(bgtm_t *bm, bcf1_t *b)
{
	int ret;
	if (bm->h_out == 0) bgtm_prepare(bm);
	while ((ret = bgtm_read_core(bm, b)) > 0);
	if ((bm->flag & BGT_F_NO_GT) == 0)
		bgt_gen_gt(bm->h_out, b, bm->n_out, (const uint8_t**)bm->a);
	return ret;
}

/**********************
 * Haplotype counting *
 **********************/

KHASH_MAP_INIT_INT64(hc, int) // WARNING: the default 64-bit hash function is bad

bgt_hapcnt_t *bgtm_hapcnt(const bgtm_t *bm, int *n_hap)
{
	int i, j, absent, n;
	khint_t k;
	khash_t(hc) *h;
	bgt_hapcnt_t *hc;
	if (bm->hap == 0 || bm->n_out == 0) return 0;
	h = kh_init(hc);
	for (i = 0; i < bm->n_out<<1; ++i) {
		k = kh_put(hc, h, bm->hap[i], &absent);
		if (absent) kh_val(h, k) = kh_size(h) - 1;
	}
	n = kh_size(h);
	hc = (bgt_hapcnt_t*)calloc(n, sizeof(bgt_hapcnt_t));
	for (i = 0; i < n; ++i)
		hc[i].cnt = (int*)calloc(bm->n_groups, sizeof(int));
	for (k = 0; k < kh_end(h); ++k)
		if (kh_exist(h, k))
			hc[kh_val(h, k)].hap = kh_key(h, k);
	for (i = 0; i < bm->n_out<<1; ++i) {
		int t;
		k = kh_get(hc, h, bm->hap[i]);
		t = kh_val(h, k);
		for (j = 0; j < bm->n_groups; ++j)
			if (bm->group[i>>1] & 1U<<j)
				++hc[t].cnt[j];
	}
	kh_destroy(hc, h);
	*n_hap = n;
	return hc;
}

char *bgtm_hapcnt_print_destroy(const bgtm_t *bm, int n_hap, bgt_hapcnt_t *hc)
{
	int i, j;
	kstring_t s = {0,0,0};
	ksprintf(&s, "NA\t%d\n", bm->n_aal);
	for (i = 0; i < bm->n_aal; ++i) {
		bgt_allele_t *a = &bm->aal[i];
		ksprintf(&s, "AA\t%s:%d:%d:%s\n", a->chr.s, a->pos+1, a->rlen, a->al);
	}
	ksprintf(&s, "NH\t%d\t%d\n", n_hap, bm->n_groups);
	for (i = 0; i < n_hap; ++i) {
		kputs("HC\t", &s);
		for (j = 0; j < bm->n_aal; ++j)
			kputc('0' + (hc[i].hap>>j&1), &s);
		for (j = 0; j < bm->n_groups; ++j)
			ksprintf(&s, "\t%d", hc[i].cnt[j]);
		kputc('\n', &s);
		free(hc[i].cnt);
	}
	free(hc);
	return s.s;
}

/******************
 * Allele parsing *
 ******************/

int bgt_al_parse(const char *al, bgt_allele_t *a) // TODO: add bgt_al_atomize()
{
	char *p = (char*)al, *ref = 0, *alt = 0;
	int sep = ':', off, tmp, i;
	a->chr.l = 0; a->al = 0; a->pos = -1; a->rlen = -1; a->rid = -1;
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
	} else if (*p == sep) { // a reference allele
		a->rlen = -1;
	}
	if (*p != sep) return -1;
	alt = ++p;
	if (a->rlen < 0) {
		for (i = 0; isalpha(alt[i]); ++i);
		a->rlen = i;
	}
	for (off = 0; *p && isalpha(*p); ++p)
		if (ref && toupper(*p) == toupper(ref[off])) ++off;
		else break;
	a->pos += off; a->rlen -= off;
	tmp = a->chr.l;
	kputs(alt + off, &a->chr);
	a->al = a->chr.s + tmp;
	if (ref) {
		int l_alt = a->chr.s + a->chr.l - a->al;
		int min_l = l_alt < a->rlen? l_alt : a->rlen;
		ref += off;
		for (off = 0; off < min_l && isalpha(ref[a->rlen - 1 - off]) && toupper(ref[a->rlen - 1 - off]) == toupper(a->al[l_alt - 1 - off]); ++off);
		a->rlen -= off;
		a->al[l_alt - off] = 0;
		a->chr.l -= off;
	}
	return 0;
}

int bgt_al_test(const bcf1_t *b, const bgt_allele_t *a) // IMPORTANT: a->rid MUST be set; FIXME: NOT WORKING FOR INDELS!!! Need minize
{
	int l_ref, l_alt, l_al;
	char *ref, *alt;
	if (b->rid != a->rid || b->pos != a->pos || b->rlen != a->rlen) return 0;
	l_al = strlen(a->al);
	bcf_get_ref_alt1(b, &l_ref, &ref, &l_alt, &alt);
	if (l_al == l_alt && strncmp(a->al, alt, l_alt) == 0) return 1;
	if (l_al == l_ref && strncmp(a->al, ref, l_ref) == 0) return 2;
	return 0;
}

void bgt_al_from_bcf(const bcf_hdr_t *h, const bcf1_t *b, bgt_allele_t *a)
{
	char *ref, *alt;
	const char *chr;
	int l_ref, l_alt, min_l;
	bcf_get_ref_alt1(b, &l_ref, &ref, &l_alt, &alt);
	min_l = l_ref < l_alt? l_ref : l_alt;
	a->rid = b->rid, a->pos = b->pos, a->rlen = b->rlen;
	a->chr.l = 0;
	chr = h->id[BCF_DT_CTG][b->rid].key;
	kputs(chr, &a->chr);
	kputc(0, &a->chr);
	kputsn(alt, l_alt, &a->chr);
	a->al = a->chr.s + strlen(chr) + 1;
}
