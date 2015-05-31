#ifndef BGT_H
#define BGT_H

#include "vcf.h"
#include "pbwt.h"

#define BGT_F_SET_AC    0x0001
#define BGT_F_NO_GT     0x0002

#define BGT_MAX_GROUPS  8

typedef struct {
	int n_samples;
	char **samples;

	pbf_t *pb;
	htsFile *bcf;
	bcf_hdr_t *h0; // site-only BCF header
	bcf1_t *b0; // site-only BCF record
	hts_idx_t *idx; // BCF index
	hts_itr_t *itr;
	const void *bed;

	int bed_excl, n_out, *out;
	bcf_hdr_t *h_out;
} bgt_t;

typedef struct { // during reading, these are all links
	const bcf1_t *b0;
	const uint8_t *a[2];
} bgt_rec_t;

typedef int (*bgt_filter_f)(bcf_hdr_t *h, bcf1_t *b, int an, int ac, int n_groups, int32_t *gan, int32_t *gac1, void *data);

typedef struct {
	int n_bgt, n_out, n_groups, flag;
	uint64_t *sample_idx;
	uint8_t *group;
	bgt_t **bgt;
	bgt_rec_t *r;
	bcf_hdr_t *h_out;
	bgt_filter_f filter_func;
	void *filter_data;
	uint8_t *a[2];
} bgtm_t;

typedef struct {
	kstring_t chr;
	char *alt;
	int pos, rlen;
} bgt_allele_t;

bgt_t *bgt_open(const char *prefix);
void bgt_close(bgt_t *bgt);
void bgt_set_samples(bgt_t *bgt, int n, char *const* samples);
void bgt_set_bed(bgt_t *bgt, const void *bed, int excl);
int bgt_set_region(bgt_t *bgt, const char *reg);
int bgt_read(bgt_t *bgt, bcf1_t *b);

bgtm_t *bgtm_open(int n_files, char *const*fns);
void bgtm_close(bgtm_t *bm);
void bgtm_set_flag(bgtm_t *bm, int flag);
void bgtm_set_samples(bgtm_t *bm, int n, char *const* samples);
void bgtm_set_filter(bgtm_t *bm, bgt_filter_f flt, void *flt_data);
void bgtm_set_bed(bgtm_t *bm, const void *bed, int excl);
int bgtm_set_region(bgtm_t *bm, const char *reg);
void bgtm_add_group(bgtm_t *bm, int n, char *const* samples);
int bgtm_read(bgtm_t *bm, bcf1_t *b);

int bgt_al_parse(const char *al, bgt_allele_t *a);
int bgt_al_test(const bcf_hdr_t *h, const bcf1_t *b, const bgt_allele_t *a);

#endif
