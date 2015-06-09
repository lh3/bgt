#ifndef BGT_H
#define BGT_H

#include "vcf.h"
#include "pbwt.h"
#include "fmf.h"

#define BGT_F_SET_AC    0x0001
#define BGT_F_NO_GT     0x0002
#define BGT_F_CNT_AL    0x0004
#define BGT_F_CNT_HAP   0x0008

#define BGT_MAX_GROUPS  32
#define BGT_MAX_ALLELES 64

#define BGT_SET_ALL_SAMPLES (-1)

typedef struct {
	char *prefix;
	fmf_t *f;
	bcf_hdr_t *h0; // site-only BCF header
	hts_idx_t *idx; // BCF index
} bgt_file_t;

typedef struct {
	const bgt_file_t *f;
	pbf_t *pb;
	BGZF *bcf;
	bcf1_t *b0; // site-only BCF record
	hts_itr_t *itr;
	const void *bed;
	int bed_excl, n_out, n_groups, *out;
	uint32_t *group, *flag;
	bcf_hdr_t *h_out;
} bgt_t;

typedef struct { // during reading, these are all links
	const bcf1_t *b0;
	const uint8_t *a[2];
} bgt_rec_t;

typedef struct {
	int32_t ac[2], an, n_groups;
	int32_t gan[BGT_MAX_GROUPS], gac[BGT_MAX_GROUPS][2];
} bgt_info_t;

typedef struct {
	kstring_t chr;
	char *al;
	int rid, pos, rlen;
} bgt_allele_t;

typedef struct {
	uint64_t hap;
	int *cnt;
} bgt_hapcnt_t;

typedef struct {
	int n_bgt, n_out, n_groups, flag;
	uint64_t *sample_idx;
	uint32_t *group;
	bgt_t **bgt;
	bgt_rec_t *r;
	kexpr_t *site_flt;
	bcf_hdr_t *h_out;
	uint8_t *a[2];

	int n_fields;
	kexpr_t **fields;
	kstring_t tbl_line;

	int n_aal;
	bgt_allele_t *aal;
	void *h_al;
	int *alcnt;
	uint64_t *hap;
} bgtm_t;

#ifdef __cplusplus
extern "C" {
#endif

bgt_file_t *bgt_open(const char *prefix);
void bgt_close(bgt_file_t *bgt);

bgt_t *bgt_reader_init(const bgt_file_t *bf);
void bgt_reader_destroy(bgt_t *bgt);
void bgt_set_bed(bgt_t *bgt, const void *bed, int excl);
int bgt_set_region(bgt_t *bgt, const char *reg);
int bgt_set_start(bgt_t *bgt, int64_t n);

int bgt_read(bgt_t *bgt, bcf1_t *b);

bgtm_t *bgtm_reader_init(int n_files, bgt_file_t *const*fns);
void bgtm_reader_destroy(bgtm_t *bm);
void bgtm_set_flag(bgtm_t *bm, int flag);
int bgtm_set_flt_site(bgtm_t *bm, const char *expr);
void bgtm_set_bed(bgtm_t *bm, const void *bed, int excl);
int bgtm_set_region(bgtm_t *bm, const char *reg);
int bgtm_set_start(bgtm_t *bm, int64_t n);
int bgtm_set_table(bgtm_t *bm, const char *fmt);
int bgtm_set_alleles(bgtm_t *bm, const char *expr, const fmf_t *f, const char *fn); // call this AFTER bgtm_set_region()
void bgtm_add_group_core(bgtm_t *bm, int n, char *const* samples, const char *expr);
void bgtm_add_group(bgtm_t *bm, const char *expr);
int bgtm_add_allele(bgtm_t *bm, const char *al);
void bgtm_prepare(bgtm_t *bm);

int bgtm_read(bgtm_t *bm, bcf1_t *b);

bgt_hapcnt_t *bgtm_hapcnt(const bgtm_t *bm, int *n_hap);
char *bgtm_hapcnt_print_destroy(const bgtm_t *bm, int n_hap, bgt_hapcnt_t *hc);
char *bgtm_alcnt_print(const bgtm_t *bm);

int bgt_al_parse(const char *al, bgt_allele_t *a);
void bgt_al_format(const bgt_allele_t *a, kstring_t *s);
void bgt_al_from_bcf(const bcf_hdr_t *h, const bcf1_t *b, bgt_allele_t *a, bgt_allele_t *r);

#ifdef __cplusplus
}
#endif

#endif
