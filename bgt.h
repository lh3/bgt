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
	bcf1_t *b;
	uint8_t *bit[2];
	int row, copied;
} bgt_raw1_t;

typedef struct {
	int rid, pos;
	int m_b, n_b;
	bgt_raw1_t *b, next;
} bgt_rawpos_t;

typedef struct {
	int rid, start, end;
	kstring_t allele;
} bgtm_allele_t;

typedef struct {
	int n_bgt;
	bgt_t **bgt;
	bgt_rawpos_t *rp;
} bgtm_t;

bgt_t *bgt_open(const char *prefix);
void bgt_close(bgt_t *bgt);
void bgt_set_samples(bgt_t *bgt, int n, char *const* samples);
int bgt_set_region(bgt_t *bgt, const char *reg);
int bgt_read(bgt_t *bgt, bcf1_t *b);

bgtm_t *bgtm_open(int n_files, char **fns);
void bgtm_close(bgtm_t *bm);

#endif
