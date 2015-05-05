#ifndef PBWT_H
#define PBWT_H

#include <stdint.h>

typedef struct { // full codec
	int32_t m, l, *S0, *S;
	uint8_t *u;
} pbc_t;

typedef struct {
	uint32_t r;
	uint32_t S:31, b:1;
} pbs_dat_t;

struct pbf_s;
typedef struct pbf_s pbf_t;

#ifdef __cplusplus
extern "C" {
#endif

/***********************************
 * File-based high-level functions *
 ***********************************/

pbf_t *pbf_open_w(const char *fn, int m, int g, int shift);
pbf_t *pbf_open_r(const char *fn);
int pbf_close(pbf_t *pb);
int pbf_write(pbf_t *pb, uint8_t *const*a);
const uint8_t **pbf_read(pbf_t *pb);
int pbf_seek(pbf_t *pb, uint64_t k);
int pbf_subset(pbf_t *fp, int n_sub, int *sub);

/***********************
 * Low-level functions *
 ***********************/

pbc_t *pbc_init(int m);
void pbc_enc(pbc_t *pb, const uint8_t *a);
void pbc_dec(pbc_t *pb, const uint8_t *b);

void pbs_dec(int m, int r, pbs_dat_t *d, const uint8_t *u);

#ifdef __cplusplus
extern "C" {
#endif

#endif
