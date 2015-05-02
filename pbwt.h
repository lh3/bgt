#ifndef PBWT_H
#define PBWT_H

#include <stdio.h>
#include <stdint.h>

typedef struct {
	int32_t *S0, *S1;
	uint8_t *u;
	int32_t m, dummy;
} pbwt_t;

pbwt_t *pb_init(int m);
void pb_enc(pbwt_t *pb, const uint8_t *a);
void pb_dec_reset(pbwt_t *pb, const int32_t *S);
void pb_dec(pbwt_t *pb, const uint8_t *b);

#endif
