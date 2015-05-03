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

#endif
