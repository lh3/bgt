#ifndef PBWT_H
#define PBWT_H

typedef struct { // full codec
	int32_t m, l, *S0, *S;
	uint8_t *u;
} pbc_f_t;

#endif
