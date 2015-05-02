#include <stdlib.h>
#include <string.h>
#include "pbwt.h"

uint32_t pbw_tbl[128];

void pbw_precal_tbl(uint32_t tbl[128])
{
	int i, j, k;
	for (i = k = 0; i < 8; ++i)
		for (j = 0; j < 16; ++j)
			tbl[k++] = (uint32_t)j<<i*4;
}

pbwt_t *pb_init(int m)
{
	int j;
	uint8_t *p;
	pbwt_t *pb;
	if (pbw_tbl[1] == 0) pbw_precal_tbl(pbw_tbl);
	p = (uint8_t*)calloc(sizeof(pbwt_t) + 2 * m * 4 + (m + 1), 1);
	pb = (pbwt_t*)p; p += sizeof(pbwt_t);
	pb->S0 = (int32_t*)p; p += m * 4;
	pb->S1 = (int32_t*)p; p += m * 4;
	pb->u = p; p += m + 1;
	pb->m = m;
	for (j = 0; j < pb->m; ++j) pb->S0[j] = pb->S1[j] = j;
	return pb;
}

// before pb_enc(a_k) is called, S1 = S_{k-1}
// after the call, u = B_k, S0 = S_{k-1} and S1 = S_k
void pb_enc(pbwt_t *pb, const uint8_t *a)
{
	int32_t *p[2], j, n1, *swap;
	swap = pb->S0, pb->S0 = pb->S1, pb->S1 = swap;
	for (j = n1 = 0; j < pb->m; ++j)
		n1 += (pb->u[j] = !!a[pb->S0[j]]);
	p[0] = pb->S1, p[1] = p[0] + (pb->m - n1);
	for (j = 0; j < pb->m; ++j)
		*p[pb->u[j]]++ = pb->S0[j];
}

static inline int pbw_rlenc1(uint8_t *p, int l, int b)
{
	uint8_t *q = p;
	while (l > 16) {
		int i;
		for (i = 17; i < 128 && pbw_tbl[i] <= l; ++i);
		if (l < pbw_tbl[i]) --i;
		*q++ = i<<1 | b;
		l -= pbw_tbl[i];
	}
	if (l > 0) *q++ = l<<1 | b;
	return q - p;
}

int pbw_rlenc(int m, uint8_t *u)
{
	int j, l;
	uint8_t *p = u;
	for (j = 1, l = 1; j < m; ++j) {
		if (u[j] == u[j-1]) ++l;
		else p += pbw_rlenc1(p, l, u[j-1]), l = 1;
	}
	p += pbw_rlenc1(p, l, u[m-1]);
	*p = 0;
	return p - u;
}

static uint8_t a[2][4] = {{0,1,0,0}, {0,1,1,0}};

int main()
{
	pbwt_t *pb;
	int j;
	pb = pb_init(4);
	pb_enc(pb, a[0]);
	pb_enc(pb, a[1]);
	for (j = 0; j < pb->m; ++j) putchar('0'+pb->u[j]);
	for (j = 0; j < pb->m; ++j) printf("\t%d", pb->S1[j]); putchar('\n');
	free(pb);
	return 0;
}
