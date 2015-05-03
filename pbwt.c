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
	pb->S_pre = (int32_t*)p; p += m * 4;
	pb->S_cur = (int32_t*)p; p += m * 4;
	pb->u = p; p += m + 1;
	pb->m = m;
	for (j = 0; j < pb->m; ++j) pb->S_pre[j] = pb->S_cur[j] = j;
	return pb;
}

// before pb_enc(a_k) is called, S_cur = S_{k-1}
// after the call, u = B_k, S_pre = S_{k-1} and S_cur = S_k
void pb_enc(pbwt_t *pb, const uint8_t *a)
{
	int32_t *p[2], j, n1, *swap;
	swap = pb->S_pre, pb->S_pre = pb->S_cur, pb->S_cur = swap;
	for (j = n1 = 0; j < pb->m; ++j)
		n1 += (pb->u[j] = !!a[pb->S_pre[j]]);
	p[0] = pb->S_cur, p[1] = p[0] + (pb->m - n1);
	for (j = 0; j < pb->m; ++j)
		*p[pb->u[j]]++ = pb->S_pre[j];
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

int pbw_rlenc(int m, const uint8_t *u, uint8_t *rle)
{
	int j, l;
	uint8_t *p = rle;
	for (j = 1, l = 1; j < m; ++j) {
		if (u[j] == u[j-1]) ++l;
		else p += pbw_rlenc1(p, l, u[j-1]), l = 1;
	}
	p += pbw_rlenc1(p, l, u[m-1]);
	*p = 0;
	return p - rle;
}

int pbw_rldec(const uint8_t *rle, uint8_t *u)
{
	const uint8_t *p;
	uint8_t *q = u;
	for (p = rle; *p; ++p) {
		int i, l = pbw_tbl[*p>>1], b = *p&1;
		for (i = 0; i < l; ++i) *q++ = b;
	}
	return q - u;
}

void pb_dec_reset(pbwt_t *pb, const int32_t *S)
{
	int i;
	for (i = 0; i < pb->m; ++i)
		pb->S_cur[i] = S? S[i] : i;
}

static uint8_t a[2][4] = {{0,1,0,0}, {0,0,1,1}};

int main()
{
	pbwt_t *pb;
	int k, j;
	pb = pb_init(4);
	for (k = 0; k < 2; ++k) {
		pb_enc(pb, a[k]);
		for (j = 0; j < pb->m; ++j) putchar('0'+pb->u[j]);
		for (j = 0; j < pb->m; ++j) printf("\t%d", pb->S_cur[j]); putchar('\n');
	}
	free(pb);
	return 0;
}
