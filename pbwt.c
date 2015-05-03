#include <stdlib.h>
#include <string.h>
#include "pbwt.h"

/********************************
 * Run-length encoding/decoding *
 ********************************/

uint32_t pbw_tbl[128];

void pbw_precal_tbl(uint32_t tbl[128])
{
	int i, j, k;
	for (i = k = 0; i < 8; ++i)
		for (j = 0; j < 16; ++j)
			tbl[k++] = (uint32_t)j<<i*4;
}

static inline int pbw_rlenc1(uint8_t *p, int l, int b)
{
	if (l >= 16) {
		uint8_t *q = p;
		uint32_t x, i;
		for (x = 0xfU<<28, i = 7<<2; x; x >>= 4, i -= 1<<2)
			if (x&l) *q++ = (i<<2 | (x&l)>>i) << 1 | b;
		*q = 0;
		return q - p;
	} else {
		*p = l<<1 | b; p[1] = 0;
		return 1;
	}
}

// rle can be the same as u. In this case, u is overwritten.
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

/******************
 * Positional BWT *
 ******************/

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
	pb->l = pbw_rlenc(pb->m, pb->u, pb->u);
}

void pb_dec_reset(pbwt_t *pb, const int32_t *S)
{
	int i;
	for (i = 0; i < pb->m; ++i)
		pb->S_cur[i] = S? S[i] : i;
}

void pb_dec_all(pbwt_t *pb, const uint8_t *b)
{
	const uint8_t *q;
	int32_t *swap, *p[2], s = 0;
	int i, l, c, n1;
	swap = pb->S_pre, pb->S_pre = pb->S_cur, pb->S_cur = swap;
	for (q = b, n1 = 0; *q; ++q) // count the number of 1 bits
		if (*q&1) n1 += pbw_tbl[*q>>1];
	p[0] = pb->S_cur, p[1] = p[0] + (pb->m - n1);
	for (q = b; *q; ++q) {
		l = pbw_tbl[*q>>1], c = *q&1;
		for (i = 0; i < l; ++i) {
			int32_t S = pb->S_pre[s + i];
			pb->u[S] = c, *p[c]++ = S;
		}
		s += l;
	}
}

const int N = 5, M = 4;
static uint8_t a[N][M] = {{0,1,0,0}, {0,0,1,1}, {1,0,1,1}, {0,1,0,1}, {1,1,0,0}};

int main()
{
	pbwt_t *in, *out;
	int k, j;
	in = pb_init(M);
	out = pb_init(M);
	for (j = 0; j < 8; ++j) {
		for (k = 0; k < 16; ++k)
			printf("%8x,", pbw_tbl[j*16+k]);
		putchar('\n');
	}
	for (k = 0; k < N; ++k) {
		pb_enc(in, a[k]);
		pb_dec_all(out, in->u);
		for (j = 0; j < out->m; ++j) putchar('0' + out->u[j]); putchar('\n');
//		for (j = 0; j < pb->m; ++j) putchar('0'+in->u[j]);
//		for (j = 0; j < pb->m; ++j) printf("\t%d", in->S_cur[j]); putchar('\n');
	}
	free(in); free(out);
	return 0;
}
