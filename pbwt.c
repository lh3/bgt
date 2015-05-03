#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "pbwt.h"

/********************************
 * Run-length encoding/decoding *
 ********************************/

static uint32_t pbr_tbl[128] = {
	0x0U,       0x1U,       0x2U,       0x3U,       0x4U,       0x5U,       0x6U,       0x7U,       0x8U,       0x9U,       0xaU,       0xbU,       0xcU,       0xdU,       0xeU,       0xfU,       
	0x0U,      0x10U,      0x20U,      0x30U,      0x40U,      0x50U,      0x60U,      0x70U,      0x80U,      0x90U,      0xa0U,      0xb0U,      0xc0U,      0xd0U,      0xe0U,      0xf0U,      
	0x0U,     0x100U,     0x200U,     0x300U,     0x400U,     0x500U,     0x600U,     0x700U,     0x800U,     0x900U,     0xa00U,     0xb00U,     0xc00U,     0xd00U,     0xe00U,     0xf00U,     
	0x0U,    0x1000U,    0x2000U,    0x3000U,    0x4000U,    0x5000U,    0x6000U,    0x7000U,    0x8000U,    0x9000U,    0xa000U,    0xb000U,    0xc000U,    0xd000U,    0xe000U,    0xf000U,    
	0x0U,   0x10000U,   0x20000U,   0x30000U,   0x40000U,   0x50000U,   0x60000U,   0x70000U,   0x80000U,   0x90000U,   0xa0000U,   0xb0000U,   0xc0000U,   0xd0000U,   0xe0000U,   0xf0000U,   
	0x0U,  0x100000U,  0x200000U,  0x300000U,  0x400000U,  0x500000U,  0x600000U,  0x700000U,  0x800000U,  0x900000U,  0xa00000U,  0xb00000U,  0xc00000U,  0xd00000U,  0xe00000U,  0xf00000U,  
	0x0U, 0x1000000U, 0x2000000U, 0x3000000U, 0x4000000U, 0x5000000U, 0x6000000U, 0x7000000U, 0x8000000U, 0x9000000U, 0xa000000U, 0xb000000U, 0xc000000U, 0xd000000U, 0xe000000U, 0xf000000U, 
	0x0U,0x10000000U,0x20000000U,0x30000000U,0x40000000U,0x50000000U,0x60000000U,0x70000000U,0x80000000U,0x90000000U,0xa0000000U,0xb0000000U,0xc0000000U,0xd0000000U,0xe0000000U,0xf0000000U
};

static inline int pbr_enc1(uint8_t *p, int l, int b)
{
	if (l >= 16) {
		uint8_t *q = p;
		uint32_t x, i;
		for (x = 0xfU<<28, i = 7<<2; x; x >>= 4, i -= 1<<2)
			if (x&l) *q++ = (i<<2 | (x&l)>>i) << 1 | b;
		return q - p;
	} else {
		*p = l<<1 | b;
		return 1;
	}
}

// rle can be the same as u. In this case, u is overwritten.
static int pbr_enc(int m, const uint8_t *u, uint8_t *rle)
{
	int j, l;
	uint8_t *p = rle, last;
	for (j = 1, l = 1, last = u[j-1]; j < m; ++j) {
		if (u[j] == last) ++l;
		else p += pbr_enc1(p, l, last), l = 1, last = u[j];
	}
	p += pbr_enc1(p, l, last);
	*p = 0;
	return p - rle;
}

/****************************************
 * Basic codec routines for full matrix *
 ****************************************/

// u MUST be at least m+1 long
int pbf_enc(int m, const int32_t *S0, const uint8_t *a, int32_t *S, uint8_t *u)
{
	int32_t *p[2], j, n1;
	for (j = n1 = 0; j < m; ++j)
		n1 += (u[j] = !!a[S0[j]]);
	p[0] = S, p[1] = p[0] + (m - n1);
	for (j = 0; j < m; ++j)
		*p[u[j]]++ = S0[j];
	return pbr_enc(m, u, u);
}

// u MUST be null terminated
void pbf_dec(int m, const int32_t *S0, const uint8_t *u, int32_t *S, uint8_t *a)
{
	const uint8_t *q;
	int32_t *p[2], n1, s;
	for (q = u, n1 = 0; *q; ++q) // count the number of 1 bits
		if (*q&1) n1 += pbr_tbl[*q>>1];
	p[0] = S, p[1] = p[0] + (m - n1);
	for (q = u, s = 0; *q; ++q) {
		int i, l = pbr_tbl[*q>>1], b = *q&1;
		for (i = 0; i < l; ++i) {
			int32_t x = S0[s + i];
			a[x] = b, *p[b]++ = x;
		}
		s += l;
	}
}

/**********************************
 * Basic routine to decode subset *
 **********************************/

#include "ksort.h"
#define pbs_key_r(x) ((x).r)
KRADIX_SORT_INIT(r, pbs_dat_t, pbs_key_r, 4)

void pbs_dec(int m, int r, pbs_dat_t *d, const uint8_t *u)
{
	const uint8_t *q;
	pbs_dat_t *p = d, *end = d + r;
	int n1, c[2], acc[2];
	radix_sort_r(d, d + r); // sort by rank
	for (q = u, n1 = 0; *q; ++q) // count the number of 1 bits
		if (*q&1) n1 += pbr_tbl[*q>>1];
	acc[0] = 0, acc[1] = m - n1; // accumulative counts
	c[0] = c[1] = 0; // running marginal counts
	for (q = u; p != end && *q; ++q) {
		int l = pbr_tbl[*q>>1], b = *q&1, s = c[0] + c[1];
		if (s <= p->r && p->r < s + l) {
			do {
				p->r = acc[b] + c[b] + (p->r - s);
				p->b = b;
				++p;
			} while (p != end && s <= p->r && p->r < s + l);
		}
		c[b] += l;
	}
}

/*****************************
 * More convenient interface *
 *****************************/

pbc_f_t *pbc_f_init(int m)
{
	int j;
	uint8_t *p;
	pbc_f_t *pb;
	p = (uint8_t*)calloc(sizeof(pbc_f_t) + 2 * m * 4 + (m + 1), 1);
	pb = (pbc_f_t*)p; p += sizeof(pbc_f_t);
	pb->S0 = (int32_t*)p; p += m * 4;
	pb->S  = (int32_t*)p; p += m * 4;
	pb->u = p; p += m + 1;
	pb->m = m;
	for (j = 0; j < pb->m; ++j) pb->S[j] = j;
	return pb;
}

void pbc_f_enc(pbc_f_t *pb, const uint8_t *a)
{
	int32_t *swap;
	swap = pb->S, pb->S = pb->S0, pb->S0 = swap;
	pb->l = pbf_enc(pb->m, pb->S0, a, pb->S, pb->u);
}

void pbc_f_dec(pbc_f_t *pb, const uint8_t *b)
{
	int32_t *swap;
	swap = pb->S, pb->S = pb->S0, pb->S0 = swap;
	pbf_dec(pb->m, pb->S0, b, pb->S, pb->u);
}


const int N = 7, M = 4, R = 2;
static uint8_t a[N][M] = {{0,1,0,0}, {0,0,1,1}, {1,0,1,1}, {0,1,0,1}, {1,1,0,0}, {1,0,1,0}, {0,1,1,1}};

#define pbs_key_S(x) ((x).S)
KRADIX_SORT_INIT(S, pbs_dat_t, pbs_key_S, 4)

int main()
{
	pbc_f_t *in, *out;
	pbs_dat_t d[R];
	int k, j;
	d[0].r = 1, d[0].S = 1;
	d[1].r = 2, d[1].S = 2;
	in = pbc_f_init(M);
	out = pbc_f_init(M);
	for (k = 0; k < N; ++k) {
		pbc_f_enc(in, a[k]);
		pbc_f_dec(out, in->u); for (j = 0; j < out->m; ++j) putchar('0' + out->u[j]); putchar('\n');
		pbs_dec(M, R, d, in->u);
		radix_sort_S(d, d + R);
		for (j = 0; j < R; ++j) putchar('0' + d[j].b); putchar('\n');
	}
	free(in); free(out);
	return 0;
}
