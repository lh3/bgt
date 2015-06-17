#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
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

// encode one run
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

// encode a binary string with RLE. $rle can be the same as $u. In this case, $u is overwritten.
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

/*****************************
 * Encode/decode all columns *
 *****************************/

// Given S_{k-1} and A_k, derive B_k and S_k. $u MUST be at least m+1 long.
int pbc_enc_core(int m, const int32_t *S0, const uint8_t *a, int32_t *S, uint8_t *u)
{
	int32_t *p[2], j, n1;
	for (j = n1 = 0; j < m; ++j) // count the number of 1 bits
		n1 += (u[j] = !!a[S0[j]]);
	p[0] = S, p[1] = p[0] + (m - n1);
	for (j = 0; j < m; ++j)
		*p[u[j]]++ = S0[j];
	return pbr_enc(m, u, u);
}

// Given S_{k-1} and B_k, derive A_k and S_k. $u MUST be null terminated.
void pbc_dec_core(int m, const int32_t *S0, const uint8_t *u, int32_t *S, uint8_t *a)
{
	const uint8_t *q;
	int32_t *p[2], n1, s;
	for (q = u, n1 = 0; *q; ++q) // count the number of 1 bits
		if (*q&1) n1 += pbr_tbl[*q>>1];
	p[0] = S, p[1] = p[0] + (m - n1);
	for (q = u, s = 0; *q; ++q) {
		int i, l = pbr_tbl[*q>>1], b = *q&1; // $l is the run length
		for (i = 0; i < l; ++i) {
			int32_t x = S0[s + i];
			a[x] = b, *p[b]++ = x;
		}
		s += l;
	}
}

pbc_t *pbc_init(int m)
{
	int j;
	uint8_t *p;
	pbc_t *pb;
	p = (uint8_t*)calloc(sizeof(pbc_t) + 2 * m * 4 + (m + 1), 1);
	pb = (pbc_t*)p; p += sizeof(pbc_t);
	pb->S0 = (int32_t*)p; p += m * 4;
	pb->S  = (int32_t*)p; p += m * 4;
	pb->u = p; p += m + 1;
	pb->m = m;
	for (j = 0; j < pb->m; ++j) pb->S[j] = j;
	return pb;
}

void pbc_enc(pbc_t *pb, const uint8_t *a)
{
	int32_t *swap;
	swap = pb->S, pb->S = pb->S0, pb->S0 = swap;
	pb->l = pbc_enc_core(pb->m, pb->S0, a, pb->S, pb->u);
}

void pbc_dec(pbc_t *pb, const uint8_t *b)
{
	int32_t *swap;
	swap = pb->S, pb->S = pb->S0, pb->S0 = swap;
	pbc_dec_core(pb->m, pb->S0, b, pb->S, pb->u);
}

/******************************
 * Decode a subset of columns *
 ******************************/

#include "ksort.h"
#define pbs_key_r(x) ((x).r)
KRADIX_SORT_INIT(r, pbs_dat_t, pbs_key_r, 4)

void pbs_dec(int m, int r, pbs_dat_t *d, const uint8_t *u) // IMPORTANT: d MUST BE sorted by d[i].r
{
	const uint8_t *q;
	pbs_dat_t *p = d, *end = d + r, *swap, *x[2];
	int n1, c[2], acc[2];
//	radix_sort_r(d, d + r); // sort by rank
	for (q = u, n1 = 0; *q; ++q) // count the number of 1 bits
		if (*q&1) n1 += pbr_tbl[*q>>1];
	acc[0] = 0, acc[1] = m - n1; // accumulative counts
	c[0] = c[1] = 0; // running marginal counts
	for (q = u, n1 = 0; p != end && *q; ++q) {
		int l = pbr_tbl[*q>>1], b = *q&1, s = c[0] + c[1];
		if (s <= p->r && p->r < s + l) {
			do {
				p->r = acc[b] + c[b] + (p->r - s);
				p->b = b;
				n1 += b;
				++p;
			} while (p != end && s <= p->r && p->r < s + l);
		}
		c[b] += l;
	}
	// one round of radix sort; this sorts d by d[i].r
	swap = (pbs_dat_t*)malloc(r * sizeof(pbs_dat_t));
	memcpy(swap, d, r * sizeof(pbs_dat_t));
	x[0] = d, x[1] = d + (r - n1);
//	The following 3 lines are the loop fission of: for (i = 0; i < r; ++i) *x[swap[i].b]++ = swap[i];
	end = swap + r;
	for (p = swap; p != end; ++p) if (p->b == 0) *x[0]++ = *p;
	for (p = swap; p != end; ++p) if (p->b != 0) *x[1]++ = *p;
	free(swap);
}

/************
 * File I/O *
 ************/

struct pbf_s {
	FILE *fp;   // PBF file handler
	int32_t m;  // number of columns
	int32_t g;  // number of bits per group
	int32_t shift; // insert S every 1<<shift rows
	int32_t is_writing; // file opend for writing
	int64_t n;  // number of rows

	pbc_t **pb; // pbwt full codecs
	const uint8_t **ret; // ret[g] points to pb[g]->u; this is for return (writing only)

	int32_t n_idx, m_idx;
	uint64_t *idx; // file offset of "S" records

	int n_sub;
	pbs_dat_t **sub;
	int *sub_pos;

	int64_t k;     // the row index just processed (reading only)
	uint8_t *buf;  // reading only
	int32_t *invS; // reading only
};

pbf_t *pbf_open_w(const char *fn, int m, int g, int shift)
{
	FILE *fp;
	pbf_t *pb;
	int32_t i, v[3];
	if (fn && strcmp(fn, "-") != 0) {
		if ((fp = fopen(fn, "wb")) == NULL)
			return 0;
	} else fp = stdout;
	pb = (pbf_t*)calloc(1, sizeof(pbf_t));
	pb->fp = fp;
	pb->m = m, pb->g = g, pb->shift = shift;
	pb->pb = (pbc_t**)calloc(g, sizeof(void*));
	for (i = 0; i < g; ++i)
		pb->pb[i] = pbc_init(m);
	v[0] = pb->m, v[1] = pb->g, v[2] = pb->shift;
	fwrite("PBF\1", 1, 4, fp);
	fwrite(v, 4, 3, fp);
	pb->is_writing = 1;
	return pb;
}

pbf_t *pbf_open_r(const char *fn)
{
	pbf_t *pb;
	FILE *fp;
	int32_t i, v[3];
	char magic[4];
	if (fn && strcmp(fn, "-") != 0) {
		if ((fp = fopen(fn, "rb")) == 0)
			return 0;
	} else fp = stdin;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, "PBF\1", 4) != 0) {
		fclose(fp);
		return 0;
	}
	pb = (pbf_t*)calloc(1, sizeof(pbf_t));
	fread(v, 4, 3, fp);
	pb->m = v[0], pb->g = v[1], pb->shift = v[2];
	pb->pb = (pbc_t**)calloc(pb->g, sizeof(void*));
	for (i = 0; i < pb->g; ++i)
		pb->pb[i] = pbc_init(pb->m);
	pb->buf = (uint8_t*)calloc(pb->m + 1, 1);
	pb->invS = (int32_t*)calloc(pb->m, 4);
	pb->ret = (const uint8_t**)calloc(pb->g, sizeof(uint8_t*));
	for (i = 0; i < pb->g; ++i) pb->ret[i] = pb->pb[i]->u;
	pb->sub = (pbs_dat_t**)calloc(pb->g, sizeof(pbs_dat_t*));
	if (fseek(fp, -8, SEEK_END) >= 0) {
		uint64_t off;
		uint8_t t;
		fread(&off, 8, 1, fp);
		fseek(fp, off, SEEK_SET);
		fread(&t, 1, 1, fp);
		fread(&pb->n, 8, 1, fp);
		fread(&pb->n_idx, 4, 1, fp);
		pb->m_idx = pb->n_idx;
		pb->idx = (uint64_t*)calloc(pb->n_idx, 8);
		fread(pb->idx, 8, pb->n_idx, fp);
		fseek(fp, 16, SEEK_SET);
	}
	pb->fp = fp;
	return pb;
}

int pbf_close(pbf_t *pb)
{
	int g;
	if (pb == 0) return 0;
	if (pb->is_writing) { // write the index
		uint64_t off;
		off = ftell(pb->fp);
		fputc('I', pb->fp);
		fwrite(&pb->n, 8, 1, pb->fp);
		fwrite(&pb->n_idx, 4, 1, pb->fp);
		fwrite(pb->idx, 8, pb->n_idx, pb->fp);
		fwrite(&off, 8, 1, pb->fp);
	}
	free(pb->idx); free(pb->ret); free(pb->invS); free(pb->buf); free(pb->sub_pos);
	for (g = 0; g < pb->g; ++g) {
		free(pb->pb[g]);
		if (pb->sub) free(pb->sub[g]);
	}
	free(pb->sub); free(pb->pb);
	fclose(pb->fp);
	free(pb);
	return 0;
}

int pbf_write(pbf_t *pb, uint8_t *const*a)
{
	int g;
	if (!pb->is_writing) return -1;
	if ((pb->n & ((1ULL<<pb->shift) - 1)) == 0) {
		if (pb->n_idx == pb->m_idx) {
			pb->m_idx = pb->m_idx? pb->m_idx<<1 : 8;
			pb->idx = (uint64_t*)realloc(pb->idx, pb->m_idx * 8);
		}
		pb->idx[pb->n_idx++] = ftell(pb->fp); // save the index offset
		fputc('S', pb->fp);
		for (g = 0; g < pb->g; ++g) // write S[]
			fwrite(pb->pb[g]->S, 4, pb->m, pb->fp);
	}
	fputc('B', pb->fp);
	for (g = 0; g < pb->g; ++g) {
		pbc_t *pbc = pb->pb[g];
		pbc_enc(pbc, a[g]);
		fwrite(&pbc->l, 4, 1, pb->fp);
		fwrite(pbc->u, 1, pbc->l, pb->fp);
	}
	++pb->n;
	return 0;
}

const uint8_t **pbf_read(pbf_t *pb)
{
	int g;
	uint8_t t;
	if (pb->is_writing) return 0;
	fread(&t, 1, 1, pb->fp);
	if (t == 'S') {
		for (g = 0; g < pb->g; ++g)
			fread(pb->pb[g]->S, 4, pb->m, pb->fp);
		fread(&t, 1, 1, pb->fp);
	}
	if (t == 'B') {
		for (g = 0; g < pb->g; ++g) {
			int32_t l, i;
			fread(&l, 4, 1, pb->fp);
			fread(pb->buf, 1, l, pb->fp);
			pb->buf[l] = 0;
			if (pb->n_sub > 0 && pb->n_sub < pb->m) { // subset decoding
				pbs_dat_t *sub = pb->sub[g];
				pbs_dec(pb->m, pb->n_sub, sub, pb->buf);
				for (i = 0; i < pb->n_sub; ++i)
					pb->pb[g]->u[pb->sub_pos[sub[i].S]] = sub[i].b;
			} else pbc_dec(pb->pb[g], pb->buf); // full decoding
		}
		++pb->k;
	} else return 0;
	return pb->ret;
}

// find the rank of a subset of columns given S
static inline void pbf_fill_sub(int m, const int32_t *S, int n_sub, pbs_dat_t *sub, int32_t *invS)
{
	int i;
	for (i = 0; i < m; ++i) invS[S[i]] = i;
	for (i = 0; i < n_sub; ++i)
		sub[i].r = invS[sub[i].S];
	radix_sort_r(sub, sub + n_sub);
}

int pbf_seek(pbf_t *pb, uint64_t k)
{
	int x, i, g;
	uint8_t t;
	if (pb->is_writing) return -1;
	if (k == pb->k) return 0;
	if (k > pb->k && k - pb->k <= 1<<pb->shift) {
		while (pb->k < k) pbf_read(pb);
		return 0;
	}
	if (pb->idx == 0 || k >= pb->n) return -1;
	fseek(pb->fp, pb->idx[k>>pb->shift], SEEK_SET);
	fread(&t, 1, 1, pb->fp);
	assert(t == 'S'); // a bug or corrupted file if it is not an "S" line
	for (g = 0; g < pb->g; ++g) {
		fread(pb->pb[g]->S, 4, pb->m, pb->fp);
		if (pb->n_sub > 0 && pb->n_sub < pb->m) // update pb->sub if needed
			pbf_fill_sub(pb->m, pb->pb[g]->S, pb->n_sub, pb->sub[g], pb->invS);
	}
	pb->k = k >> pb->shift << pb->shift;
	x = k & ((1<<pb->shift) - 1);
	for (i = 0; i < x; ++i) pbf_read(pb);
	return 0;
}

int pbf_subset(pbf_t *pb, int n_sub, int *sub)
{
	int i, g;
	if (n_sub <= 0 || n_sub >= pb->m || sub == 0) n_sub = 0;
	if ((pb->n_sub = n_sub) != 0) {
		if (pb->sub_pos == 0)
			pb->sub_pos = (int*)calloc(pb->m, sizeof(int));
		for (i = 0; i < n_sub; ++i) pb->sub_pos[sub[i]] = i;
		for (g = 0; g < pb->g; ++g) {
			pb->sub[g] = (pbs_dat_t*)realloc(pb->sub[g], n_sub * sizeof(pbs_dat_t));
			for (i = 0; i < n_sub; ++i) pb->sub[g][i].S = sub[i];
			pbf_fill_sub(pb->m, pb->pb[g]->S, n_sub, pb->sub[g], pb->invS);
		}
	}
	return 0;
}

int pbf_get_g(const pbf_t *pb) { return pb->g; }
int pbf_get_m(const pbf_t *pb) { return pb->m; }
int pbf_get_n(const pbf_t *pb) { return pb->n; }
int pbf_get_shift(const pbf_t *pb) { return pb->shift; }
