#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "fmf.h"
#include "kexpr.h"
#include "kseq.h"
#include "khash.h"

#if defined(FMF_HAVE_HTS)
KSTREAM_DECLARE(gzFile, gzread)
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)
#else
KSTREAM_INIT2(static, gzFile, gzread, 16384)
KHASH_INIT2(s2i, static, kh_cstr_t, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)
#endif

int fmf_verbose = 3;

fmf_t *fmf_read(const char *fn)
{
	kstream_t *ks;
	gzFile fp;
	fmf_t *fmf = 0;
	kstring_t s = {0,0,0};
	int dret;
	khash_t(s2i) *kh, *rh;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	fmf = (fmf_t*)calloc(1, sizeof(fmf_t));
	kh = kh_init(s2i);
	rh = kh_init(s2i);
	while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
		char *p, *q, *r;
		int i, absent, n_meta;
		khint_t k;
		fmf1_t *u = 0;
		if (s.l == 0) continue;
		for (p = s.s, n_meta = 0; *p; ++p)
			if (*p == '\t') ++n_meta;
		for (p = q = s.s, i = 0;; ++p) {
			if (*p == 0 || *p == '\t') {
				int c = *p, c2;
				*p = 0;
				if (i == 0) { // row name
					if (fmf->n_rows == fmf->m_rows) {
						fmf->m_rows = fmf->m_rows? fmf->m_rows<<1 : 16;
						fmf->rows = (fmf1_t*)realloc(fmf->rows, fmf->m_rows * sizeof(fmf1_t));
					}
					u = &fmf->rows[fmf->n_rows++];
					u->name = strdup(q);
					k = kh_put(s2i, rh, u->name, &absent);
					if (!absent) {
						if (fmf_verbose >= 2)
							fprintf(stderr, "[W::%s] row '%s' is duplicated and not query-able.\n", __func__, u->name);
					} else kh_val(rh, k) = fmf->n_rows - 1;
					u->m_meta = n_meta, u->n_meta = 0;
					u->meta = (fmf_meta_t*)calloc(u->m_meta, sizeof(fmf_meta_t));
				} else { // metadata
					fmf_meta_t *m;
					for (r = q; *r && *r != ':'; ++r);
					c2 = *r; *r = 0;
					k = kh_put(s2i, kh, q, &absent);
					if (absent) {
						if (fmf->n_keys == fmf->m_keys) {
							fmf->m_keys = fmf->m_keys? fmf->m_keys<<1 : 8;
							fmf->keys = (char**)realloc(fmf->keys, fmf->m_keys * sizeof(char*));
						}
						kh_val(kh, k) = fmf->n_keys;
						kh_key(kh, k) = fmf->keys[fmf->n_keys++] = strdup(q);
					}
					m = &u->meta[u->n_meta++];
					m->key = kh_val(kh, k), m->v.i = 0;
					if (c2 == ':' && p - r >= 3) {
						if (r[1] == 'i') m->type = FMF_INT, m->v.i = strtol(r + 3, &r, 0);
						else if (r[1] == 'f') m->type = FMF_REAL, m->v.r = strtod(r + 3, &r);
						else m->type = FMF_STR, m->v.s = strdup(r + 3);
					} else m->type = FMF_FLAG;
				}
				q = p + 1; ++i;
				if (c == 0) break; // end-of-line
			}
		}
	}
	free(s.s);
	ks_destroy(ks);
	gzclose(fp);
	fmf->rh = rh, fmf->kh = kh;
	return fmf;
}

void fmf_destroy(fmf_t *f)
{
	int i, j;
	if (f == 0) return;
	for (i = 0; i < f->n_keys; ++i) free(f->keys[i]);
	for (i = 0; i < f->n_rows; ++i) {
		fmf1_t *r = &f->rows[i];
		for (j = 0; j < r->n_meta; ++j)
			if (r->meta[j].type == FMF_STR)
				free(r->meta[j].v.s);
		free(r->name); free(r->meta);
	}
	free(f->rows); free(f->keys);
	kh_destroy(s2i, (khash_t(s2i)*)f->kh);
	kh_destroy(s2i, (khash_t(s2i)*)f->rh);
	free(f);
}

#ifdef FMF_MAIN
int main(int argc, char *argv[])
{
	fmf_t *f;
	if (argc == 1) {
		fprintf(stderr, "Usage: fmf <in.fmf> [condition]\n");
		return 1;
	}
	f = fmf_read(argv[1]);
	fmf_destroy(f);
	return 0;
}
#endif
