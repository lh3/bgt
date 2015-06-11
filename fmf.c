#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "fmf.h"
#include "kseq.h"
#include "khash.h"
#include "kstring.h"

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
	khash_t(s2i) *kh, *vh;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	fmf = (fmf_t*)calloc(1, sizeof(fmf_t));
	kh = kh_init(s2i);
	vh = kh_init(s2i);
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
						else {
							char *s = r + 3;
							m->type = FMF_STR;
							k = kh_put(s2i, vh, s, &absent);
							if (absent) {
								if (fmf->n_vals == fmf->m_vals) {
									fmf->m_vals = fmf->m_vals? fmf->m_vals<<1 : 8;
									fmf->vals = (char**)realloc(fmf->vals, fmf->m_vals * sizeof(char*));
								}
								kh_val(vh, k) = fmf->n_vals;
								kh_key(vh, k) = fmf->vals[fmf->n_vals++] = strdup(s);
							}
							m->v.s = kh_val(vh, k);
						}
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
	kh_destroy(s2i, kh);
	kh_destroy(s2i, vh);
	return fmf;
}

void fmf_destroy(fmf_t *f)
{
	int i;
	if (f == 0) return;
	for (i = 0; i < f->n_keys; ++i) free(f->keys[i]);
	for (i = 0; i < f->n_vals; ++i) free(f->vals[i]);
	for (i = 0; i < f->n_rows; ++i) {
		fmf1_t *r = &f->rows[i];
		free(r->name); free(r->meta);
	}
	free(f->rows); free(f->keys); free(f->vals);
	free(f);
}

char *fmf_write(const fmf_t *f, int r)
{
	static char *type_str = "\0ifZ";
	kstring_t s = {0,0,0};
	fmf1_t *u;
	int i;
	if (r >= f->n_rows) return 0;
	u = &f->rows[r];
	kputs(u->name, &s);
	for (i = 0; i < u->n_meta; ++i) {
		fmf_meta_t *m = &u->meta[i];
		kputc('\t', &s);
		kputs(f->keys[m->key], &s);
		if (m->type != FMF_FLAG) {
			kputc(':', &s); kputc(type_str[m->type], &s);
			kputc(':', &s);
			if (m->type == FMF_INT) ksprintf(&s, "%lld", (long long)m->v.i);
			else if (m->type == FMF_REAL) ksprintf(&s, "%g", m->v.r);
			else kputs(f->vals[m->v.s], &s);
		}
	}
	return s.s;
}

int fmf_test(const fmf_t *f, int r, kexpr_t *ke) // FIXME: a quadratic implementation!!! optimize later if too slow
{
	fmf1_t *u;
	int err, i, is_true;
	if (r >= f->n_rows) return 0;
	u = &f->rows[r];
	ke_unset(ke);
	for (i = 0; i < u->n_meta; ++i) {
		fmf_meta_t *m = &u->meta[i];
		ke_set_str(ke, "_ROW_", u->name);
		if (m->type == FMF_STR) ke_set_str(ke, f->keys[m->key], f->vals[m->v.s]);
		else if (m->type == FMF_INT) ke_set_int(ke, f->keys[m->key], m->v.i);
		else if (m->type == FMF_REAL) ke_set_int(ke, f->keys[m->key], m->v.r);
	}
	is_true = !!ke_eval_int(ke, &err);
	return !(err || !is_true);
}

struct fms_s {
	kstream_t *ks;
	kstring_t s;
	gzFile fp;
};

fms_t *fms_open(const char *fn)
{
	fms_t *f;
	gzFile fp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	f = (fms_t*)calloc(1, sizeof(fms_t));
	f->ks = ks_init(fp);
	f->fp = fp;
	return f;
}

void fms_close(fms_t *f)
{
	if (f == 0) return;
	free(f->s.s);
	ks_destroy(f->ks);
	gzclose(f->fp);
	free(f);
}

static int fms_read_and_test(fms_t *f, kexpr_t *ke, char **end0)
{
	char *p, *q, *r, *rr;
	int i, err = 0, is_true, dret, ret;
	ret = ks_getuntil(f->ks, KS_SEP_LINE, &f->s, &dret);
	if (ret < 0) return ret;
	if (f->s.l == 0) return 0;
	if (ke) ke_unset(ke);
	for (p = q = f->s.s, i = 0;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p, c2;
			*p = 0;
			if (i == 0) { // row name
				if (ke) ke_set_str(ke, "_ROW_", q);
				*end0 = p;
			} else { // metadata
				for (r = q; *r && *r != ':'; ++r);
				c2 = *r; *r = 0; rr = r;
				if (c2 == ':' && p - r >= 3 && ke) {
					if (r[1] == 'i') ke_set_int(ke, q, strtol(r + 3, &r, 0));
					else if (r[1] == 'f') ke_set_real(ke, q, strtod(r + 3, &r));
					else ke_set_str(ke, q, r + 3);
				}
				*rr = c2;
			}
			q = p + 1; ++i;
			if (c == 0) break; // end-of-line
			*p = c;
		}
	}
	is_true = ke == 0 || !!ke_eval_int(ke, &err);
	return (!err && is_true);
}

const char *fms_read(fms_t *f, kexpr_t *ke, int name_only)
{
	int ret;
	char *end0 = 0;
	while ((ret = fms_read_and_test(f, ke, &end0)) == 0);
	if (ret < 0) return 0;
	if (name_only) *end0 = 0;
	return f->s.s;
}

#ifndef FMF_LIB_ONLY
#include <unistd.h>

int main_fmf(int argc, char *argv[])
{
	kexpr_t *ke = 0;
	int i, c, err, in_mem = 0, name_only = 0;
	while ((c = getopt(argc, argv, "mn")) >= 0)
		if (c == 'm') in_mem = 1;
		else if (c == 'n') name_only = 1;
	if (argc == optind) {
		fprintf(stderr, "Usage: fmf [-mn] <in.fmf> [condition]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m   load the entire FMF into RAM\n");
		fprintf(stderr, "  -n   only output the row name (the 1st column)\n");
		return 1;
	}
	if (argc - optind >= 2) ke = ke_parse(argv[optind+1], &err);
	if (in_mem) {
		fmf_t *f;
		f = fmf_read(argv[optind]);
		for (i = 0; i < f->n_rows; ++i) {
			char *s;
			if (ke && !fmf_test(f, i, ke)) continue;
			if (!name_only) {
				s = fmf_write(f, i);
				puts(s);
				free(s);
			} else puts(f->rows[i].name);
		}
		fmf_destroy(f);
	} else {
		fms_t *f;
		const char *s;
		f = fms_open(argv[optind]);
		while ((s = fms_read(f, ke, name_only)) != 0)
			puts(s);
		fms_close(f);
	}
	if (ke) ke_destroy(ke);
	return 0;
}

#ifdef FMF_MAIN
int main(int argc, char *argv[]) { return main_fmf(argc, argv); }
#endif

#endif
