#ifndef FMF_H
#define FMF_H

#include <stdint.h>
#include <limits.h>
#include "kexpr.h"

#define FMF_FLAG  0
#define FMF_INT   1
#define FMF_REAL  2
#define FMF_STR   3

typedef struct {
	uint32_t key:28, type:4;
	union {
		int32_t i;
		float r;
		uint32_t s;
	} v;
} fmf_meta_t;

typedef struct {
	char *name; // row name
	int n_meta, m_meta;
	fmf_meta_t *meta;
} fmf1_t;

typedef struct {
	int n_keys, m_keys;
	char **keys;
	int n_vals, m_vals;
	char **vals;
	int n_rows, m_rows;
	fmf1_t *rows;
} fmf_t;

struct fms_s;
typedef struct fms_s fms_t;

#ifdef __cplusplus
extern "C" {
#endif

fmf_t *fmf_read(const char *fn);
void fmf_destroy(fmf_t *f);
char *fmf_write(const fmf_t *f, int r);
int fmf_test(const fmf_t *f, int r, kexpr_t *ke);

fms_t *fms_open(const char *fn);
void fms_close(fms_t *f);
const char *fms_read(fms_t *f, kexpr_t *ke, int name_only);

#ifdef __cplusplus
}
#endif

#endif
