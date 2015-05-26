#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include "bgt.h"
#include "kexpr.h"
#include "fmf.h"

#include "khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)

void *bed_read(const char *fn);
void bed_destroy(void *_h);

typedef struct {
	kexpr_t *ke;
} flt_aux_t;

static inline int read1(bgt_t *bgt, bgtm_t *bm, bcf1_t *b)
{
	if (bgt) return bgt_read(bgt, b);
	if (bm)  return bgtm_read(bm, b);
	return -1;
}

static int filter_func(bcf_hdr_t *h, bcf1_t *b, int an, int ac1, int n_groups, int32_t *gan, int32_t *gac1, void *data)
{
	flt_aux_t *flt = (flt_aux_t*)data;
	int is_true, err, i;
	char key[4];
	if (flt->ke == 0) return 0;
	ke_set_int(flt->ke, "AN", an);
	ke_set_int(flt->ke, "AC", ac1);
	key[3] = 0;
	for (i = 0; i <= n_groups; ++i) {
		key[0] = 'A', key[1] = 'N', key[2] = '0' + i; ke_set_int(flt->ke, key, gan[i]);
		key[0] = 'A', key[1] = 'C', key[2] = '0' + i; ke_set_int(flt->ke, key, gac1[i]);
	}
	is_true = !!ke_eval_int(flt->ke, &err);
	if (err) return 0;
	return !is_true;
}

static char **get_samples(const char *expr, int *n, int n_pre, fmf_t *const*fmf)
{
	int err, i, j;
	khint_t k;
	char **s;
	kexpr_t *ke;
	khash_t(s2i) *h;

	*n = 0;
	if (*expr != '?') return hts_readlines(expr, n);
	ke = ke_parse(expr+1, &err);
	if (err) return 0;
	h = kh_init(s2i);
	for (j = 0; j < n_pre; ++j) {
		fmf_t *f = fmf[j];
		int absent;
		if (f == 0) continue;
		for (i = 0; i < f->n_rows; ++i)
			if (fmf_test(f, i, ke)) {
				k = kh_put(s2i, h, f->rows[i].name, &absent);
				if (absent) kh_key(h, k) = strdup(f->rows[i].name);
			}
	}
	ke_destroy(ke);
	*n = kh_size(h);
	s = (char**)malloc(*n * sizeof(char*));
	for (k = i = 0; k < kh_end(h); ++k)
		if (kh_exist(h, k)) s[i++] = (char*)kh_key(h, k);
	kh_destroy(s2i, h);
	return s;
}

int main_view(int argc, char *argv[])
{
	int i, c, out_bcf = 0, clevel = -1, is_multi = 0, multi_flag = 0, excl = 0;
	bgt_t *bgt = 0;
	bgtm_t *bm = 0;
	bcf1_t *b;
	htsFile *out;
	char modew[8], *reg = 0, *sexpr = 0;
	flt_aux_t flt;
	void *bed = 0;
	int n_groups = 0;
	char *gexpr[BGT_MAX_GROUPS];

	memset(&flt, 0, sizeof(flt_aux_t));
	assert(strcmp(argv[0], "view") == 0 || strcmp(argv[0], "sview") == 0);
	is_multi = strcmp(argv[0], "sview")? 1 : 0;
	while ((c = getopt(argc, argv, "bs:r:l:aGB:ef:g:")) >= 0) {
		if (c == 'b') out_bcf = 1;
		else if (c == 'r') reg = optarg;
		else if (c == 'l') clevel = atoi(optarg);
		else if (c == 'e') excl = 1;
		else if (c == 'B') bed = bed_read(optarg);
		else if (is_multi && c == 's') sexpr = optarg;
		else if (is_multi && c == 'a') multi_flag |= BGT_F_SET_AC;
		else if (is_multi && c == 'G') multi_flag |= BGT_F_NO_GT;
		else if (is_multi && c == 'f') {
			int err = 0;
			flt.ke = ke_parse(optarg, &err);
			assert(err == 0 && flt.ke != 0);
		} else if (c == 'g' && n_groups < BGT_MAX_GROUPS) gexpr[n_groups++] = optarg;
	}
	if (clevel > 9) clevel = 9;
	if (argc - optind < 1) {
		fprintf(stderr, "Usage: bgt %s [options] <bgt-prefix>", argv[0]);
		if (is_multi) fprintf(stderr, " [...]");
		fputc('\n', stderr);
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b           BCF output\n");
		fprintf(stderr, "  -r STR       region [all]\n");
		fprintf(stderr, "  -l INT       compression level for BCF [detault]\n");
		fprintf(stderr, "  -B FILE      extract variants overlapping BED FILE [null]\n");
		fprintf(stderr, "  -e           exclude variants overlapping BED FILE (effective with -B)\n");
		if (is_multi) {
			fprintf(stderr, "  -a           write AC/AN to the INFO field\n");
			fprintf(stderr, "  -G           don't output sample genotype\n");
			fprintf(stderr, "  -f STR       frequency filters [null]\n");
			fprintf(stderr, "  -s EXPR      list of samples (see Notes below) [all]\n");
			fprintf(stderr, "  -g EXPR      define a sample group (see Notes below) [null]\n");
		}
		return 1;
	}

	if (!is_multi) {
		bgt = bgt_open(argv[optind]);
		if (reg) bgt_set_region(bgt, reg);
		if (bed) bgt_set_bed(bgt, bed, excl);
	} else {
		fmf_t **fmf = 0;
		int n_fmf = argc - optind, n_samples = 0;
		char **samples = 0;
		bm = bgtm_open(argc - optind, &argv[optind]);
		bgtm_set_flag(bm, multi_flag);
		if (flt.ke) bgtm_set_filter(bm, filter_func, &flt);
		if (reg) bgtm_set_region(bm, reg);
		if (bed) bgtm_set_bed(bm, bed, excl);

		fmf = (fmf_t**)malloc(n_fmf * sizeof(fmf_t*));
		for (i = 0; i < n_fmf; ++i) {
			char *tmpfn;
			tmpfn = (char*)calloc(strlen(argv[optind+i]) + 5, 1);
			strcat(strcpy(tmpfn, argv[optind+i]), ".spl");
			fmf[i] = fmf_read(tmpfn);
			free(tmpfn);
		}
		if (sexpr) {
			samples = get_samples(sexpr, &n_samples, n_fmf, fmf);
			bgtm_set_samples(bm, n_samples, samples);
			for (i = 0; i < n_samples; ++i) free(samples[i]);
			free(samples);
		}
		for (i = 0; i < n_groups; ++i) {
			int j;
			samples = get_samples(gexpr[i], &n_samples, n_fmf, fmf);
			bgtm_add_group(bm, n_samples, samples);
			for (j = 0; j < n_samples; ++j) free(samples[j]);
			free(samples);
		}
		for (i = 0; i < n_fmf; ++i) fmf_destroy(fmf[i]);
		free(fmf);
	}

	strcpy(modew, "w");
	if (out_bcf) strcat(modew, "b");
	sprintf(modew + strlen(modew), "%d", clevel);
	out = hts_open("-", modew, 0);
	vcf_hdr_write(out, bgt? bgt->h_out : bm->h_out);

	b = bcf_init1();
	while (read1(bgt, bm, b) >= 0)
		vcf_write1(out, bgt? bgt->h_out : bm->h_out, b);
	bcf_destroy1(b);

	hts_close(out);
	if (bgt) bgt_close(bgt);
	if (bm) bgtm_close(bm);
	if (bed) bed_destroy(bed);
	if (flt.ke) ke_destroy(flt.ke);
	return 0;
}
