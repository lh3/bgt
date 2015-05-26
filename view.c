#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include "bgt.h"
#include "kexpr.h"
#include "fmf.h"

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

static int filter_func(bcf_hdr_t *h, bcf1_t *b, int an, int ac, int n, const int *sidx, uint8_t *a[2], void *data)
{
	flt_aux_t *flt = (flt_aux_t*)data;
	int is_true, err;
	if (flt->ke == 0) return 0;
	ke_set_real(flt->ke, "AN", an);
	ke_set_real(flt->ke, "AC", ac);
	is_true = !!ke_eval_int(flt->ke, &err);
	if (err) return 0;
	return !is_true;
}

static char **get_samples(const char *expr, const fmf_t *fmf, int *_n)
{
	int err, i, n, m;
	char **s;
	kexpr_t *ke;
	*_n = 0;
	if (*expr != '?') return hts_readlines(expr, _n);
	ke = ke_parse(expr+1, &err);
	if (err) return 0;
	for (i = n = m = 0, s = 0; i < fmf->n_rows; ++i) {
		if (fmf_test(fmf, i, ke)) {
			if (n == m) {
				m = m? m<<1 : 16;
				s = (char**)realloc(s, m * sizeof(char*));
			}
			s[n++] = strdup(fmf->rows[i].name);
		}
	}
	s = (char**)realloc(s, n * sizeof(char*));
	*_n = n;
	return s;
}

int main_view(int argc, char *argv[])
{
	int i, c, out_bcf = 0, n_samples = 0, clevel = -1, is_multi = 0, multi_flag = 0, excl = 0;
	bgt_t *bgt = 0;
	bgtm_t *bm = 0;
	bcf1_t *b;
	htsFile *out;
	fmf_t *fmf = 0;
	char modew[8], *reg = 0, *sexpr, **samples = 0, *tmpfn;
	flt_aux_t flt;
	void *bed = 0;

	memset(&flt, 0, sizeof(flt_aux_t));
	assert(strcmp(argv[0], "view") == 0 || strcmp(argv[0], "mview") == 0);
	is_multi = strcmp(argv[0], "mview") == 0? 1 : 0;
	while ((c = getopt(argc, argv, "bs:r:l:aGB:ef:")) >= 0) {
		if (c == 'b') out_bcf = 1;
		else if (c == 'r') reg = optarg;
		else if (c == 'l') clevel = atoi(optarg);
		else if (c == 's') sexpr = optarg;
		else if (c == 'e') excl = 1;
		else if (c == 'B') bed = bed_read(optarg);
		else if (is_multi && c == 'a') multi_flag |= BGT_F_SET_AC;
		else if (is_multi && c == 'G') multi_flag |= BGT_F_NO_GT;
		else if (is_multi && c == 'f') {
			int err = 0;
			flt.ke = ke_parse(optarg, &err);
			assert(err == 0 && flt.ke != 0);
		}
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
		fprintf(stderr, "  -s STR/FILE  list of samples (STR if started with ':'; FILE otherwise) [all]\n");
		fprintf(stderr, "  -B FILE      extract variants overlapping BED FILE [null]\n");
		fprintf(stderr, "  -e           exclude variants overlapping BED FILE (effective with -B) [null]\n");
		if (is_multi) {
			fprintf(stderr, "  -a           write AC/AN to the INFO field\n");
			fprintf(stderr, "  -G           don't output sample genotype\n");
			fprintf(stderr, "  -f STR       frequency filters [null]\n");
		}
		return 1;
	}

	tmpfn = (char*)malloc(strlen(argv[optind]) + 5);
	sprintf(tmpfn, "%s.spl", argv[optind]);
	fmf = fmf_read(tmpfn);
	free(tmpfn);

	if (sexpr) samples = get_samples(sexpr, fmf, &n_samples);

	if (!is_multi) {
		bgt = bgt_open(argv[optind]);
		if (n_samples > 0) bgt_set_samples(bgt, n_samples, samples);
		if (reg) bgt_set_region(bgt, reg);
		if (bed) bgt_set_bed(bgt, bed, excl);
	} else {
		bm = bgtm_open(argc - optind, &argv[optind]);
		bgtm_set_flag(bm, multi_flag);
		if (flt.ke) bgtm_set_filter(bm, filter_func, &flt);
		if (n_samples > 0) bgtm_set_samples(bm, n_samples, samples);
		if (reg) bgtm_set_region(bm, reg);
		if (bed) bgtm_set_bed(bm, bed, excl);
	}
	for (i = 0; i < n_samples; ++i) free(samples[i]);
	free(samples);

	strcpy(modew, "w");
	if (out_bcf) strcat(modew, "b");
	sprintf(modew + strlen(modew), "%d", clevel);
	out = hts_open("-", modew, 0);
	vcf_hdr_write(out, bgt? bgt->h_out : bm->h_out);

	b = bcf_init1();
	while (read1(bgt, bm, b) >= 0)
		vcf_write1(out, bgt? bgt->h_out : bm->h_out, b);
	bcf_destroy1(b);

	fmf_destroy(fmf);
	hts_close(out);
	if (bgt) bgt_close(bgt);
	if (bm) bgtm_close(bm);
	if (bed) bed_destroy(bed);
	if (flt.ke) ke_destroy(flt.ke);
	return 0;
}
