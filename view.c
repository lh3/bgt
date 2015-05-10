#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include "bgt.h"

typedef struct {
	int min_ac, min_an;
	float min_af;
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
	if (an < flt->min_an) return 1;
	if (ac < flt->min_ac) return 1;
	if (an > 0 && (float)ac/an < flt->min_af) return 1;
	return 0;
}

int main_view(int argc, char *argv[])
{
	int i, c, out_bcf = 0, n_samples = 0, clevel = -1, is_multi = 0, multi_flag = 0, set_flt = 0;
	bgt_t *bgt = 0;
	bgtm_t *bm = 0;
	bcf1_t *b;
	htsFile *out;
	char modew[8], *reg = 0, **samples = 0;
	flt_aux_t flt;

	memset(&flt, 0, sizeof(flt_aux_t));
	assert(strcmp(argv[0], "view") == 0 || strcmp(argv[0], "mview") == 0);
	is_multi = strcmp(argv[0], "mview") == 0? 1 : 0;
	while ((c = getopt(argc, argv, "bs:r:l:aGC:F:N:")) >= 0) {
		if (c == 'b') out_bcf = 1;
		else if (c == 'r') reg = optarg;
		else if (c == 'l') clevel = atoi(optarg);
		else if (c == 's') samples = hts_readlines(optarg, &n_samples);
		else if (is_multi && c == 'a') multi_flag |= BGT_F_SET_AC;
		else if (is_multi && c == 'G') multi_flag |= BGT_F_NO_GT;
		else if (is_multi && c == 'N') flt.min_an = atoi(optarg), set_flt = 1;
		else if (is_multi && c == 'C') flt.min_ac = atoi(optarg), set_flt = 1;
		else if (is_multi && c == 'F') flt.min_af = atof(optarg), set_flt = 1;
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
		if (is_multi) {
			fprintf(stderr, "  -a           write AC/AN to the INFO field\n");
			fprintf(stderr, "  -G           don't output sample genotype\n");
			fprintf(stderr, "  -N INT       min AN [0]\n");
			fprintf(stderr, "  -C INT       min AC [0]\n");
			fprintf(stderr, "  -F FLOAT     min AF [0]\n");
		}
		return 1;
	}

	if (!is_multi) {
		bgt = bgt_open(argv[optind]);
		if (n_samples > 0) bgt_set_samples(bgt, n_samples, samples);
		if (reg) bgt_set_region(bgt, reg);
	} else {
		bm = bgtm_open(argc - optind, &argv[optind]);
		bgtm_set_flag(bm, multi_flag);
		if (set_flt) bgtm_set_filter(bm, filter_func, &flt);
		if (n_samples > 0) bgtm_set_samples(bm, n_samples, samples);
		if (reg) bgtm_set_region(bm, reg);
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

	hts_close(out);
	if (bgt) bgt_close(bgt);
	if (bm) bgtm_close(bm);
	return 0;
}
