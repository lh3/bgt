#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include "bgt.h"

static inline int read1(bgt_t *bgt, bgtm_t *bm, bcf1_t *b)
{
	if (bgt) return bgt_read(bgt, b);
	if (bm) return bgtm_read(bm, b);
	return -1;
}

int main_view(int argc, char *argv[])
{
	int i, c, out_bcf = 0, n_samples = 0, clevel = -1;
	bgt_t *bgt = 0;
	bgtm_t *bm = 0;
	bcf1_t *b;
	htsFile *out;
	char modew[8], *reg = 0, **samples = 0;

	while ((c = getopt(argc, argv, "bs:r:l:")) >= 0) {
		if (c == 'b') out_bcf = 1;
		else if (c == 'r') reg = optarg;
		else if (c == 'l') clevel = atoi(optarg);
		else if (c == 's') samples = hts_readlines(optarg, &n_samples);
	}
	if (clevel > 9) clevel = 9;
	if (argc - optind < 1) {
		fprintf(stderr, "Usage: bgt %s [options] <bgt-prefix>\n", argv[0]);
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b           BCF output\n");
		fprintf(stderr, "  -r STR       region [all]\n");
		fprintf(stderr, "  -l INT       compression level for BCF [detault]\n");
		fprintf(stderr, "  -s STR/FILE  list of samples (STR if started with ':'; FILE otherwise) [all]\n");
		return 1;
	}

	if (strcmp(argv[0], "view") == 0) {
		bgt = bgt_open(argv[optind]);
		if (n_samples > 0) bgt_set_samples(bgt, n_samples, samples);
		if (reg) bgt_set_region(bgt, reg);
	} else if (strcmp(argv[0], "mview") == 0) {
		bm = bgtm_open(argc - optind, &argv[optind]);
		if (n_samples > 0) bgtm_set_samples(bm, n_samples, samples);
		if (reg) bgtm_set_region(bm, reg);
	} else abort();
	for (i = 0; i < n_samples; ++i) free(samples[i]);
	free(samples);

	strcpy(modew, "w");
	if (out_bcf) strcat(modew, "b");
	sprintf(modew + strlen(modew), "%d", clevel);
	out = hts_open("-", modew, 0);
	vcf_hdr_write(out, bgt? bgt->h_sub : bm->h);

	b = bcf_init1();
	while (read1(bgt, bm, b) >= 0) {
		vcf_write1(out, bgt? bgt->h_sub : bm->h, b);
	}
	bcf_destroy1(b);

	hts_close(out);
	if (bgt) bgt_close(bgt);
	if (bm) bgtm_close(bm);
	return 0;
}
