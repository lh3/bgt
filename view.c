#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include "bgt.h"

int main_view(int argc, char *argv[])
{
	int i, c, out_bcf = 0, n_samples = 0, clevel = -1;
	bgt_t *bgt;
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
		fprintf(stderr, "Usage: bgt view [options] <bgt-prefix>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b           BCF output\n");
		fprintf(stderr, "  -r STR       region [all]\n");
		fprintf(stderr, "  -l INT       compression level for BCF [detault]\n");
		fprintf(stderr, "  -s STR/FILE  list of samples (STR if started with ':'; FILE otherwise) [all]\n");
		return 1;
	}

	bgt = bgt_open(argv[optind]);
	if (n_samples > 0) bgt_set_samples(bgt, n_samples, samples);
	if (reg) bgt_set_region(bgt, reg);

	strcpy(modew, "w");
	if (out_bcf) strcat(modew, "b");
	sprintf(modew + strlen(modew), "%d", clevel);
	out = hts_open("-", modew, 0);
	vcf_hdr_write(out, bgt->h_sub);

	b = bcf_init1();
	while (bgt_read(bgt, b) >= 0) {
		vcf_write1(out, bgt->h_sub, b);
	}
	bcf_destroy1(b);

	hts_close(out);
	bgt_close(bgt);

	for (i = 0; i < n_samples; ++i) free(samples[i]);
	free(samples);
	return 0;
}
