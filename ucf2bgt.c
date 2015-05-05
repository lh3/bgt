#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>
#include "vcf.h"

int main_ucf2bgt(int argc, char *argv[])
{
	int c, clevel = -1, flag = 0, id_GT = -1;
	char *fn_ref = 0, *fn_out = 0, moder[8], modew[8];
	uint8_t *bits[2];
	htsFile *in, *out;
	bcf_hdr_t *h;
	bcf1_t *b;

	while ((c = getopt(argc, argv, "l:bSt:o:")) >= 0) {
		switch (c) {
		case 'l': clevel = atoi(optarg); flag |= 2; break;
		case 'S': flag |= 1; break;
		case 'b': flag |= 2; break;
		case 't': fn_ref = optarg; flag |= 1; break;
		case 'o': fn_out = optarg; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: bgt ucf2bgt [options] <in.bcf>|<in.vcf>|<in.vcf.gz> <out-prefix>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b           output in BCF\n");
		fprintf(stderr, "  -S           input is VCF\n");
		fprintf(stderr, "  -o FILE      output file name [stdout]\n");
		fprintf(stderr, "  -l INT       compression level [%d]\n", clevel);
		fprintf(stderr, "  -t FILE      list of reference names and lengths [null]\n");
		return 1;
	}
	strcpy(moder, "r");
	if ((flag&1) == 0) strcat(moder, "b");

	in = hts_open(argv[optind], moder, fn_ref);
	assert(in);
	h = vcf_hdr_read(in);
	assert(h);
	assert(h->n[BCF_DT_SAMPLE] > 0);
	id_GT = bcf_id2int(h, BCF_DT_ID, "GT");
	assert(id_GT >= 0);
	bits[0] = (uint8_t*)calloc(h->n[BCF_DT_SAMPLE]*2, 1);
	bits[1] = (uint8_t*)calloc(h->n[BCF_DT_SAMPLE]*2, 1);

	strcpy(modew, "w");
	if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
	if (flag&2) strcat(modew, "b");
	out = hts_open(fn_out? fn_out : "-", modew, 0);
	vcf_hdr_write(out, h);
	b = bcf_init1();
	while (vcf_read1(in, h, b) >= 0) {
		int i, k, j;
		bcf_fmt_t *gt;
		bcf_unpack(b, BCF_UN_FMT);
		for (i = 0; i < b->n_fmt; ++i)
			if (b->d.fmt[i].id == id_GT) break;
		if (i == b->n_fmt) continue; // no GT field
		gt = &b->d.fmt[i];
		assert(gt->type == BCF_BT_INT8);
		for (i = k = 0; i < b->n_sample; ++i, k += gt->n) {
			for (j = 0; j < gt->n; ++j) {
				int a = (int)(gt->p[k+j] >> 1) - 1;
				if (a < 0) bits[0][j+k] = 0, bits[1][j+k] = 1;
				else if (a >= 2) bits[0][j+k] = bits[1][j+k] = 1;
				else bits[0][j+k] = a, bits[1][j+k] = 0;
			}
		}
		for (i = 0; i < b->n_sample*gt->n; ++i) putchar("01"[bits[0][i]]); putchar('\n');
		for (i = 0; i < b->n_sample*gt->n; ++i) putchar("01"[bits[1][i]]); putchar('\n');
		vcf_write1(out, h, b);
	}
	bcf_destroy1(b);
	hts_close(out);

	free(bits[0]); free(bits[1]);
	bcf_hdr_destroy(h);
	hts_close(in);
	return 0;
}
