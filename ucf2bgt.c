#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include "vcf.h"

int main_ucf2bgt(int argc, char *argv[])
{
	int c, clevel = -1, flag = 0;
	char *fn_ref = 0, *fn_out = 0, moder[8], modew[8];
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
	h = vcf_hdr_read(in);
	b = bcf_init1();

	strcpy(modew, "w");
	if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
	if (flag&2) strcat(modew, "b");
	out = hts_open(fn_out? fn_out : "-", modew, 0);
	vcf_hdr_write(out, h);
	while (vcf_read1(in, h, b) >= 0) {
		vcf_write1(out, h, b);
	}
	hts_close(out);

	bcf_destroy1(b);
	bcf_hdr_destroy(h);
	hts_close(in);
	return 0;
}
