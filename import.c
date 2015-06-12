#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>
#include "atomic.h"
#include "pbwt.h"

int main_import(int argc, char *argv[])
{
	int i, c, clevel = -1, flag = 0, id_GT = -1;
	char *fn_ref = 0, moder[8], modew[8];
	char *prefix, *fn;
	uint8_t *bits[2], *bit1;
	int64_t n = 0;
	htsFile *in, *out;
	bcf_hdr_t *h0;
	bcf1_t *b;
	FILE *fp;
	pbf_t *pb, *pb1;
	bcf_atombuf_t *ab;
	const bcf_atom_t *a;

	while ((c = getopt(argc, argv, "l:SFt:")) >= 0) {
		switch (c) {
		case 'l': clevel = atoi(optarg); flag |= 2; break;
		case 'S': flag |= 1; break;
		case 't': fn_ref = optarg; flag |= 1; break;
		case 'F': flag |= 4; break;
		}
	}
	if (argc - optind < 2) {
		fprintf(stderr, "Usage: bgt import [options] <out-prefix> <in.bcf>|<in.vcf>|<in.vcf.gz>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -S           input is VCF\n");
		fprintf(stderr, "  -t FILE      list of reference names and lengths [null]\n");
		fprintf(stderr, "  -F           keep filtered variants\n");
		return 1;
	}
	prefix = argv[optind];
	fn = (char*)malloc(strlen(prefix) + 9);
	strcpy(moder, "r");
	if ((flag&1) == 0) strcat(moder, "b");

	in = hts_open(argv[optind+1], moder, fn_ref);
	assert(in);
	ab = bcf_atombuf_init(in, flag&4);
	assert(ab->h->n[BCF_DT_SAMPLE] > 0);
	h0 = bcf_hdr_subset(ab->h, 0, 0, 0);
	id_GT = bcf_id2int(h0, BCF_DT_ID, "GT");
	if (id_GT < 0) {
		bcf_hdr_append(h0, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		id_GT = bcf_id2int(h0, BCF_DT_ID, "GT");
	}
	bcf_hdr_append(h0, "##INFO=<ID=_row,Number=1,Type=Integer,Description=\"row number\">");

	// write sample list
	sprintf(fn, "%s.spl", prefix);
	fp = fopen(fn, "wb");
	for (i = 0; i < ab->h->n[BCF_DT_SAMPLE]; ++i) {
		fputs(ab->h->id[BCF_DT_SAMPLE][i].key, fp);
		fputc('\n', fp);
	}
	fclose(fp);

	// prepare PBF to write
	sprintf(fn, "%s.pbf", prefix);
	pb = pbf_open_w(fn, ab->h->n[BCF_DT_SAMPLE]*2, 2, 13);
	bits[0] = (uint8_t*)calloc(ab->h->n[BCF_DT_SAMPLE]*2, 1);
	bits[1] = (uint8_t*)calloc(ab->h->n[BCF_DT_SAMPLE]*2, 1);

	// prepare PBF to write
	sprintf(fn, "%s.pb1", prefix);
	pb1 = pbf_open_w(fn, ab->h->n[BCF_DT_SAMPLE]*2, 1, 13);
	bit1 = (uint8_t*)calloc(ab->h->n[BCF_DT_SAMPLE]*2, 1);

	strcpy(modew, "wb");
	if (clevel >= 0 && clevel <= 9) sprintf(modew + 2, "%d", clevel);
	sprintf(fn, "%s.bcf", prefix);
	out = hts_open(fn, modew, 0);
	vcf_hdr_write(out, h0);
	b = bcf_init1();

	while ((a = bcf_atom_read(ab)) != 0) {
		int32_t i, val = n;
		bcf_atom2bcf(a, b, 1, -1);
		bcf_append_info_ints(h0, b, "_row", 1, &val);
		for (i = 0; i < a->n_gt; ++i) {
			bits[0][i] = a->gt[i]&1, bits[1][i] = a->gt[i]>>1&1;
			bit1[i] = (a->gt[i] == 1);
		}
		pbf_write(pb, bits);
		pbf_write(pb1, &bit1);
		bcf_subset(h0, b, 0, 0);
		vcf_write1(out, h0, b);
		++n;
	}

	bcf_destroy1(b);
	hts_close(out);

	pbf_close(pb1); free(bit1);
	pbf_close(pb);  free(bits[0]); free(bits[1]);

	bcf_hdr_destroy(h0);
	bcf_atombuf_destroy(ab);
	hts_close(in);

	bcf_index_build(fn, 14);
	free(fn);
	return 0;
}

int main_bcfidx(int argc, char *argv[])
{
	int c, min_shift = 14;
	while ((c = getopt(argc, argv, "s:")) >= 0)
		if (c == 's') min_shift = atoi(optarg);
	if (optind == argc) {
		fprintf(stderr, "Usage: bgt bcfidx [-s minShift] <in.bcf>\n");
		return 1;
	}
	bcf_index_build(argv[optind], min_shift);
	return 0;
}
