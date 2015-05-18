#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>
#include "vcf.h"
#include "pbwt.h"

int main_ucf2bgt(int argc, char *argv[])
{
	int i, c, clevel = -1, flag = 0, id_GT = -1;
	char *fn_ref = 0, moder[8], modew[8];
	char *prefix, *fn;
	uint8_t *bits[2];
	int64_t n = 0;
	htsFile *in, *out;
	bcf_hdr_t *h, *h0;
	bcf1_t *b;
	FILE *fp;
	pbf_t *pb;

	while ((c = getopt(argc, argv, "l:St:")) >= 0) {
		switch (c) {
		case 'l': clevel = atoi(optarg); flag |= 2; break;
		case 'S': flag |= 1; break;
		case 't': fn_ref = optarg; flag |= 1; break;
		}
	}
	if (argc - optind < 2) {
		fprintf(stderr, "Usage: bgt ucf2bgt [options] <in.bcf>|<in.vcf>|<in.vcf.gz> <out-prefix>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -S           input is VCF\n");
		fprintf(stderr, "  -l INT       compression level [%d]\n", clevel);
		fprintf(stderr, "  -t FILE      list of reference names and lengths [null]\n");
		return 1;
	}
	prefix = argv[optind+1];
	fn = (char*)malloc(strlen(prefix) + 9);
	strcpy(moder, "r");
	if ((flag&1) == 0) strcat(moder, "b");

	in = hts_open(argv[optind], moder, fn_ref);
	assert(in);
	h = vcf_hdr_read(in);
	assert(h);
	assert(h->n[BCF_DT_SAMPLE] > 0);
	id_GT = bcf_id2int(h, BCF_DT_ID, "GT");
	if (id_GT < 0) {
		bcf_hdr_append(h, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		id_GT = bcf_id2int(h, BCF_DT_ID, "GT");
	}
	bcf_hdr_append(h, "##INFO=<ID=_row,Number=1,Type=Integer,Description=\"row number\">");
	h0 = bcf_hdr_subset(h, 0, 0, 0);

	// write sample list
	sprintf(fn, "%s.spl", prefix);
	fp = fopen(fn, "wb");
	for (i = 0; i < h->n[BCF_DT_SAMPLE]; ++i) {
		fputs(h->id[BCF_DT_SAMPLE][i].key, fp);
		fputc('\n', fp);
	}
	fclose(fp);

	// prepare PBF to write
	sprintf(fn, "%s.pbf", prefix);
	pb = pbf_open_w(fn, h->n[BCF_DT_SAMPLE]*2, 2, 13);
	bits[0] = (uint8_t*)calloc(h->n[BCF_DT_SAMPLE]*2, 1);
	bits[1] = (uint8_t*)calloc(h->n[BCF_DT_SAMPLE]*2, 1);

	strcpy(modew, "wb");
	if (clevel >= 0 && clevel <= 9) sprintf(modew + 2, "%d", clevel);
	sprintf(fn, "%s.bcf", prefix);
	out = hts_open(fn, modew, 0);
	vcf_hdr_write(out, h0);
	b = bcf_init1();

	bcf_atom_v atoms = {0,0,0};

	while (vcf_read1(in, h, b) >= 0) {
		int i, k, j;
		int32_t val = n;
		bcf_fmt_t *gt;

		// insert "_row" to INFO
		bcf_append_info_ints(h, b, "_row", 1, &val);

		atoms.n = 0;
		bcf_atomize(h, b, &atoms);
		for (i = 0; i < atoms.n; ++i)
			fprintf(stderr, "%d\t%d\t%s\t%s\n", atoms.a[i].rid, atoms.a[i].pos, atoms.a[i].ref.s, atoms.a[i].alt);

		// write genotypes
		bcf_unpack(b, BCF_UN_FMT);
		for (i = 0; i < b->n_fmt; ++i)
			if (b->d.fmt[i].id == id_GT) break;
		if (i == b->n_fmt) continue; // no GT field
		gt = &b->d.fmt[i];
		assert(gt->type == BCF_BT_INT8);
		for (i = k = 0; i < b->n_sample; ++i, k += 2) {
			assert(gt->n <= 2);
			for (j = 0; j < gt->n; ++j) {
				int a = (int)(gt->p[k+j] >> 1) - 1;
				int b = a < 0? 2 : a >= 2? 3 : a;
				if (gt->n == 2) bits[0][k+j] = b&1, bits[1][k+j] = b>>1;
				else bits[0][k] = bits[0][k+1] = b&1, bits[1][k] = bits[1][k+1] = b>>1;
			}
		}
		pbf_write(pb, bits);

		// write BCF
		bcf_subset(h0, b, 0, 0);
		vcf_write1(out, h0, b);
		++n;
	}
	bcf_destroy1(b);
	hts_close(out);

	pbf_close(pb);
	free(bits[0]); free(bits[1]);

	bcf_hdr_destroy(h0);
	bcf_hdr_destroy(h);
	hts_close(in);

	bcf_index_build(fn, 14);
	free(fn);
	return 0;
}
