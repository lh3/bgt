#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include "bgt.h"
#include "kexpr.h"
#include "fmf.h"

#include "khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)

void *bed_read(const char *fn);
void bed_destroy(void *_h);
char **hts_readlines(const char *fn, int *_n);

typedef struct {
	kexpr_t *ke;
} flt_aux_t;

static int filter_func(bcf_hdr_t *h, bcf1_t *b, int an, int ac1, int n_groups, int32_t *gan, int32_t *gac1, void *data)
{
	flt_aux_t *flt = (flt_aux_t*)data;
	int is_true, err, i;
	char key[4];
	if (flt->ke == 0) return 0;
	ke_set_int(flt->ke, "AN", an);
	ke_set_int(flt->ke, "AC", ac1);
	key[3] = 0;
	for (i = 1; i <= n_groups; ++i) {
		key[0] = 'A', key[1] = 'N', key[2] = '0' + i; ke_set_int(flt->ke, key, gan[i]);
		key[0] = 'A', key[1] = 'C', key[2] = '0' + i; ke_set_int(flt->ke, key, gac1[i]);
	}
	is_true = !!ke_eval_int(flt->ke, &err);
	if (err) return 0;
	return !is_true;
}

int main_view(int argc, char *argv[])
{
	int i, c, n_files = 0, out_bcf = 0, clevel = -1, multi_flag = 0, excl = 0;
	long seekn = -1, n_rec = LONG_MAX, n_read = 0;
	bgtm_t *bm = 0;
	bcf1_t *b;
	htsFile *out;
	char modew[8], *reg = 0;
	flt_aux_t flt;
	void *bed = 0;
	int n_groups = 0;
	char *gexpr[BGT_MAX_GROUPS];
	int n_a = 0, m_a = 0;
	bgt_allele_t *a = 0;
	bgt_file_t **files = 0;

	memset(&flt, 0, sizeof(flt_aux_t));
	while ((c = getopt(argc, argv, "bs:r:l:AGB:ef:g:a:i:n:")) >= 0) {
		if (c == 'b') out_bcf = 1;
		else if (c == 'r') reg = optarg;
		else if (c == 'l') clevel = atoi(optarg);
		else if (c == 'e') excl = 1;
		else if (c == 'B') bed = bed_read(optarg);
		else if (c == 'A') multi_flag |= BGT_F_SET_AC;
		else if (c == 'G') multi_flag |= BGT_F_NO_GT;
		else if (c == 'i') seekn = atol(optarg) - 1;
		else if (c == 'n') n_rec = atol(optarg);
		else if (c == 'f') {
			int err = 0;
			flt.ke = ke_parse(optarg, &err);
			assert(err == 0 && flt.ke != 0);
		} else if (c == 's' && n_groups < BGT_MAX_GROUPS) {
			gexpr[n_groups++] = optarg;
		} else if (c == 'a') {
			bgt_allele_t *p;
			if (n_a == m_a) {
				int old = m_a;
				m_a = m_a? m_a<<1 : 2;
				a = (bgt_allele_t*)realloc(a, m_a * sizeof(bgt_allele_t));
				memset(a + old, 0, (m_a - old) * sizeof(bgt_allele_t));
			}
			p = &a[n_a++];
			if (bgt_al_parse(optarg, p) < 0) {
				fprintf(stderr, "[E::%s] failed to parse allele '%s'\n", __func__, optarg);
				return 1; // FIXME: memory leak
			}
			if (n_a > 1 && strcmp(p->chr.s, a[0].chr.s) != 0) {
				fprintf(stderr, "[E::%s] for now, -a doesn't support alleles on different CHROM. Abort!\n", __func__);
				return 1; // FIXME: memory leak
			}
		}
	}
	if (seekn < 0) seekn = 0;
	if (clevel > 9) clevel = 9;
	if (n_groups > 1) multi_flag |= BGT_F_SET_AC;
	if (reg && n_a > 0) {
		fprintf(stderr, "[W::%s] -r and -a can't be specified at the same time. -r is igored.\n", __func__);
		reg = 0;
	}
	if (argc - optind < 1) {
		fprintf(stderr, "Usage: bgt %s [options] <bgt-prefix> [...]", argv[0]);
		fputc('\n', stderr);
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b           BCF output\n");
		fprintf(stderr, "  -r STR       region [all]\n");
		fprintf(stderr, "  -i INT       process from the INT-th record (1-based) [null]\n");
		fprintf(stderr, "  -n INT       process at most INT records [null]\n");
		fprintf(stderr, "  -l INT       compression level for BCF [default]\n");
		fprintf(stderr, "  -B FILE      extract variants overlapping BED FILE [null]\n");
		fprintf(stderr, "  -e           exclude variants overlapping BED FILE (effective with -B)\n");
		fprintf(stderr, "  -a           write AC/AN to the INFO field\n");
		fprintf(stderr, "  -G           don't output sample genotype\n");
		fprintf(stderr, "  -f STR       frequency filters [null]\n");
		fprintf(stderr, "  -s EXPR      list of samples (see Notes below) [all]\n");
		fprintf(stderr, "Notes:\n");
		fprintf(stderr, "  For option -s, EXPR can be one of:\n");
		fprintf(stderr, "    1) comma-delimited sample list following a colon. e.g. -s:NA12878,NA12044\n");
		fprintf(stderr, "    2) space-delimited file with the first column giving a sample name. e.g. -s samples.txt\n");
		fprintf(stderr, "    3) expression if .spl file contains metadata. e.g.: -s\"gender=='M'&&population!='CEU'\"\n");
		fprintf(stderr, "  If multiple -s is specified, the AC/AN of the first group will be written to VCF INFO AC1/AN1,\n");
		fprintf(stderr, "  the second to AC2/AN2, etc.\n");
		return 1;
	}

	// set reg
	if (n_a > 0) {
		int min = INT_MAX, max = INT_MIN;
		for (i = 0; i < n_a; ++i) {
			min = min < a[i].pos? min : a[i].pos;
			max = max > a[i].pos? max : a[i].pos;
		}
		reg = (char*)malloc(strlen(a[0].chr.s) + 24);
		sprintf(reg, "%s:%d-%d", a[0].chr.s, min+1, max+1);
	}

	n_files = argc - optind;
	files = (bgt_file_t**)calloc(n_files, sizeof(bgt_file_t*));
	for (i = 0; i < n_files; ++i) files[i] = bgt_open(argv[optind+i]);

	bm = bgtm_reader_init(n_files, files);
	bgtm_set_flag(bm, multi_flag);
	if (flt.ke) bgtm_set_filter(bm, filter_func, &flt);
	if (reg) bgtm_set_region(bm, reg);
	if (bed) bgtm_set_bed(bm, bed, excl);
	if (seekn >= 0) bgtm_set_start(bm, seekn);
	for (i = 0; i < n_groups; ++i)
		bgtm_add_group(bm, gexpr[i]);
	bgtm_prepare(bm); // bgtm_prepare() generates the VCF header

	strcpy(modew, "w");
	if (out_bcf) strcat(modew, "b");
	sprintf(modew + strlen(modew), "%d", clevel);
	out = hts_open("-", modew, 0);
	vcf_hdr_write(out, bm->h_out);

	b = bcf_init1();
	while (bgtm_read(bm, b) >= 0 && n_read < n_rec) {
		vcf_write1(out, bm->h_out, b);
		++n_read;
	}
	bcf_destroy1(b);

	hts_close(out);
	bgtm_reader_destroy(bm);
	if (bed) bed_destroy(bed);
	if (flt.ke) ke_destroy(flt.ke);
	if (n_a) free(reg);
	for (i = 0; i < n_a; ++i) free(a[i].chr.s);
	free(a);
	for (i = 0; i < n_files; ++i) bgt_close(files[i]);
	free(files);
	return 0;
}

int main_getalt(int argc, char *argv[])
{
	int c;
	char *fn;
	BGZF *fp;
	bcf1_t *b;
	bcf_hdr_t *h;
	kstring_t s = {0,0,0};

	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (argc - optind == 0) {
		fprintf(stderr, "Usage: bgt getalt <bgt-base>\n");
		return 1;
	}

	fn = (char*)calloc(strlen(argv[optind]) + 5, 1);
	sprintf(fn, "%s.bcf", argv[optind]);
	fp = bgzf_open(fn, "r");
	free(fn);
	assert(fp);

	h = bcf_hdr_read(fp);
	b = bcf_init1();
	while (bcf_read1(fp, b) >= 0) {
		char *ref, *alt;
		int l_ref, l_alt, i, min_l;
		bcf_get_ref_alt1(b, &l_ref, &ref, &l_alt, &alt);
		min_l = l_ref < l_alt? l_ref : l_alt;
		for (i = 0; i < min_l && ref[i] == alt[i]; ++i);
		s.l = 0;
		kputs(h->id[BCF_DT_CTG][b->rid].key, &s);
		kputc(':', &s); kputw(b->pos + 1 + i, &s);
		kputc(':', &s); kputw(b->rlen - i, &s);
		kputc(':', &s); kputsn(alt + i, l_alt - i, &s);
		puts(s.s);
	}
	bcf_destroy1(b);
	bcf_hdr_destroy(h);

	bgzf_close(fp);
	free(s.s);
	return 0;
}
