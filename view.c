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
	for (i = 1; i <= n_groups; ++i) {
		key[0] = 'A', key[1] = 'N', key[2] = '0' + i; ke_set_int(flt->ke, key, gan[i]);
		key[0] = 'A', key[1] = 'C', key[2] = '0' + i; ke_set_int(flt->ke, key, gac1[i]);
	}
	is_true = !!ke_eval_int(flt->ke, &err);
	if (err) return 0;
	return !is_true;
}

static char **get_samples(const char *expr, int *n, int n_pre, char *const* prefix)
{
	int err, i, j, is_file = 0;
	khint_t k;
	char **s;
	kexpr_t *ke;
	khash_t(s2i) *h;
	FILE *fp;

	*n = 0;
	if ((fp = fopen(expr, "r")) != 0) { // test if expr is a file
		is_file = 1;
		fclose(fp);
	}
	if (*expr == ':' || (*expr != '?' && is_file))
		return hts_readlines(expr, n);
	ke = ke_parse(*expr == '?'? expr+1 : expr, &err);
	if (err) return 0;
	h = kh_init(s2i);
	for (j = 0; j < n_pre; ++j) {
		fmf_t *f;
		char *tmpfn;
		int absent;
		tmpfn = (char*)calloc(strlen(prefix[j]) + 5, 1);
		strcat(strcpy(tmpfn, prefix[j]), ".spl");
		f = fmf_read(tmpfn);
		free(tmpfn);
		if (f == 0) continue;
		for (i = 0; i < f->n_rows; ++i)
			if (fmf_test(f, i, ke)) {
				k = kh_put(s2i, h, f->rows[i].name, &absent);
				if (absent) kh_key(h, k) = strdup(f->rows[i].name);
			}
		fmf_destroy(f);
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
	char modew[8], *reg = 0;
	flt_aux_t flt;
	void *bed = 0;
	int n_groups = 0;
	char *gexpr[BGT_MAX_GROUPS];
	int n_a = 0, m_a = 0;
	bgt_allele_t *a = 0;

	memset(&flt, 0, sizeof(flt_aux_t));
	assert(strcmp(argv[0], "view") == 0 || strcmp(argv[0], "sview") == 0 || strcmp(argv[0], "mview") == 0);
	is_multi = strcmp(argv[0], "sview")? 1 : 0;
	while ((c = getopt(argc, argv, "bs:r:l:AGB:ef:g:a:")) >= 0) {
		if (c == 'b') out_bcf = 1;
		else if (c == 'r') reg = optarg;
		else if (c == 'l') clevel = atoi(optarg);
		else if (c == 'e') excl = 1;
		else if (c == 'B') bed = bed_read(optarg);
		else if (is_multi && c == 'A') multi_flag |= BGT_F_SET_AC;
		else if (is_multi && c == 'G') multi_flag |= BGT_F_NO_GT;
		else if (is_multi && c == 'f') {
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
	if (clevel > 9) clevel = 9;
	if (n_groups > 1) multi_flag |= BGT_F_SET_AC;
	if (reg && n_a > 0) {
		fprintf(stderr, "[W::%s] -r and -a can't be specified at the same time. -r is igored.\n", __func__);
		reg = 0;
	}
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
			fprintf(stderr, "Notes:\n");
			fprintf(stderr, "  For option -s, EXPR can be one of:\n");
			fprintf(stderr, "    1) comma-delimited sample list following a colon. e.g. -s:NA12878,NA12044\n");
			fprintf(stderr, "    2) space-delimited file with the first column giving a sample name. e.g. -s samples.txt\n");
			fprintf(stderr, "    3) expression if .spl file contains metadata. e.g.: -s\"gender=='M'&&population!='CEU'\"\n");
			fprintf(stderr, "  If multiple -s is specified, the AC/AN of the first group will be written to VCF INFO AC1/AN1,\n");
			fprintf(stderr, "  the second to AC2/AN2, etc.\n");
		}
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

	if (!is_multi) {
		bgt = bgt_open(argv[optind]);
		if (reg) bgt_set_region(bgt, reg);
		if (bed) bgt_set_bed(bgt, bed, excl);
	} else {
		bm = bgtm_open(argc - optind, &argv[optind]);
		bgtm_set_flag(bm, multi_flag);
		if (flt.ke) bgtm_set_filter(bm, filter_func, &flt);
		if (reg) bgtm_set_region(bm, reg);
		if (bed) bgtm_set_bed(bm, bed, excl);

		if (n_groups > 0) {
			char **g[BGT_MAX_GROUPS], **samples;
			int j, k, ng[BGT_MAX_GROUPS], n_samples = 0;
			for (i = 0; i < n_groups; ++i) {
				g[i] = get_samples(gexpr[i], &ng[i], argc - optind, argv + optind);
				n_samples += ng[i];
			}
			samples = (char**)malloc(n_samples * sizeof(char*));
			for (i = k = 0; i < n_groups; ++i)
				for (j = 0; j < ng[i]; ++j)
					samples[k++] = g[i][j];
			bgtm_set_samples(bm, n_samples, samples);
			free(samples);
			if (n_groups > 1)
				for (i = 0; i < n_groups; ++i)
					bgtm_add_group(bm, ng[i], g[i]);
			for (i = 0; i < n_groups; ++i) {
				for (j = 0; j < ng[i]; ++j) free(g[i][j]);
				free(g[i]);
			}
		}
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
	if (n_a) free(reg);
	for (i = 0; i < n_a; ++i) free(a[i].chr.s);
	free(a);
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
