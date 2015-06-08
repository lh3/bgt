#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include "bgt.h"
#include "kexpr.h"
#include "fmf.h"

void *bed_read(const char *fn);
void bed_destroy(void *_h);
char **hts_readlines(const char *fn, int *_n);

int main_view(int argc, char *argv[])
{
	int i, c, n_files = 0, out_bcf = 0, clevel = -1, multi_flag = 0, excl = 0, not_vcf = 0, u_set = 0;
	long seekn = -1, n_rec = LONG_MAX, n_read = 0;
	bgtm_t *bm = 0;
	bcf1_t *b;
	htsFile *out = 0;
	char modew[8], *reg = 0, *site_flt = 0;
	void *bed = 0;
	int n_groups = 0;
	char *gexpr[BGT_MAX_GROUPS], **al = 0, *fmt = 0;
	int n_al = 0;
	bgt_file_t **files = 0;

	while ((c = getopt(argc, argv, "ubs:r:l:AGB:ef:g:a:i:n:SHt:")) >= 0) {
		if (c == 'b') out_bcf = 1;
		else if (c == 'r') reg = optarg;
		else if (c == 'l') clevel = atoi(optarg);
		else if (c == 'e') excl = 1;
		else if (c == 'u') u_set = 1;
		else if (c == 'B') bed = bed_read(optarg);
		else if (c == 'A') multi_flag |= BGT_F_SET_AC;
		else if (c == 'G') multi_flag |= BGT_F_NO_GT;
		else if (c == 'S') multi_flag |= BGT_F_NO_GT | BGT_F_CNT_AL, not_vcf = 1;
		else if (c == 'H') multi_flag |= BGT_F_NO_GT | BGT_F_CNT_HAP, not_vcf = 1;
		else if (c == 'i') seekn = atol(optarg) - 1;
		else if (c == 'n') n_rec = atol(optarg);
		else if (c == 'f') site_flt = optarg;
		else if (c == 't') fmt = optarg, not_vcf = 1;
		else if (c == 's' && n_groups < BGT_MAX_GROUPS) gexpr[n_groups++] = optarg;
		else if (c == 'a') {
			al = hts_readlines(optarg, &n_al);
			if (n_al > BGT_MAX_ALLELES) {
				for (i = BGT_MAX_ALLELES; i < n_al; ++i) free(al[i]);
				n_al = BGT_MAX_ALLELES;
				fprintf(stderr, "[W::%s] only the first %d alleles are parsed\n", __func__, BGT_MAX_ALLELES);
			}
		}
	}
	if (seekn < 0) seekn = 0;
	if (clevel > 9) clevel = 9;
	if (u_set) clevel = 0, out_bcf = 1;
	if (n_groups > 1) multi_flag |= BGT_F_SET_AC;
	if (argc - optind < 1) {
		fprintf(stderr, "Usage: bgt %s [options] <bgt-prefix> [...]", argv[0]);
		fputc('\n', stderr);
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Sample selection:\n");
		fprintf(stderr, "    -s EXPR      samples list (,sample1,sample2 or a file or expr; see Notes below) [all]\n");
		fprintf(stderr, "  Site selection:\n");
		fprintf(stderr, "    -r STR       region [all]\n");
		fprintf(stderr, "    -B FILE      extract variants overlapping BED FILE []\n");
		fprintf(stderr, "    -e           exclude variants overlapping BED FILE (effective with -B)\n");
		fprintf(stderr, "    -i INT       process from the INT-th record (1-based) []\n");
		fprintf(stderr, "    -n INT       process at most INT records []\n");
		fprintf(stderr, "    -a STR       alleles list chr:1basedPos:refLen:seq (,allele1,allele2 or a file) []\n");
		fprintf(stderr, "    -f STR       frequency filters []\n");
		fprintf(stderr, "  VCF output:\n");
		fprintf(stderr, "    -b           BCF output (effective without -S/-H)\n");
		fprintf(stderr, "    -l INT       compression level for BCF [default]\n");
		fprintf(stderr, "    -u           equivalent to -bl0 (overriding -b and -l)\n");
		fprintf(stderr, "    -G           don't output sample genotypes\n");
		fprintf(stderr, "    -A           write AC/AN to the INFO field (auto applied with -f or multipl -s)\n");
		fprintf(stderr, "  Non-VCF output:\n");
		fprintf(stderr, "    -S           show samples with a set of alleles (with -a)\n");
		fprintf(stderr, "    -H           count of haplotypes with a set of alleles (with -a)\n");
		fprintf(stderr, "    -t STR       comma-delimited list of fields to output. Accepted variables:\n");
		fprintf(stderr, "                 AC, AN, AC#, AN#, CHROM, POS, END, REF, ALT (# for a group number)\n");
		fprintf(stderr, "Notes:\n");
		fprintf(stderr, "  For option -s, EXPR can be one of:\n");
		fprintf(stderr, "    1) comma-delimited sample list following a colon/comma. e.g. -s,NA12878,NA12044\n");
		fprintf(stderr, "    2) space-delimited file with the first column giving a sample name. e.g. -s samples.txt\n");
		fprintf(stderr, "    3) expression if .spl file contains metadata. e.g.: -s\"gender=='M'&&population!='CEU'\"\n");
		fprintf(stderr, "  If multiple -s is specified, the AC/AN of the first group will be written to VCF INFO AC1/AN1,\n");
		fprintf(stderr, "  the second to AC2/AN2, etc.\n");
		return 1;
	}

	if ((multi_flag&(BGT_F_CNT_AL|BGT_F_CNT_HAP)) && n_al == 0) {
		fprintf(stderr, "[E::%s] at least one -a must be specified when -S/-H is in use.\n", __func__);
		return 1;
	}

	n_files = argc - optind;
	files = (bgt_file_t**)calloc(n_files, sizeof(bgt_file_t*));
	for (i = 0; i < n_files; ++i) files[i] = bgt_open(argv[optind+i]);

	bm = bgtm_reader_init(n_files, files);
	bgtm_set_flag(bm, multi_flag);
	if (site_flt) bgtm_set_flt_site(bm, site_flt);
	if (reg) bgtm_set_region(bm, reg);
	if (bed) bgtm_set_bed(bm, bed, excl);
	if (fmt) bgtm_set_table(bm, fmt);
	if (seekn >= 0) bgtm_set_start(bm, seekn);
	for (i = 0; i < n_al; ++i)
		bgtm_add_allele(bm, al[i]);
	if (multi_flag&(BGT_F_CNT_AL|BGT_F_CNT_HAP)) {
		for (i = 1; i < bm->n_al; ++i)
			if (strcmp(bm->al[0].chr.s, bm->al[i].chr.s) != 0)
				fprintf(stderr, "[W::%s] alleles on different chr; unexpected errors may happen!\n", __func__);
	}
	for (i = 0; i < n_groups; ++i)
		bgtm_add_group(bm, gexpr[i]);
	bgtm_prepare(bm); // bgtm_prepare() generates the VCF header

	if (!not_vcf) {
		strcpy(modew, "w");
		if (out_bcf) strcat(modew, "b");
		sprintf(modew + strlen(modew), "%d", clevel);
		out = hts_open("-", modew, 0);
		vcf_hdr_write(out, bm->h_out);
	}

	b = bcf_init1();
	while (bgtm_read(bm, b) >= 0 && n_read < n_rec) {
		if (out) vcf_write1(out, bm->h_out, b);
		if (fmt && bm->n_fields > 0) puts(bm->tbl_line.s);
		++n_read;
	}
	bcf_destroy1(b);

	if (not_vcf && bm->n_al > 0) {
		if (bm->flag & BGT_F_CNT_HAP) {
			bgt_hapcnt_t *hc;
			int n_hap;
			char *s;
			hc = bgtm_hapcnt(bm, &n_hap);
			s = bgtm_hapcnt_print_destroy(bm, n_hap, hc);
			fputs(s, stdout);
			free(s);
		}
		if (bm->flag & BGT_F_CNT_AL) {
			char *s;
			s = bgtm_alcnt_print(bm);
			fputs(s, stdout);
			free(s);
		}
	}

	if (out) hts_close(out);
	bgtm_reader_destroy(bm);
	if (bed) bed_destroy(bed);
	for (i = 0; i < n_files; ++i) bgt_close(files[i]);
	free(files);
	for (i = 0; i < n_al; ++i) free(al[i]);
	free(al);
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
