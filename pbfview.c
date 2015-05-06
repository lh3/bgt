#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "pbwt.h"

int main(int argc, char *argv[])
{
	int c, in_txt = 0, out_pbf = 0, m_sub = 0, n_sub = 0, *sub = 0, shift = 13;
	int64_t row_start = 0, n_rec = -1;
	pbf_t *out = 0;

	while ((c = getopt(argc, argv, "Sbc:r:n:s:")) >= 0) {
		if (c == 'S') in_txt = 1;
		else if (c == 'b') out_pbf = 1;
		else if (c == 'r') row_start = atol(optarg);
		else if (c == 'n') n_rec = atol(optarg);
		else if (c == 's') shift = atoi(optarg);
		else if (c == 'c') {
			if (n_sub == m_sub) {
				m_sub = m_sub? m_sub<<1 : 4;
				sub = (int*)realloc(sub, m_sub * sizeof(int));
			}
			sub[n_sub++] = atol(optarg);
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: pbfview [options] <in.pbf>|<in.pim>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -S       input is PIM (portable integer matrix format)\n");
		fprintf(stderr, "  -b       output PBF (positional BWT format)\n");
		fprintf(stderr, "  -s INT   write S array every 1<<INT rows (effective with -b) [%d]\n", shift);
		fprintf(stderr, "  -r INT   start decoding from row INT (effective w/o -S) [0]\n");
		fprintf(stderr, "  -n INT   read INT rows starting from -r (effective w/o -S) [inf]\n");
		fprintf(stderr, "  -c INT   decode column INT (there can be multiple -c; effective w/o -S) [inf]\n");
		return 1;
	}

	if (n_rec < 0) n_rec = INT64_MAX;
	if (in_txt) {
		char magic[256];
		FILE *fp;
		int i, j, m, g;
		uint8_t **a;
		fp = strcmp(argv[optind], "-")? fopen(argv[optind], "r") : stdin;
		fscanf(fp, "%s%d%d", magic, &m, &g); // unsafe!!!
		if (out_pbf) out = pbf_open_w(0, m, g, shift);
		else printf("PIM1 %d %d\n", m, g);
		a = (uint8_t**)calloc(g, sizeof(void*));
		for (j = 0; j < g; ++j) a[j] = (uint8_t*)calloc(m, 1);
		while (!feof(fp)) {
			for (i = 0; i < m; ++i) {
				long x;
				if (feof(fp)) break;
				fscanf(fp, "%ld", &x);
				if (out) {
					for (j = 0; j < g; ++j)
						a[j][i] = x>>j & 1;
				} else {
					if (i) putchar(' ');
					printf("%ld", x);
				}
			}
			if (i < m) break;
			if (out) pbf_write(out, a);
			else putchar('\n');
		}
		for (j = 0; j < g; ++j) free(a[j]);
		free(a);
		fclose(fp);
		fflush(stdout);
	} else {
		pbf_t *in;
		int64_t i;
		int m, g, j, k;
		const uint8_t **a;

		in = pbf_open_r(argv[optind]);
		m = n_sub > 0? n_sub : pbf_get_m(in);
		g = pbf_get_g(in);
		if (out_pbf) out = pbf_open_w(0, m, g, shift);
		else printf("PIM1 %d %d\n", m, g);
		if (row_start > 0) pbf_seek(in, row_start);
		if (n_sub > 0) pbf_subset(in, n_sub, sub);
		for (i = 0; i < n_rec && (a = pbf_read(in)) != 0; ++i) {
			if (!out) {
				for (j = 0; j < m; ++j) {
					uint64_t x = 0;
					if (j) putchar(' ');
					for (k = 0; k < g; ++k)
						x |= (uint64_t)a[k][j]<<k;
					printf("%llu", (unsigned long long)x);
				}
				putchar('\n');
			} else pbf_write(out, (uint8_t**)a);
		}
		pbf_close(in);
	}
	if (out) pbf_close(out);
	free(sub);
	return 0;
}
