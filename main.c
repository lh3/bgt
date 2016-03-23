#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#define BGT_VERSION "1.0-r282-dirty"

int main_import(int argc, char *argv[]);
int main_view(int argc, char *argv[]);
int main_getalt(int argc, char *argv[]);
int main_bcfidx(int argc, char *argv[]);
int main_fmf(int argc, char *argv[]);
int main_atomize(int argc, char *argv[]);

static int usage()
{
	fprintf(stderr, "Usage: bgt <command> <argument>\n");
	fprintf(stderr, "Commands:\n");
	fprintf(stderr, "  import       convert VCF to BGT\n");
	fprintf(stderr, "  atomize      atomize VCF\n");
	fprintf(stderr, "  view         extract from BGT\n");
	fprintf(stderr, "  fmf          manipulate FMF files\n");
	fprintf(stderr, "  bcfidx       (re)index BCF with record number index\n");
	fprintf(stderr, "  version      show version number\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "import") == 0) return main_import(argc-1, argv+1);
	else if (strcmp(argv[1], "atomize") == 0) return main_atomize(argc-1, argv+1);
	else if (strcmp(argv[1], "view") == 0 || strcmp(argv[1], "mview") == 0 ) return main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "fmf") == 0 ) return main_fmf(argc-1, argv+1);
	else if (strcmp(argv[1], "getalt") == 0) return main_getalt(argc-1, argv+1);
	else if (strcmp(argv[1], "bcfidx") == 0) return main_bcfidx(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		puts(BGT_VERSION);
		return 0;
	} else {
		fprintf(stderr, "[E::%s] unrecognized command '%s'\n", __func__, argv[1]);
		return 1;
	}
}
