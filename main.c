#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#define BGT_VERSION "r130"

int main_ucf2bgt(int argc, char *argv[]);
int main_view(int argc, char *argv[]);

static int usage()
{
	fprintf(stderr, "Usage: bgt <command> <argument>\n");
	fprintf(stderr, "Commands:\n");
	fprintf(stderr, "  import       convert unary VCF to BGT\n");
//	fprintf(stderr, "  sview        single-BGT view (obsolete)\n"); // now this is for debugging only. use mview instead
	fprintf(stderr, "  view         extract from BGT\n");
	fprintf(stderr, "  version      show version number\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "ucf2bgt") == 0) return main_ucf2bgt(argc-1, argv+1);
	else if (strcmp(argv[1], "import") == 0) return main_ucf2bgt(argc-1, argv+1);
	else if (strcmp(argv[1], "sview") == 0) return main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "view") == 0) return main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "mview") == 0) return main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		puts(BGT_VERSION);
		return 0;
	} else {
		fprintf(stderr, "[E::%s] unrecognized command '%s'\n", __func__, argv[1]);
		return 1;
	}
}
