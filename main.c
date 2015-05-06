#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

int main_ucf2bgt(int argc, char *argv[]);
int main_view(int argc, char *argv[]);

static int usage()
{
	fprintf(stderr, "Usage: bgt <command> <argument>\n");
	fprintf(stderr, "Commands:\n");
	fprintf(stderr, "  import       convert unary VCF to BGT\n");
	fprintf(stderr, "  view         extract from BGT\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "ucf2bgt") == 0) return main_ucf2bgt(argc-1, argv+1);
	if (strcmp(argv[1], "import") == 0) return main_ucf2bgt(argc-1, argv+1);
	if (strcmp(argv[1], "view") == 0) return main_view(argc-1, argv+1);
	else {
		fprintf(stderr, "[E::%s] unrecognized command '%s'\n", __func__, argv[1]);
		return 1;
	}
}
