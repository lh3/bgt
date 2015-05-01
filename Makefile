CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat -Wno-unused-function
DFLAGS=
LOBJS=		bgzf.o hts.o vcf.o
INCLUDES=
LIBPATH=

.SUFFIXES:.c .o
.PHONY:lib

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

libhts.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

bgzf.o:bgzf.c bgzf.h khash.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -DBGZF_MT -DBGZF_CACHE $(INCLUDES) $< -o $@

hts.o:hts.c bgzf.h bgzf.h khash.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

vcf.o:vcf.h bgzf.h kstring.h khash.h hts.h

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ *.a *.so *.dylib
