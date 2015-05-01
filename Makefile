CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat -Wno-unused-function
CPPFLAGS=
OBJS=		bgzf.o hts.o vcf.o bgt.o ucf2bgt.o
INCLUDES=
LIBS=		-lpthread -lz
PROG=		bgt

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bgt:$(OBJS) main.o
		$(CC) $^ -o $@ $(LIBS)

bgzf.o:bgzf.c bgzf.h khash.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DBGZF_MT -DBGZF_CACHE $(INCLUDES) $< -o $@

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ *.a *.so *.dylib

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bgzf.o: bgzf.h
hts.o: bgzf.h hts.h kseq.h khash.h ksort.h
ucf2bgt.o: vcf.h bgzf.h hts.h kstring.h
vcf.o: kstring.h bgzf.h vcf.h hts.h khash.h kseq.h
