CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat -Wno-unused-function
CPPFLAGS=
OBJS=		bgzf.o hts.o vcf.o pbwt.o bgt.o import.o view.o
INCLUDES=
LIBS=		-lpthread -lz
PROG=		bgt pbfview

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bgt:$(OBJS) main.o
		$(CC) $^ -o $@ $(LIBS)

pbfview:pbfview.o pbwt.o
		$(CC) $^ -o $@

bgzf.o:bgzf.c bgzf.h khash.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DBGZF_MT -DBGZF_CACHE $(INCLUDES) $< -o $@

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ *.a *.so *.dylib $(PROG) pbwt.aux pbwt.pdf pbwt.log

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bgt.o: bgt.h vcf.h bgzf.h hts.h kstring.h pbwt.h khash.h ksort.h
bgzf.o: bgzf.h
hts.o: bgzf.h hts.h kseq.h khash.h ksort.h
import.o: vcf.h bgzf.h hts.h kstring.h pbwt.h
pbfview.o: pbwt.h
pbwt.o: pbwt.h ksort.h
vcf.o: kstring.h bgzf.h vcf.h hts.h khash.h kseq.h
view.o: bgt.h vcf.h bgzf.h hts.h kstring.h pbwt.h
