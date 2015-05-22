CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat -Wno-unused-function
CPPFLAGS=
OBJS=		kexpr.o bgzf.o hts.o fmf.o vcf.o atomic.o bedidx.o pbwt.o bgt.o import.o view.o
INCLUDES=
LIBS=		-lpthread -lz
PROG=		bgt
PROG_EXTRA= pbfview kexpr fmf

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

extra:$(PROG_EXTRA)

bgt:$(OBJS) main.o
		$(CC) $^ -o $@ $(LIBS)

pbfview:pbfview.o pbwt.o
		$(CC) $^ -o $@

kexpr:kexpr.c kexpr.h
		$(CC) $(CFLAGS) -DKE_MAIN $< -o $@ -lm

fmf:fmf.c kexpr.o fmf.h
		$(CC) $(CFLAGS) -DFMF_MAIN fmf.c kexpr.o -o $@ -lz -lm

fmf.o:fmf.c fmf.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DFMF_HAVE_HTS $< -o $@

bgzf.o:bgzf.c bgzf.h khash.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DBGZF_MT -DBGZF_CACHE $(INCLUDES) $< -o $@

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ *.a *.so *.dylib $(PROG) $(PROG_EXTRA) pbwt.aux pbwt.pdf pbwt.log

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

atomic.o: atomic.h vcf.h bgzf.h hts.h kstring.h ksort.h
bedidx.o: ksort.h kseq.h khash.h
bgt.o: bgt.h vcf.h bgzf.h hts.h kstring.h pbwt.h khash.h ksort.h
bgzf.o: bgzf.h
fmf.o: fmf.h kexpr.h kseq.h khash.h kstring.h
hts.o: bgzf.h hts.h kseq.h khash.h ksort.h
import.o: atomic.h vcf.h bgzf.h hts.h kstring.h pbwt.h
kexpr.o: kexpr.h
pbfview.o: pbwt.h
pbwt.o: pbwt.h ksort.h
vcf.o: kstring.h bgzf.h vcf.h hts.h khash.h kseq.h
view.o: bgt.h vcf.h bgzf.h hts.h kstring.h pbwt.h kexpr.h
