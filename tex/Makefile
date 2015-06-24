.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

.pdf.eps:
		pdftops -eps $< $@

all:bgt.pdf

bgt.pdf:bgt.tex bgt.bib
		pdflatex bgt; bibtex bgt; pdflatex bgt; pdflatex bgt;

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ bgt.pdf
