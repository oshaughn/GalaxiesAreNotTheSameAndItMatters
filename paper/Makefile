
SRC=paper.tex


DVIS    = $(SRC:.tex=.dvi) $(ART:.tex=.dvi)
PDFS	= $(SRC:.tex=.pdf) $(ART:.tex=.dvi)


pdf:	 $(PDFS)
	echo "Done processing" `pwd`;


dvi:	 $(DVIS)
	echo "Done processing" `pwd`;


%.pdf:	%.dvi links
	pdflatex -interaction=batchmode $(<:.dvi=.tex)

%.dvi:  %.tex 
	rm -f *.bbl
	latex -src-specials -interaction=batchmode ${TEXOPTIONS} $<; exit 0   # May fail
	(bibtex $(<:.tex=); ./bib-nocomments.sh $(<:.tex=.bbl);    exit 0)
	echo `for i in bu[0-9]*.aux; do bibtex $$i;   done`
#	bibexport.sh $(<:.tex=.aux) -o $(<:.tex=export.bib)
	latex -src-specials -interaction=batchmode ${TEXOPTIONS}  $<; exit 0   # may fail
	latex -src-specials -interaction=batchmode ${TEXOPTIONS} $<; exit 0   # force run to completion -- for some reason this tends not to do well
