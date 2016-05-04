#
# Makefile for Thesis
#
# Using AASTex and BibTeX
#

AGR_FILES:=$(shell ls *.agr */*.agr)
EPS_FIGURES2=$(shell ls *.eps */*.eps)
EPS_FIGURES=$(AGR_FILES:%.agr=%.eps)
PDF_FIGURES=$(AGR_FILES:%.agr=%.pdf)

all : thesis.pdf

thesis.pdf : $(PDF_FIGURES) thesis.tex thesis.bib */*.tex 
	pdflatex -draftmode -file-line-error thesis
	bibtex thesis
	pdflatex -draftmode -file-line-error thesis
	pdflatex -file-line-error thesis

thesis.dvi : $(EPS_FIGURES) $(EPS_FIGURES2) thesis.tex thesis.bib */*.tex 
	latex thesis
	bibtex thesis
	latex thesis
	latex thesis

thesis.ps : thesis.dvi
	dvips thesis

clean:
	rm -fr thesis.dvi thesis.ps thesis.pdf thesis.aux thesis.log \
	thesis.lof thesis.lot thesis.toc thesis.blg thesis.bbl */*.aux *~ 

%.eps:%.agr
	DISPLAY=:0.0 xmgrace -hardcopy -hdevice EPS -printfile $@ $<

%.pdf:%.agr
	DISPLAY=:0.0 xmgrace -hardcopy -hdevice EPS -printfile epsfig.eps $<
	epstopdf --outfile=$@ epsfig.eps

%.eps:%.fig
	fig2dev -L eps $< $@

