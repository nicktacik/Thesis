SHELL := /bin/bash

# Makefile for papers
PAPER := paper

MANUAL_FIGURES :=
AGR_FILES   := $(shell find . -name '*.agr')
AGR_FIGURES = $(AGR_FILES:%.agr=%.pdf)

AUTO_FIGURES   := $(AGR_FIGURES)

.DEFAULT_GOAL := $(PAPER).pdf

# Get References SVN repo if it does not yet exist
#References/References.bib:
#	@echo "Checking out *SVN* References repository"
#	svn checkout svn+ssh://svn@sxs-archive.tapir.caltech.edu/home/SVNroot/References/trunk References

nick: paper.pdf
	open /Applications/Adobe*Reader*.app paper.pdf

$(PAPER).pdf: $(PAPER).tex $(AUTO_FIGURES) $(MANUAL_FIGURES) ../References/References.bib
	pdflatex -draftmode -file-line-error $<
	bibtex $(PAPER)
	pdflatex -draftmode -file-line-error $<
	pdflatex -file-line-error $<

# Automatically make pdf files from .gpi (gnuplot) scripts
%.pdf: %.gpi
	gnuplot $<

# Automatically make pdf files from .fig files
%.pdf: %.fig
	fig2mpdf $<

# Automatically make pdf files from .svg files
%.pdf: %.svg
	inkscape -D -z --file=$< --export-pdf=$@

# Automatically make pdf files from .agr files
%.pdf:%.agr
	#DISPLAY=:0.0 xmgrace -hardcopy -hdevice EPS -printfile epsfig.eps $<
	xmgrace -hardcopy -hdevice EPS -printfile epsfig.eps $<
	epstopdf --outfile=$@ epsfig.eps
	rm -f epsfig.eps

# Automatically make pdf files from .py files
%.pdf:%.py
	python $<

clean:
	@rm -f *.aux
	@rm -f *.out
	@rm -f *.log
	@rm -f *.bbl
	@rm -f *.blg
	@rm -f *.bib
	@rm -f $(AUTO_FIGURES)
	@rm -f $(PAPER).pdf
