
# Mainfile:
m=main

all:
	latex -interaction=nonstopmode -halt-on-error $m.tex
draft:
	latex -interaction=nonstopmode -draftmode -halt-on-error $m.tex
full:
	make
	dvips -o $m.ps $m.dvi
	ps2pdf $m.ps
clean:
	for f in $(find . -name '*.aux' -o -name '*.dvi' -o -name '*.lof' -o -name '*.log' -o -name '*.lot' -o -name '*.out' -o -name '*.ps' -o -name '*.toc' -o -name '*.bbl' -o -name '*.blg' -o -name '*.lox' -o -name '*.nav' -o -name '*.snm' -o -name '*.pdf'); do rm $f; done
bibtex:
	bibtex main1.aux
release:
	make draft
	make bibtex
	make draft
	make draft
	make full


