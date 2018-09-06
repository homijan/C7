#! /bin/bash
pdflatex C7.tex
bibtex C7
pdflatex C7.tex
pdflatex C7.tex

pdflatex AWBS.tex
bibtex AWBS
pdflatex AWBS.tex
pdflatex AWBS.tex
