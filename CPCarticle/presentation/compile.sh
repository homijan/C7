#! /bin/sh
rm -r build
mkdir build && cd build

cp ../nonlocalMHD.tex document.tex

pdflatex document.tex
pdflatex document.tex

mv document.pdf ../nonlocalMHD.pdf
cd ..
rm -r build
