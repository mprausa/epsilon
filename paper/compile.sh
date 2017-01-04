#!/bin/bash

pdflatex epsilon.tex || exit $?
bibtex epsilon.aux || exit $?
pdflatex epsilon.tex || exit $?
pdflatex epsilon.tex || exit $?
exit 0


