#!/bin/bash

echo "Clearing/creating the plotting directory"
if [ ! -d "plots" ]; then
  mkdir plots
fi
rm -rf plots/*

echo "Check if the sine movie exist"
if [ -e sinemovie.mp4 ]; then
  echo "Remove mp4 file"
  rm sinemovie.mp4
fi

echo "Download Vandermonde points for exercise 2 ..."
if [ ! -e Vandermonde.txt ]; then
  wget home.strw.leidenuniv.nl/~daalen/Handin_files/Vandermonde.txt
fi

# Script that returns a plot
echo "Run the Poisson script ..."
python3 Poisson.py > poisson_output.txt

# Script that pipes output to a file
echo "Run the Vandermonde script ..."
python3 Vandermonde.py > Vandermonde_output.txt

echo "Generating the pdf"

pdflatex Nieuwhof.tex
bibtex Nieuwhof.aux
pdflatex Nieuwhof.tex
pdflatex Nieuwhof.tex


