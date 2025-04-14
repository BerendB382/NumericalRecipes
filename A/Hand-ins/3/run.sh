#!/bin/bash

echo "Downloading satgals..."
if [ ! -e satgals_m11.txt ]; then
  wget home.strw.leidenuniv.nl/~daalen/Handin_files/Vandermonde.txt
fi

if [ ! -e satgals_m12.txt ]; then
  wget home.strw.leidenuniv.nl/~daalen/Handin_files/Vandermonde.txt
fi

if [ ! -e satgals_m13.txt ]; then
  wget home.strw.leidenuniv.nl/~daalen/Handin_files/Vandermonde.txt
fi

if [ ! -e satgals_m14.txt ]; then
  wget home.strw.leidenuniv.nl/~daalen/Handin_files/Vandermonde.txt
fi

if [ ! -e satgals_m15.txt ]; then
  wget home.strw.leidenuniv.nl/~daalen/Handin_files/Vandermonde.txt
fi

# Script that returns a plot
echo "Run the Satellites script ..."
python3 Satellites2.py > satellite_output.txt

echo "Generating the pdf"

pdflatex Nieuwhof.tex
bibtex Nieuwhof.aux




