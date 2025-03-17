#!/bin/bash

# Script that returns a plot
echo "Run the Satellites script ..."
python3 Satellites.py > satellite_output.txt

# Script that pipes output to a file
echo "Run the HII script ..."
python3 HII-regions.py > HII-regions_output.txt

echo "Generating the pdf"

pdflatex Nieuwhof.tex
bibtex Nieuwhof.aux



