#!/bin/sh

mkdir -p display-droso

# Motifs in pdf
imogene display -s droso -a ../extract/align-droso -m ../genmot/motifs-droso.txt --logos --pdf -n 10
mv Motif*.pdf display-droso

# Motifs in png
imogene display -s droso -a ../extract/align-droso -m ../genmot/motifs-droso.txt --logos --png -n 10
mv Motif*.png display-droso

# .tex file displaying motifs on reference sequence (textual)
imogene display -s droso -a ../extract/align-droso -m ../genmot/motifs-droso.txt -n 10 --tex-ref -t 10
mv result.tex display-droso

# .tex files displaying motifs on reference and orthologous sequences (textual)
imogene display -s droso -a ../extract/align-droso -m ../genmot/motifs-droso.txt -n 10 --tex-align -t 10
mv *.tex display-droso

# .svg files of motifs on sequences (graphical) 
imogene display -s droso -a ../extract/align-droso -m ../genmot/motifs-droso.txt -n 10 --svg -t 10
mv *.svg display-droso

# .html summary of both previous commands. Motifs need to be produced in png and pdf first.
imogene display -s droso -a ../extract/align-droso -m ../genmot/motifs-droso.txt -n 10 --html-ref -t 10
mv *.svg display-droso
mv css display-droso
mv results_genmot.html display-droso

# .html summary of scangen result
imogene display -e ../scangen/scangen-droso/result5.dat
mv results_scangen.html display-droso
