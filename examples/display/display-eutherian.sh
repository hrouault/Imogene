#!/bin/sh

mkdir -p display-eutherian

# Motifs in pdf
imogene display -s eutherian -a ../extract/align-eutherian -m ../genmot/motifs-eutherian.txt --logos --pdf -n 5
mv Motif*.pdf display-eutherian

# Motifs in png
imogene display -s eutherian -a ../extract/align-eutherian -m ../genmot/motifs-eutherian.txt --logos --png -n 5
mv Motif*.png display-eutherian

# .tex file displaying motifs on reference sequence (textual)
imogene display -s eutherian -a ../extract/align-eutherian -m ../genmot/motifs-eutherian.txt -n 5 --tex-ref -t 10
mv result.tex display-eutherian

# .tex files displaying motifs on reference and orthologous sequences (textual)
imogene display -s eutherian -a ../extract/align-eutherian -m ../genmot/motifs-eutherian.txt -n 5 --tex-align -t 10
mv *.tex display-eutherian

# .svg files of motifs on sequences (graphical) 
imogene display -s eutherian -a ../extract/align-eutherian -m ../genmot/motifs-eutherian.txt -n 5 --svg -t 10
mv *.svg display-eutherian

# .html summary of both previous commands. Motifs need to be produced in png and pdf first.
imogene display -s eutherian -a ../extract/align-eutherian -m ../genmot/motifs-eutherian.txt -n 5 --html-ref -t 10
mv *.svg display-eutherian
mv css display-eutherian
mv results_genmot.html display-eutherian

# .html summary of scangen result
imogene display -e ../scangen/scangen-eutherian/result5.dat
mv results_scangen.html display-eutherian
