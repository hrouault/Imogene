imogene extract -s eutherian -i coord_mouse.bed
imogene genmot -s eutherian -w 10 -t 11.0 -x 20 -e 2 -a align
imogene display -t 8.0 -s eutherian -n 10 -m motifs.txt -a align --html-ref
imogene display -t 8.0 -s eutherian -n 10 -m motifs.txt -a align --logos --pdf --png

imogene scangen -s eutherian -t 8.0 -x 20 -m motifs.txt --scanwidth=1000
