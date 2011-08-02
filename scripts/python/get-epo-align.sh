#!/bin/sh

echo  "[ `date` ] Getting EMF files from ENSEMBL...";echo
#wget -rc -nH --cut-dirs=4 ftp://ftp.ensembl.org/pub/release-63/emf/ensembl-compara/epo_12_eutherian
wget -c ftp://ftp.ensembl.org/pub/release-63/emf/ensembl-compara/epo_12_eutherian/Compara.12_eutherian_mammals_EPO.chr10_1.emf.gz
mv *chr10_1.emf.gz epo_12_eutherian

cd epo_12_eutherian

   echo  "[ `date` ] Unzipping .gz files...";echo
   gunzip -v *
   echo  "[ `date` ] Creating index file for emf files...";echo
   ls -1 `pwd`/*.emf > files.dat

cd ..

echo  "[ `date` ] Converting emf files to fasta...";echo
./align

cd epo_12_eutherian/
   
   echo  "[ `date` ] Creating index file for fasta files...";echo
   find $PWD -name "*.fa" > foo
   sed "s/.*chr//;s/\// /;s/-/ /;s/\.fa//" <foo > bar
   paste bar foo| tr " " "\t" > files.dat
   rm foo
   rm bar
   
   echo  "[ `date` ] Removing .emf files...";echo
   rm *.emf

cd ..



