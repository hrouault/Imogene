#!/bin/sh
#    
# Copyright (C) 2006-2011 Hervé Rouault <rouault@lps.ens.fr>
# Copyright (C) 2009-2011 Marc Santolini <santolin@lps.ens.fr>
# 
# This file is part of Imogene.
# 
# Imogene is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Imogene is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Imogene; see the file COPYING  If not, see <http://www.gnu.org/licenses/>.

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



