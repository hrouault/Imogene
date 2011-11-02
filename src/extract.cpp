/*    
 * Copyright (C) 2006-2011 Hervé Rouault <rouault@lps.ens.fr>
 * Copyright (C) 2009-2011 Marc Santolini <santolin@lps.ens.fr>
 *
 * This file is part of Imogene.
 *
 * Imogene is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Imogene is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Imogene.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <sys/stat.h>

using namespace std;

#include "extract_cmdline.h"
#include "extract.hpp"
#include "const.hpp"
#include "sequence.hpp"
#include "tree.hpp"

extract_args_info extract_args;

   void
seq2fasta(Sequence &seq,string folder)
{
   ofstream outf;
   stringstream file;
   file << folder;
   file << seq.name << "_";
   file << chromfromint(seq.chrom) << "_";
   file << seq.start << "_";
   file << seq.stop << ".fa";
   outf.open(file.str().c_str());
   if (outf.fail()){
      cerr << "Cannot open file for fasta recording: " << strerror(errno) << endl;
      exit(EXIT_FAILURE);
   }
   Sequence & s=seq;
   for (unsigned int i=0;i<nbspecies;i++){
      if (s.species[i]){
         if (i==0) outf << ">" << numtospecies(i) << " " <<
            chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
         else  outf << ">" << numtospecies(i) << endl;
         outf << s.seqsrealigned[i] << endl;
      }; 
   }
   outf.close();
}

   void
extractfromcoord(const char * coordfile)
{

   // INPUT COORDINATES
   ifstream coordinates(coordfile);
   if (coordinates.fail()){
      cerr << "Cannot open coordinate file for reading: " << strerror(errno) << endl;
      exit(EXIT_FAILURE);
   }

   vcoord coords;
   back_insert_iterator<vcoord> dest(coords);
   copy(iiscoord(coordinates),iiscoord(),dest);

   // ALIGNMENTS COORDINATES
   ifstream align;

   if (species=="droso"){
      cout << "Reading droso alignments..." << endl;
      align.open( (extract_datapath+"/droso/align.dat").c_str() );
   } else if (species=="eutherian"){
      cout << "Reading eutherian alignments..." << endl;
      align.open( (extract_datapath+"/eutherian/align.dat").c_str() );
   }
   if (align.fail()){
      cerr << "Alignment file opening failed: " << strerror(errno) << endl;
      exit(EXIT_FAILURE);
   }

   alignscoord=loadcoordconserv(align);

   align.close();

   stringstream basename;
   if (extract_args.background_given){
      if (species=="droso"){
         mkdir( (extract_datapath+"/droso/background").c_str() ,S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);      
         basename << extract_datapath+"/droso/background/";
      } else if (species=="eutherian"){
         mkdir( (extract_datapath+"/eutherian/background").c_str() ,S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);      
         basename << extract_datapath+"/eutherian/background/";
      }
   } else {
      mkdir("align",S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH); 
      basename << "align/";
   }

   // SEQUENCE EXTRACTION
   cout << "Extraction start" << endl;
   cout << "Writing sequences in align/..." << endl;
   for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){

      Sequence seqtoimport=coordtoseq(*ivc);
      if (seqtoimport.species[0] && seqtoimport.nbtb>0)  {
         cout << seqtoimport.name << endl;
         seq2fasta(seqtoimport,basename.str());
      }
      //      cout << chromfromint(ivc->chrom) << endl;
      //      cout << ivc->start << endl;
      //      cout << ivc->stop << endl;
      //      cout << seqtoimport.seqs[0] << endl;
   }

   cout << "Extraction stop" << endl;

}

   void
extract_args_init()
{
   if (!strcmp(extract_args.species_arg,"droso")){
      species="droso";
      nbspecies=12;
   } else if (!strcmp(extract_args.species_arg,"eutherian")){
      species="eutherian";
      nbspecies=12;
   }
}

string extract_datapath;

/** 
 * ===  FUNCTION  ======================================================================
 *         Name:  cmd_extract
 *  Description:  Alignment extraction
 * =====================================================================================
 */
   int
cmd_extract(int argc, char **argv)
{

   if ( extract_cmdline_parser(argc, argv, & extract_args)!=0)
      exit(EXIT_FAILURE);

   extract_args_init();

   const char * imo_extract_datapath = getenv( "IMOGENE_DATA" );
   if (imo_extract_datapath==NULL){
      extract_datapath = DATA_PATH;
   } else {
      extract_datapath=imo_extract_datapath;
   }

   extractfromcoord(extract_args.input_arg);

   cout << "exit normally" << endl;
   return EXIT_SUCCESS;

}		/* -----  end of function extract  ----- */
