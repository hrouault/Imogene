/*
 * =====================================================================================
 *
 *       Filename:  extract.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06.08.2011 13:03:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <sys/stat.h>

using namespace std;

#include "extract_cmdline.h"
#include "const.hpp"
#include "sequence.hpp"
#include "tree.hpp"

extract_args_info extract_args;

   void
extractfromcoord(const char * coordfile)
{

   // INPUT COORDINATES
   ifstream coordinates(coordfile);

   vcoord coords;
   back_insert_iterator<vcoord> dest(coords);
   copy(iiscoord(coordinates),iiscoord(),dest);
   
   // ALIGNMENTS COORDINATES
   ifstream align;

   if (species=="droso"){
      cout << "Reading droso alignments..." << endl;
      align.open(DATA_PATH"/droso/align.dat");
   } else if (species=="eutherian"){
      cout << "Reading eutherian alignments..." << endl;
      align.open(DATA_PATH"/eutherian/align.dat");
   }

   alignscoord=loadcoordconserv(align);

   align.close();

   // SEQUENCE EXTRACTION
   vseq seqs;
   cout << "Extraction start" << endl;
   for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
//      cout << chromfromint(ivc->chrom) << endl;
//      cout << ivc->start << endl;
//      cout << ivc->stop << endl;
      Sequence seqtoimport=coordtoseq(*ivc);
      cout << seqtoimport.name << endl;
//      cout << seqtoimport.seqs[0] << endl;
      if (seqtoimport.species[0] && seqtoimport.nbtb>0)  
         seqs.push_back(seqtoimport);
   }
   
   cout << "Extraction stop" << endl;

   cout << "Writing sequences in align/..." << endl;
   
   stringstream basename;
   if (extract_args.background_given){
      if (species=="droso"){
         mkdir(DATA_PATH"/droso/background",S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);      
         basename << DATA_PATH"/droso/background/";
      } else if (species=="eutherian"){
         mkdir(DATA_PATH"/eutherian/background",S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);      
         basename << DATA_PATH"/eutherian/background/";
      }
   } else {
      mkdir("align",S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH); 
      basename << "align/";
   }

   ofstream outf;
   for (ivseq iv=seqs.begin();iv!=seqs.end();iv++){
      Sequence seq=*iv;
      stringstream file;
      file << basename;
      file << seq.name << "_";
      file << chromfromint(seq.chrom) << "_";
      file << seq.start << "_";
      file << seq.stop << ".fa";
      outf.open(file.str().c_str());
      Sequence & s=seq;
      for (int i=0;i<nbspecies;i++){
         if (s.species[i]){
            if (i==0) outf << ">" << numtospecies(i) << " " <<
               "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
            else  outf << ">" << numtospecies(i) << endl;
            outf << s.seqsrealigned[i] << endl;
         }; 
      }
      outf.close();
   }
}

// this one is using specified folder and uses short names for fasta files

//vseq seqs;
//
//   void
//extracttofasta(string folder)
//{
//   ofstream outf;
//   string pname("");
//   int pnum(1);
//   for (ivseq iv=seqs.begin();iv!=seqs.end();iv++){
//      Sequence seq=*iv;
//      if (seq.species[0] && seq.nbtb>0){ 
//         stringstream file;
//         file << folder;
//         file << seq.name;
//         if (seq.name==pname){ 
//            file << "_";
//            file << pnum;
//            pnum++;
//         }
//         else {
//            pnum=1;
//            pname=seq.name;
//         }
//         file << ".fa";
//         outf.open(file.str().c_str());
//         Sequence & s=seq;
//         for (int i=0;i<nbspecies;i++){
//            if (s.species[i]){
//               if (i==0) outf << ">" << numtospecies(i) << " " <<
//                  "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
//               else  outf << ">" << numtospecies(i) << endl;
//               outf << s.seqsrealigned[i] << endl;
//            }; 
//         }
//         outf.close();
//      }
//   }
//}

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
      exit(1);
   if (!strcmp(extract_args.species_arg,"droso")){
      species="droso";
      nbspecies=12;
   } else if (!strcmp(extract_args.species_arg,"eutherian")){
      species="eutherian";
      nbspecies=12;
   }

   extractfromcoord(extract_args.input_arg);

   cout << "exit normally" << endl;
   return 1;

}		/* -----  end of function extract  ----- */
