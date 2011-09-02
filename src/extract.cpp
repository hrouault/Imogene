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
#include <string>

using namespace std;

#include "extract_cmdline.h"
#include "const.hpp"
#include "sequence.hpp"

extract_args_info extract_args;

   void
extractfromcoord(const char * coordfile)
{
   cout << "extraction start" << endl;
   ifstream coordinates(coordfile);

   vcoord coords;
   back_insert_iterator<vcoord> dest(coords);
   copy(iiscoord(coordinates),iiscoord(),dest);

   if (species=="drosos"){
      cout << "extract droso alignments" << endl;
      for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
         cout << ivc->chrom << endl;
         cout << ivc->start << endl;
         cout << ivc->stop << endl;
         Sequence seqtoimport=coordtoseq(*ivc);
         cout << seqtoimport.name << endl;
         cout << seqtoimport.seqs[0] << endl;
      }

   } else if (species=="eutherian"){

   }

//   cout << coords << endl;

   // *** The following is incorrect
//   vcoord aligns=loadcoordconserv(alignmelano);
//   cout << "Align set size : " << melalignscoord.size() << endl;
//   alignmelano.close();


//   interestmelano.close();
   cout << "inputcoords finished" << endl;

//
//
//   ofstream outf;
//   for (ivseq iv=regints.begin();iv!=regints.end();iv++){
//      Sequence seq=*iv;
//      if (seq.species[0] && seq.nbtb>0){ 
//         stringstream file;
//         file << folder;
//         file << seq.name << "_";
//         file << chromfromint(seq.chrom) << "_";
//         file << seq.start << "_";
//         file << seq.stop << ".fa";
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
}

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
   if (strcmp(extract_args.species_arg,"drosos")){
      species="drosos";
   } else if (strcmp(extract_args.species_arg,"drosos")){
      species="eutherian";
   }

   extractfromcoord(extract_args.input_arg);

   extract_cmdline_parser_free(&extract_args);

   return 1;
}		/* -----  end of function extract  ----- */
