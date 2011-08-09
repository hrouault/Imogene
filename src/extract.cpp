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

using namespace std;

#include "extract_cmdline.h"

extract_args_info extract_args;

   void
extractfromcoord(string coordfile)
{
   ifstream interestmelano;
   interestmelano.open(args_info.coord_file_arg);
   ifstream alignmelano;
   alignmelano.open("/home/rouault/these/sequence/genomes/seqs.dat"); // list all the alignment files
   melalignscoord=loadcoordconserv(alignmelano);
   cout << "Align set size : " << melalignscoord.size() << endl;
   alignmelano.close();

   vcoord coords;
   back_insert_iterator<vcoord> dest(coords);
   copy(iiscoord(interestmelano),iiscoord(),dest);
   for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
      Sequence seqtoimport=coordtoseq(*ivc);
      if (seqtoimport.drosos[0] && seqtoimport.nbtb>0){
         melints.push_back(seqtoimport);
      }
   }

   interestmelano.close();
   cout << "inputcoords finished" << endl;







   ofstream outf;
   for (ivseq iv=regints.begin();iv!=regints.end();iv++){
      Sequence seq=*iv;
      if (seq.species[0] && seq.nbtb>0){ 
         stringstream file;
         file << folder;
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

   ifstream fcoord(extract_args.input_arg);

   extract_cmdline_parser_free(&extract_args);

   return 1;
}		/* -----  end of function extract  ----- */
