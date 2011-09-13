/*
 * =====================================================================================
 *
 *       Filename:  genmot.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06.08.2011 13:05:22
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
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>

using namespace std;

#include "genmot_cmdline.h"
#include "const.hpp"
#include "random.hpp"
#include "motif.hpp"
#include "tree.hpp"
#include "distinfo.hpp"

genmot_args_info genmot_args;

vseq regtests;
vseq regints;


/** 
 * ===  FUNCTION  ======================================================================
 *         Name:  motscoreorder
 *  Description:  Compare motif overrepresentation
 * =====================================================================================
 */
   bool
motscoreorder ( Motif mot1, Motif mot2 )
{
   return mot1.pvalue<mot2.pvalue;
}		/* -----  end of function motscoreorder  ----- */

   bool
motchi2order ( Motif mot1, Motif mot2 )
{
   return mot1.scorepoiss<mot2.scorepoiss;
}		/* -----  end of function motscoreorder  ----- */

   double
distcv(vvd& mat1, vvd& mat2)
{
   double max(0);
   int col(0);
   for (ivvd iv=mat1.begin();iv!=mat1.end();iv++){
      int line(0);
      for (ivd il=(*iv).begin();il!=(*iv).end();il++){
         double diff;
         diff=fabs((*il)-mat2[col][line]);
         if (diff>max) max=diff;
         line++;
      }
      col++;
   }
   return max;
}

   void
motiftomat(vint & seq,Motif & mot)
{
   vd dum(4,0);
   mot.matrice=vvd(width,dum);
   vvd::iterator imat=mot.matrice.begin();
   for (vint::iterator iseq=seq.begin();iseq!=seq.end();iseq++){
      int base=*iseq;
      double a,t,c,g;
      a=t=alpha/(1+2*alpha+2*beta);
      c=g=beta/(1+2*alpha+2*beta);
      if (base==0){
         a=(1+alpha)/(1+2*alpha+2*beta);
      }
      else if (base==1){
         t=(1+alpha)/(1+2*alpha+2*beta);
      }
      else if (base==2){
         c=(1+beta)/(1+2*alpha+2*beta);
      }
      else if (base==3){
         g=(1+beta)/(1+2*alpha+2*beta);
      }
      (*imat)[0]=log(a/conca);
      (*imat)[1]=log(t/conct);
      (*imat)[2]=log(c/concc);
      (*imat)[3]=log(g/concg);
      imat++;
   }
}

   void
seqanalysis(Sequence & currseq,vmot & genmots)
{
   unsigned int i=0;
//      for (int j=0;j<nbspecies;j++){
//      cout << currseq.iseqs[j] << endl;
//      }
   for (vint::iterator istr=currseq.iseqs[0].begin();istr!=currseq.iseqs[0].end()-width+1;istr++){
      //cout << "\r" << i+1 << "/" << currseq.iseqs[0].size()-width+1 ; 
      vint bs(istr,istr+width);
     cout << i << " " << bs << endl;
      if (compN(bs)>0) continue;
      Motif currmot;
      currmot.bsinit=vinttostring(bs);
      currmot.seqinit=currseq.name;
      currmot.pos=i;
      motiftomat(bs,currmot);
      currmot.matricerevcomp=reversecomp(currmot.matrice);
      currmot.matprec=currmot.matrice;
      currmot.matprecrevcomp=currmot.matricerevcomp;
      vvd pmat=currmot.matprec;
      unsigned int nbconv(0);
      for (int nb=1;nb<=nbiter;nb++){
         double max=0.01;
         int iter(0);
         while(max>0){
            if (nb>2) currmot.matinit(scorethr2);
            else currmot.matinit(scorethr);
            if (currmot.nbmot<1) break;
            
            currmot.compprec();
            max=distcv(currmot.matprec,pmat);
            pmat=currmot.matprec;
            iter++;
            if (iter==20) break;
            nbconv++;
         }
         if (nb==1){
            currmot.corrprec();
            pmat=currmot.matprec;
            currmot.matprecrevcomp=reversecomp(currmot.matprec);
         }
      }
      //cout << currmot.nbmot << " " ; 
      //cout.flush();

      currmot.matinit(scorethr2);
      if (currmot.nbmot>2){
         currmot.pvaluecomp();
         //currmot.display(streamfile);
         //		for (ivma ima=currmot.seqs.begin();ima!=currmot.seqs.end();ima++){
         //  cout << (*ima).alignseq[0] << endl;
         //		}
         //cout << currmot.matprec << endl;
         currmot.matprecrevcomp=reversecomp(currmot.matprec);
         currmot.matfreq=mattofreq(currmot.matprec);
         currmot.motscorethr=scorethr2;
         currmot.motwidth=width;
         genmots.push_back(currmot);
      }
      i++;
   }
   cout << endl;
}


//Loads align-file (fasta) or coord-file (name/chrom/start/stop)
   vseq
loadseqs()
{
   vseq seqs;

   if (genmot_args.align_file_given){

      ifstream inf;
      inf.open(genmot_args.align_file_arg);
      seqs=loadsequencesconserv(inf);
      inf.close();
   }
   else if (genmot_args.coord_file_given){

      ifstream inf;
      inf.open(genmot_args.coord_file_arg);

      ifstream align;
      if (species=="droso") align.open(DATA_PATH"/droso/align.dat");
      else if (species=="eutherian") align.open(DATA_PATH"/eutherian/align.dat");

      alignscoord=loadcoordconserv(align);

      align.close();

      vcoord coords;
      back_insert_iterator<vcoord> dest(coords);
      copy(iiscoord(inf),iiscoord(),dest);
      for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
         Sequence seqtoimport=coordtoseq(*ivc);
         //         cout << (*ivc).name << " " << (*ivc).start << endl;
         //         cout << seqtoimport.name << " " << seqtoimport.start << endl;
         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
            seqs.push_back(seqtoimport);
         }
         //         for (int i=0;i<nbspecies;i++){
         //         cout << seqtoimport.iseqs[i] << endl; 
         //         }
         //         exit(9);
      }
      inf.close();
   }
   else {
      cout << "Error in loadseqs: please give a coord/alignment file. Exiting..." << endl;
      exit(1);
   }
   return seqs;
}

void
args_init()
{
   if (!strcmp(genmot_args.species_arg,"droso")){
      species="droso";
      nbspecies=12;
      conca=0.3; 
   } else if (!strcmp(genmot_args.species_arg,"eutherian")){
      species="eutherian";
      nbspecies=12;
      conca=0.263; 
   }
   concc=0.5-conca;
   conct=conca;
   concg=concc;

   width=genmot_args.width_arg;
   
   // *** It would be nice to set the threshold by bp, in bits.
   scorethr2=width*genmot_args.threshold_arg/10;
   scorethr=width*(scorethr2-1.0)/10;
   scorethrcons=width*(scorethr2-1.0)/10;

   evolutionary_model=genmot_args.evolutionary_model_arg;

   neighbext=genmot_args.neighbext_arg;

}


/** 
 * ===  FUNCTION  ======================================================================
 *         Name:  cmd_genmot
 *  Description:  Motif generation
 * =====================================================================================
 */
int
cmd_genmot(int argc, char **argv)
{

   if ( genmot_cmdline_parser(argc,argv, & genmot_args)!=0)
      exit(1);

   args_init();

   rnginit();

   //   printconfig(); *** to be written so that one can rerun exacltly the same instance

   compalpha();
   //   printpriorsandthrs(); *** to be written
   cout << alpha << endl;


   cout << "Loading background set..." << endl;
   ifstream backreg;

   // files put manually: should be created in the early building process.
   if (species=="droso") backreg.open(DATA_PATH"/droso/background2000.fa");
   else if (species=="eutherian") backreg.open(DATA_PATH"/eutherian/background2000.fa");
   regtests=loadsequences(backreg);
   backreg.close();
   cout << "Background set size : " << regtests.size() << endl;

   cout << "Loading training set..." << endl;
   regints=loadseqs();
   cout << "Training set size : " << regints.size() << endl;

   inittreedist();

   cout << "Generating motifs..." << endl;
   vmot genmots;
   for (vseq::iterator iseq=regints.begin();iseq!=regints.end();iseq++){
      cout << (*iseq).name << endl;
      seqanalysis(*iseq,genmots);
      cout << endl;
   }
   cout << "Sorting motifs..." << endl;
   // Sort on chi2
   sort(genmots.begin(),genmots.end(),motchi2order);
   // find 3rd quartile
   int shiftchi2=genmots.size()*3/4;
   genmots.erase(genmots.begin()+shiftchi2,genmots.end());
   // Sort 
   sort(genmots.begin(),genmots.end(),motscoreorder);

   // Cluster...
   cout << "Clustering motifs..." << endl;
   compmotsdist(genmots);
   //
   ofstream motmeldb("motifs.txt");
   for ( ivmot ivm=genmots.begin();ivm!=genmots.end();ivm++ ) {
      if ( ivm->check ) {
         ivm->display(motmeldb);
      }
   }
   motmeldb.close();

   gsl_rng_free(gslran);
   genmot_cmdline_parser_free(&genmot_args);
   cout << "exit normally" << endl;
   return 1;
}		/* -----  end of function genmot  ----- */


