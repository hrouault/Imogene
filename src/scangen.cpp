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
 * =====================================================================================
 *
 *       Filename:  scangen.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06.08.2011 13:06:06
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
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>

using namespace std;

#include "scangen_cmdline.h"
#include "scangen.hpp"
#include "const.hpp"
#include "motif.hpp"
#include "tree.hpp"
#include "sequence.hpp"

vmot motsdef;
vinst allinstances;
vginst potregs;
vvginst groupedinst;
vstring phenos;
vstring gbacks;
vTSS TSSall;

scangen_args_info scangen_args;

   void
loadannots()
{
   ifstream annots;
   annots.open(scangen_args.phenotype_arg);
   if (annots.fail()){
      cerr << "Cannot open phenotype file: " << strerror(errno) << endl;
      exit(EXIT_FAILURE);
   }

   back_insert_iterator<vstring> dest(phenos);
   copy(iisstring(annots),iisstring(),dest);

   annots.close();

   ifstream glist;
   if (species=="droso"){
      cout << "Reading droso genes list..." << endl;
      glist.open( (scangen_datapath+"/droso/annot/genelist.dat").c_str() );
   } else if (species=="eutherian"){
      cout << "Reading eutherian genes list..." << endl;
      glist.open( (scangen_datapath+"/eutherian/annot/genelist.dat").c_str() );
   }
   if (glist.fail()){
      cerr << "Cannot open gene list file: " << strerror(errno) << endl;
      exit(EXIT_FAILURE);
   }

   back_insert_iterator<vstring> destg(gbacks);
   copy(iisstring(glist),iisstring(),destg);

   glist.close();

   for (ivstring ivs=phenos.begin();ivs!=phenos.end();ivs++){
      for (ivstring ivs2=gbacks.begin();ivs2!=gbacks.end();ivs2++){
         if (*ivs2==*ivs){
            gbacks.erase(ivs2);
            break;
         }
      }
   }

   //   for (ivstring ivs=phenos.begin();ivs!=phenos.end();ivs++){
   //      cout << *ivs << endl;
   //   }
}

   void
scanseqforconsinstances(Sequence &seq,vmot & mots)
{
   seq.instances.clear();
   for (ivmot im=mots.begin();im!=mots.end();im++){
      if (seq.iseqs[0].size()>im->motwidth){
         im->matinitforscanmots(seq);
      }
   }
   return;
}

   void
scanmots()
{
   ifstream align;
   if (species=="droso"){
      cout << "Reading droso alignments..." << endl;
      align.open( (scangen_datapath+"/droso/align.dat").c_str() );
   } else if (species=="eutherian"){
      cout << "Reading eutherian alignments..." << endl;
      align.open( (scangen_datapath+"/eutherian/align.dat").c_str() );
   }
   if (align.fail()){
      cerr << "Alignment file opening failed: " << strerror(errno) << endl;
      exit(EXIT_FAILURE);
   }

   alignscoord=loadcoordconserv(align);

   align.close();

   string pchrom(""); // for display purpose
   unsigned int totlen(0),totlentb(0);
   //ivcoord ivcstart=alignscoord.begin()+100;
   //ivcoord ivcstop=alignscoord.begin()+150;
   ivcoord ivcstart=alignscoord.begin();
   ivcoord ivcstop=alignscoord.end();
   for (ivcoord ivc=ivcstart;ivc!=ivcstop;ivc++){
      Sequence seq;
      seq=coordtoseq(*ivc);
      string chrom=chromfromint(seq.chrom);
      if (chrom!=pchrom){
         if (seq.chrom!=0) cout << "\n";
         cout << "chromosome " << chrom << "\n";
         pchrom=chrom;
      }
      cout << "\r" << (double)(ivc-ivcstart+1)/(ivcstop-ivcstart)*100 << "%\t";
      cout.flush();

      scanseqforconsinstances(seq,motsdef);

      totlen+=seq.nbtb+seq.nbN;
      totlentb+=seq.nbtb;
   }
   cout << "\n";

   for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      allinstances.insert(allinstances.end(),ivm->refinstances_short.begin(),ivm->refinstances_short.end());
   }

   sort(allinstances.begin(),allinstances.end());

   //   for (ivinst ivi=allinstances.begin();ivi!=allinstances.end();ivi++){
   //      cout << *ivi;
   //   }

   cout << "Total length of the alignement: " << (double)totlen/1000000 << " Mb, including " 
      << (double)totlentb/1000000 << " unmasked Mb" << endl;

   return;
}

   void
compgroupedinst()
{
   vginst dumvginst;
   groupedinst=vvginst(nbchrom,dumvginst);

   for (ivinst ivi=allinstances.begin();ivi!=allinstances.end();ivi++){
      ivi->isassigned=0;
   }

   // TSS importation
   ifstream annots;
   if (species=="droso"){
      cout << "Reading droso TSS annot..." << endl;
      annots.open( (scangen_datapath+"/droso/annot/TSS-coord.dat").c_str() );
   } else if (species=="eutherian"){
      cout << "Reading eutherian TSS annot..." << endl;
      annots.open( (scangen_datapath+"/eutherian/annot/TSS-coord.dat").c_str() );
   }
   if (annots.fail()){
      cerr << "TSS annotation file opening failed: " << strerror(errno) << endl;
      exit(EXIT_FAILURE);
   }
   importTSS(TSSall,annots);
   annots.close();

   // in the following, allinstances supposed to be sorted (done in scanmots)

   cout << "Defining CRMs and assigning to TSS in annotextent region..." << endl;
   int counter=1;
   for (ivinst ivi=allinstances.begin();ivi!=allinstances.end();ivi++){
      //cout << "\r" << counter << " " << allinstances.size();
      cout.flush();
      counter++;
      if (ivi->isassigned == 0){
         GroupInstance ginst;
         ginst.chrom=ivi->chrom;
         ginst.start=ivi->coord;
         ginst.stop=ivi->coord+width-1;
         ginst.instances.push_back(*ivi);
         ivi->isassigned=1;
         if (ivi!=allinstances.end()-1){
            for (ivinst ivi2=ivi+1;ivi2!=allinstances.end();ivi2++){
               if ( (ivi2->chrom == ivi->chrom) && ( ivi2->coord - ivi->coord < scanwidth - width +1 ) ){
                  ginst.instances.push_back(*ivi2);
                  ginst.stop=ivi2->coord+width-1;
                  ivi2->isassigned=1;
               }
               else break;
            }
         }
         int instlen=ginst.stop-ginst.start+1;
         ginst.start-=ceil((double)(scanwidth-instlen)/2);
         ginst.stop+=floor((double)(scanwidth-instlen)/2);

         for (ivinst iv=ginst.instances.begin();iv!=ginst.instances.end();iv++){
            ginst.nbmots[iv->motindex]++;
            ginst.totmots++;
         }

         ginst.compscore(motsdef,nbmots_for_score);

         for (ivTSS ivt=TSSall.begin();ivt!=TSSall.end();ivt++){
            if (ivt->chrom == ginst.chrom && abs(ivt->coord-(ginst.start+ginst.stop)/2)<=annotextent){
               ginst.TSSs.push_back(*ivt);
            }
         }


         //         cout << ginst;
         //         cout << "->";
         //         for (ivTSS ivt=ginst.TSSs.begin();ivt!=ginst.TSSs.end();ivt++){
         //            cout << ivt->gene << ",";
         //         }
         //         cout << "\n\n";

         groupedinst[ivi->chrom].push_back(ginst);
      }
   }
   return;
}

   void
compTSSs()
{  
   // assign nearest TSS to CRMs
   for (unsigned int i=0;i<nbchrom;i++){
      for (ivginst ivg=groupedinst[i].begin();ivg!=groupedinst[i].end();ivg++){
         (*ivg).compbestannot();
      }
   }
   return;
}

   void 
comppheno()
{
   // attribute phenotype to CRMs
   for (unsigned int i=0;i<nbchrom;i++){
      for (ivginst ivg=groupedinst[i].begin();ivg!=groupedinst[i].end();ivg++){
         (*ivg).isdiscarded();
         if (!(*ivg).discarded){
            for (ivstring ivs=phenos.begin();ivs!=phenos.end();ivs++){
               if ((*ivg).besttss.gene==*ivs){
                  (*ivg).goodpheno=1;
                  break;
               } else {
                  (*ivg).goodpheno=0;
               }
            }
         } else {
            (*ivg).goodpheno=0;
         }
      }
   }
   return;
}

   void 
outputresults()
{
   vginst finginst;
   for (unsigned int i=0;i<nbchrom;i++){
      for (ivginst ivg=groupedinst[i].begin();ivg!=groupedinst[i].end();ivg++){
         finginst.push_back(*ivg);
      }
   }
   cout << "Sorting" << endl;
   sort(finginst.begin(),finginst.end());

   cout << "Displaying results" << endl;
   stringstream filename;
   filename << "result" << nbmots_for_score << ".dat";
   ofstream res(filename.str().c_str());
   if (res.fail()){
      cerr << "Cannot open scangen result file for writing: " << strerror(errno) << endl;
      exit(EXIT_FAILURE);
   }
   for (ivginst ivg=finginst.begin();ivg!=finginst.end();ivg++){
      res << (*ivg).score << " " << chromfromint((*ivg).chrom) << ":" << (*ivg).start << ".." << (*ivg).stop << " ";
      res << (*ivg).besttss.gene << " ";
      for (ivTSS ivt=(*ivg).TSSs.begin();ivt!=(*ivg).TSSs.end();ivt++){
         res << (*ivt).gene << ";";
      }
      res << " ";
      for (ivint ivi=(*ivg).nbmots.begin();ivi!=(*ivg).nbmots.end();ivi++){
         res << *ivi << " ";
      }
      res << "\n";
   }
   res.close();

   stringstream filenameseqs;
   filenameseqs << "mots" << nbmots_for_score << ".dat";
   ofstream seqs(filenameseqs.str().c_str());
   if (seqs.fail()){
      cerr << "Cannot open sequence file for writing: " << strerror(errno) << endl;
      exit(EXIT_FAILURE);
   }
   seqs << "chrom" << " ";
   seqs << "start" << " ";
   seqs << "stop" << " ";
   seqs << "motif_index:position(0-based)" << "\n";
   unsigned int count=0;
   for (ivginst ivg=finginst.begin();ivg!=finginst.end();ivg++){
      seqs << chromfromint(ivg->chrom) << " ";
      seqs << ivg->start << " ";
      seqs << ivg->stop << " ";
      for (ivinst ivi=ivg->instances.begin();ivi!=ivg->instances.end();ivi++){
         seqs << ivi->motindex+1 << ":" << ivi->coord-ivg->start << " ";
      }
      seqs << "\n";
      count++;
      if (count>100) break;
   }
   seqs << "\n";
   seqs.close();

   if (scangen_args.phenotype_given){
      stringstream filehist;
      filehist << "hist" << nbmots_for_score << ".dat";
      ofstream histo(filehist.str().c_str());
      if (histo.fail()){
         cerr << "Cannot open histogram file for writing: " << strerror(errno) << endl;
         exit(EXIT_FAILURE);
      }
      // positions of phenotype genes in the reverse-sorted result file
      displayhist(finginst,histo);
      histo.close();
      if (scangen_args.print_histo_sets_given){
         stringstream filehist_back;
         filehist_back << "hist-back" << nbmots_for_score << ".dat";
         ofstream histo_back(filehist_back.str().c_str());
         if (histo_back.fail()){
            cerr << "Cannot open histogram background file for writing: " << strerror(errno) << endl;
            exit(EXIT_FAILURE);
         }
         // best CRM score for phenotype genes
         displayhist_set(finginst,gbacks,histo_back);
         histo_back.close();

         stringstream filehist_interest;
         filehist_interest << "hist-interest" << nbmots_for_score << ".dat";
         ofstream histo_interest(filehist_interest.str().c_str());
         if (histo_interest.fail()){
            cerr << "Cannot open histogram interest file for writing: " << strerror(errno) << endl;
            exit(EXIT_FAILURE);
         }
         // best CRM score for other genes
         displayhist_set(finginst,phenos,histo_interest);
         histo_interest.close();
      }
   }

   return;

}

   void
loadmotsforscangen()
{
   loadmots(scangen_args.motifs_arg,motsdef); 
   cout << "Loaded " << motsdef.size() << " motifs." << endl;
   if ( nbmots_for_score < motsdef.size() ) 
      motsdef.erase( motsdef.begin() + nbmots_for_score , motsdef.end() );
   cout << "Nb mots for score: " << nbmots_for_score  << endl;

   // *** It would be nice to set the threshold by bp, in bits.
   width=motsdef[0].bsinit.size();
   scorethr2=width*scangen_args.threshold_arg/10;
   scorethr=width*(scorethr2-1.0)/10;
   scorethrcons=width*(scorethr2-1.0)/10;
   cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;

   for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      ivm->motscorethr2=scorethr2;
      ivm->motscorethr=scorethr;
      ivm->motscorethrcons=scorethrcons;
   }
   return;
}

   void
scangen_args_init()
{
   if (!strcmp(scangen_args.species_arg,"droso")){
      species="droso";
      nbspecies=12;
      conca=0.3; 
      nbchrom=6;
      annotextent=10000;
   } else if (!strcmp(scangen_args.species_arg,"eutherian")){
      species="eutherian";
      nbspecies=12;
      conca=0.263;
      nbchrom=21; 
      annotextent=1000000;
   }
   concc=0.5-conca;
   conct=conca;
   concg=concc;

   scanwidth=scangen_args.scanwidth_arg;
   if (scangen_args.annotextent_given){
      annotextent=scangen_args.annotextent_arg;
   }

   nbmots_for_score=scangen_args.nbmots_arg;

   neighbext=scangen_args.neighbext_arg;

}

string scangen_datapath;

   int
cmd_scangen(int argc, char **argv)
{

   if ( scangen_cmdline_parser(argc,argv, & scangen_args)!=0)
      exit(EXIT_FAILURE);

   scangen_args_init();

   const char * imo_scangen_datapath = getenv( "IMOGENE_DATA" );
   if (imo_scangen_datapath==NULL){
      scangen_datapath = DATA_PATH;
   } else {
      scangen_datapath=imo_scangen_datapath;
   }

   cout << "annotextent=" << annotextent << endl;

   cout << "Loading Motifs" << endl;
   loadmotsforscangen();

   if (scangen_args.phenotype_given){
      cout << "Loading phenotypes" << endl;
      loadannots();
   }
   else{
      cout << "No phenotype file given" << endl;
   }

   cout << "Scanning genome for conserved motif instances" << endl;
   scanmots();

   cout << "Defining instances" << endl;
   compgroupedinst();

   cout << "Assign CRMs to nearest gene..." << endl;
   compTSSs();

   if (scangen_args.phenotype_given){
      cout << "Attribute phenotype to CRMs..." << endl;
      comppheno();
   }

   cout << "Output result..." << endl;
   outputresults();

   return EXIT_SUCCESS;
}
