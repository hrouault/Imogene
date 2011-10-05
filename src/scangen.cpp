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

   back_insert_iterator<vstring> dest(phenos);
   copy(iisstring(annots),iisstring(),dest);

   annots.close();

   ifstream glist;
   if (species=="droso"){
      cout << "Reading droso genes list..." << endl;
      glist.open(DATA_PATH"/droso/annot/genelist.dat");
   } else if (species=="eutherian"){
      cout << "Reading eutherian genes list..." << endl;
      glist.open(DATA_PATH"/eutherian/annot/genelist.dat");
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
      align.open(DATA_PATH"/droso/align.dat");
   } else if (species=="eutherian"){
      cout << "Reading eutherian alignments..." << endl;
      align.open(DATA_PATH"/eutherian/align.dat");
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
      annots.open(DATA_PATH"/droso/annot/TSS-coord.dat");
   } else if (species=="eutherian"){
      cout << "Reading eutherian TSS annot..." << endl;
      annots.open(DATA_PATH"/eutherian/annot/TSS-coord.dat");
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
               if ( (ivi2->chrom == ivi->chrom) && (ivi2->coord - ivi->coord < scanwidth - width +1 ) ){
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
      // positions of phenotype genes in the reverse-sorted result file
      displayhist(finginst,histo);
      histo.close();
      if (scangen_args.print_histo_sets_given){
         stringstream filehist_back;
         filehist_back << "hist-back" << nbmots_for_score << ".dat";
         ofstream histo_back(filehist_back.str().c_str());
         // best CRM score for phenotype genes
         displayhist_set(finginst,gbacks,histo_back);
         histo_back.close();

         stringstream filehist_interest;
         filehist_interest << "hist-interest" << nbmots_for_score << ".dat";
         ofstream histo_interest(filehist_interest.str().c_str());
         // best CRM score for other genes
         displayhist_set(finginst,phenos,histo_interest);
         histo_interest.close();
      }
   }

   return;

}


   void
outputresultsfornbmots(unsigned int nbmots_score)
{
   cout << "Shuffling" << endl;
   random_shuffle(potregs.begin(),potregs.end());
   for (ivginst ivg=potregs.begin();ivg!=potregs.end();ivg++){
      (*ivg).compscore(motsdef,nbmots_score);
   }
   cout << "Sorting" << endl;
   sort(potregs.begin(),potregs.end());
   cout << "Discarding" << endl;
   vstring gnames;
   vginst potregs_def;

   // Keep best scoring enhancer per gene
   // OR
   // Keep all enhancers per gene, removing overlapping low score ones
   if (scangen_args.discard_on_gene_names_given){
      cout << "discard on gene" << endl;
      int test=0;
      for (ivginst ivg=potregs.begin();ivg!=potregs.end();ivg++){
         test=0;
         for (ivstring ivs=gnames.begin();ivs!=gnames.end();ivs++){
            if ((*ivg).besttss.gene==*ivs){
               (*ivg).discarded=1;
               test=1;
               break;
            }
         }
         if (test==0){
            gnames.push_back((*ivg).besttss.gene);
         }
      }
   } 
   else {
      for (ivginst ivg=potregs.begin();ivg!=potregs.end();ivg++){
         int test=0;
         for (ivginst ivg2=potregs_def.begin();ivg2!=potregs_def.end();ivg2++){
            if ((*ivg).distance(*ivg2)<scanwidth){
               test=1;
               break;
            }
         }
         if (test==0){
            potregs_def.push_back(*ivg);
         }
      }
   }
   cout << "Displaying" << endl;
   stringstream filename;
   filename << "result" << nbmots_score << ".dat";
   ofstream res(filename.str().c_str());
   for (ivginst ivg=potregs_def.begin();ivg!=potregs_def.end();ivg++){
      if (!(*ivg).discarded){
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
         //         for (ivinst ivi=(*ivg).instances.begin();ivi!=(*ivg).instances.end();ivi++){
         //            cout << "\t" << chromfromint((*ivi).chrom) << ":" << (*ivi).coord << "\n";
         //         }
      }
   }
   res << "\n";
   res.close();

   //	stringstream filenameseqs;
   //	filenameseqs << "seqs-" << nbmots_score << ".dat";
   //	ofstream seqs(filenameseqs.str().c_str());
   //	unsigned int count=0;
   //	for (ivginst ivg=potregs_def.begin();ivg!=potregs_def.end();ivg++){
   //		if (!(*ivg).discarded){
   //			unsigned int mini=scanwidth;
   //			unsigned int max=0;
   //			seqs << "chrom : " << (*ivg).chrom << "\n";
   //			seqs << "start : " << (*ivg).start << "\n";
   //			seqs << "stop : " << (*ivg).stop << "\n";
   //			for (ivinst ivi=(*ivg).instances.begin();ivi!=(*ivg).instances.end();ivi++){
   //				if ((*ivi).coord<mini) mini=(*ivi).coord;
   //				if ((*ivi).coord>max) max=(*ivi).coord;
   //			}
   //			seqs << "mini : " << mini << "\n";
   //			seqs << "max : " << max << "\n";
   //			unsigned int start=(mini+max)/2-scanwidth/2;
   //			unsigned int stop=(mini+max)/2+scanwidth/2;
   //			for (ivinst ivi=(*ivg).instances.begin();ivi!=(*ivg).instances.end();ivi++){
   //				seqs << "motif " << (*ivi).motindex << " at position " << (*ivi).coord-start << "\n";;
   //			}
   //			seqs << ">seq_" << count << " annotated " << (*ivg).besttss.gene;
   //			seqs << " ";
   //			seqs << "\n";
   //			seqs << chromints[(*ivg).chrom].seq.substr(start,stop-start) << "\n";
   //			count++;
   //			if (count>100) break;
   //		}
   //	}
   //	seqs << "\n";
   //	seqs.close();

   if (scangen_args.phenotype_given){
      stringstream filehist;
      filehist << "hist" << nbmots_score << ".dat";
      ofstream histo(filehist.str().c_str());
      displayhist(potregs_def,histo);
      histo.close();
      if (scangen_args.print_histo_sets_given){
         stringstream filehist_back;
         filehist_back << "hist-back" << nbmots_score << ".dat";
         ofstream histo_back(filehist_back.str().c_str());
         displayhist_set(potregs_def,gbacks,histo_back);
         histo_back.close();

         stringstream filehist_interest;
         filehist_interest << "hist-interest" << nbmots_score << ".dat";
         ofstream histo_interest(filehist_interest.str().c_str());
         displayhist_set(potregs_def,phenos,histo_interest);
         histo_interest.close();
      }
   }

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

   int
cmd_scangen(int argc, char **argv)
{

   if ( scangen_cmdline_parser(argc,argv, & scangen_args)!=0)
      exit(1);

   scangen_args_init();

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
