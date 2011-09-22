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
vchrom chromints;
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
   if (species=="droso") glist.open("/home/rouault/these/sequence/genomes/genelist.dat");
   else if (species=="eutherian") glist.open("/home/santolin/these/files/mus/biomart/genelist-protein-coding+miRNA.dat");

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
   
   vint 
loadlengthchrom()
{

   vint lchr(nbchrom,0);

   if (species=="droso"){

      lchr[0]=22410834;//chr2L
      lchr[1]=20769785;//chr2R
      lchr[2]=23774897;//chr3L
      lchr[3]=27908053;//chr3R
      lchr[4]=1284640;//chr4
      lchr[5]=22227390;//chrX


   } else if (species=="eutherian"){

      lchr[0]=197195432;
      lchr[1]=181748087;
      lchr[2]=159599783;
      lchr[3]=155630120;
      lchr[4]=152537259;
      lchr[5]=149517037;
      lchr[6]=152524553;
      lchr[7]=131738871;
      lchr[8]=124076172;
      lchr[9]=129993255;
      lchr[10]=121843856;
      lchr[11]=121257530;
      lchr[12]=120284312;
      lchr[13]=125194864;
      lchr[14]=103494974;
      lchr[15]=98319150;
      lchr[16]=95272651;
      lchr[17]=90772031;
      lchr[18]=61342430;
      lchr[19]=166650296;//X
      lchr[20]=15902555;//Y
   }

   return lchr;

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
   //for (ivcoord ivc=alignscoord.begin();ivc!=alignscoord.end();ivc++){
   for (ivcoord ivc=alignscoord.begin();ivc!=alignscoord.begin()+10;ivc++){
      Sequence seq;
      seq=coordtoseq(*ivc);
      string chrom=chromfromint(seq.chrom);
      if (chrom!=pchrom){
         cout << "chromosome " << chrom << endl;
         pchrom=chrom;
      }

      scanseqforconsinstances(seq,motsdef);

      totlen+=seq.nbtb+seq.nbN;
      totlentb+=seq.nbtb;
   }
   
   for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      allinstances.insert(allinstances.end(),ivm->refinstances_short.begin(),ivm->refinstances_short.end());
   }

   sort(allinstances.begin(),allinstances.end());

   for (ivinst ivi=allinstances.begin();ivi!=allinstances.end();ivi++){
      cout << *ivi;
   }

   cout << "Total length of the alignement: " << (double)totlen/1000000 << " Mb, including " 
      << (double)totlentb/1000000 << " unmasked Mb" << endl;

   return;
}

   
   void
initgroupedinst()
{
   cout << "Loading length chroms..." << endl;
   lengthchrom=loadlengthchrom();

   // *** This step is the time-consuming one.
   //  We change it to a per motif search (see scanmots)
   vginst dumvginst;
   groupedinst=vvginst(nbchrom,dumvginst);
   for (unsigned int i=0;i<nbchrom;i++){
      for (unsigned int j=0;j<lengthchrom[i]/scanstep+1;j++){
         cout << "\r" << i << " " << j;
         cout.flush();
         groupedinst[i].push_back(GroupInstance(j*scanstep,j*scanstep+scanwidth,i));
      }
      cout << "\n";
   }
   exit(9);
   ifstream annots;
   if (species=="droso"){
      annots.open("/home/rouault/these/sequence/genomes/regres-wellform-all.dat");
   } else if (species=="eutherian"){
      annots.open("/home/santolin/these/files/mus/biomart/genes-n-strand-protein-coding+miRNA-no-MT-n-NT-TSS.dat");
   }
   importTSS(TSSall,annots);
   annots.close();
   // assign CRMs to genes in an annotextent region
   cout << "Assign CRMs to genes in annotextend..." << endl;
   for (ivTSS ivt=TSSall.begin();ivt!=TSSall.end();ivt++){
      int start=0;
      if (((*ivt).coord-annotextent-scanwidth/2)/scanstep+1>0) start=((*ivt).coord-annotextent-scanwidth/2)/scanstep;
      int stop=((*ivt).coord+annotextent-scanwidth/2)/scanstep+1;
      if (stop>lengthchrom[(*ivt).chrom]/scanstep+1) stop=lengthchrom[(*ivt).chrom]/scanstep+1;
      for (unsigned int i=start;i<stop;i++){
         groupedinst[(*ivt).chrom][i].TSSs.push_back(*ivt);
         //       cout << (*ivt).gene << endl;
      }
   }
   // assign nearest TSS to genes
   cout << "Assign CRMs to nearest gene..." << endl;
   for (unsigned int i=0;i<nbchrom;i++){
      for (ivginst ivg=groupedinst[i].begin();ivg!=groupedinst[i].end();ivg++){
         (*ivg).compbestannot();
      }
   }

   // attribute phenotype to CRMs
   cout << "Attribute phenotype to CRM..." << endl;
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
         }

         ginst.compscore(motsdef,nbmots_for_score);
   
         for (ivTSS ivt=TSSall.begin();ivt!=TSSall.end();ivt++){
            if (ivt->chrom == ginst.chrom && abs(ivt->coord-(ginst.start+ginst.stop)/2)<=annotextent){
               ginst.TSSs.push_back(*ivt);
            }
         }
   

         cout << ginst << "\n";
         for (ivTSS ivt=ginst.TSSs.begin();ivt!=ginst.TSSs.end();ivt++){
            cout << ivt->gene << " ";
         }
         cout << "\n";

         groupedinst[ivi->chrom].push_back(ginst);
      }
   }
   exit(9);
   // assign nearest TSS to genes
   cout << "Assign CRMs to nearest gene..." << endl;
   for (unsigned int i=0;i<nbchrom;i++){
      for (ivginst ivg=groupedinst[i].begin();ivg!=groupedinst[i].end();ivg++){
         (*ivg).compbestannot();
      }
   }

   // attribute phenotype to CRMs
   cout << "Attribute phenotype to CRM..." << endl;
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

}
   
   
   void
outputresults(unsigned int nbmots_score)
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
scanseqs(vstring & regs)
{
   vcoord pcoords;
   
   string pchrom("");
   for (ivstring is=regs.begin();is!=regs.end();is++){
      //cout << *is << endl;
      Sequence seq;
      ifstream fseq;
      fseq.open((*is).c_str());
      string fseqline;
      getline(fseq,fseqline);
      seq.name=fseqline;
      stringstream firstline(fseqline);
      string speciesname;
      firstline >> speciesname;
      string chrom;
      firstline >> chrom;
      if (chrom!=pchrom){
         cout << chrom << endl;
         pchrom=chrom;
      }
      seq.chrom=intfromchrom(chrom.substr(3));
      if (seq.chrom==-1) continue;
      //cout << "chromosome \"" << seq.chrom << "\"\n";
      firstline >> seq.start;
      //cout << "start \"" << seq.start << "\"\n";
      firstline >> seq.stop;
      seq.finame=*is;
      string dum;
      seq.species=vint(nbspecies,0);
      vint dumi;
      seq.iseqs=vvint(nbspecies,dumi);
      seq.imaps=vvint(nbspecies,dumi);
      seq.imapsinv=vvint(nbspecies,dumi);
      while (!fseq.eof()){
         //                cout << fseqline << endl;
         int seqnum=speciestonum(fseqline.substr(1,6));
         getline(fseq,fseqline);
         seq.imaps[seqnum]=alignedtomap(fseqline);
         seq.imapsinv[seqnum]=alignedtorevmap(fseqline);
         string seqwogap=remgaps(fseqline);
         if (seqwogap.size()>30){
            seq.species[seqnum]=1;
         }
         seq.iseqs[seqnum]=stringtoint(seqwogap);
         getline(fseq,fseqline);
      }
      seq.nbN=compN(seq.iseqs[0]);
      seq.nbtb=seq.iseqs[0].size()-seq.nbN;
      
      for (ivcoord ivc=pcoords.begin();ivc!=pcoords.end();ivc++){
         if (ivc->chrom==seq.chrom){ 
            if ((ivc->start>=seq.start && ivc->start<=seq.stop) ||
                  (ivc->stop>=seq.start && ivc->stop <=seq.stop)){
                  int maskstart=max(ivc->start,seq.start);
                  int maskstop=min(ivc->stop,seq.stop);
      //            cout << ivc->name << endl;
      //            cout << vinttostring(seq.iseqs[0]) << endl;
                  for (unsigned int base=maskstart-seq.start;base<=maskstop-seq.start;base++){
                     seq.iseqs[0][base]=4;
                  }
      //            cout << vinttostring(seq.iseqs[0]) << endl;
            }
         }
      }

      //check that ref species contains at least 30bp
      if (seq.species[0]){
         scoreseq(seq,motsdef);
      }

      fseq.close();
   }
   unsigned int i=0;

   cout << "Finding Instances" << endl;
   for (ivmot im=motsdef.begin();im!=motsdef.begin()+nbmots_for_score;im++){
      for (unsigned int j=0;j<nbchrom;j++){
         //           cout << "Instances size : " << (*im).instances.size() << "\n";
         sort((*im).instances[j].begin(),(*im).instances[j].end());
         for (ivinst iins=(*im).instances[j].begin();iins!=(*im).instances[j].end();iins++){
            Instance & curinst=*iins;
            int start=0;
            if ((curinst.coord-scanwidth)/scanstep+2>0) start=(curinst.coord-scanwidth)/scanstep+1;
            for (unsigned int k=start;k<curinst.coord/scanstep+1;k++){
               groupedinst[curinst.chrom][k].nbmots[i]++;
               groupedinst[curinst.chrom][k].totmots++;
               groupedinst[curinst.chrom][k].instances.push_back(curinst);
            }
         }
      }
      i++;
   }

   for (unsigned int j=0;j<nbchrom;j++){
      for (ivginst ivg=groupedinst[j].begin();ivg!=groupedinst[j].end();ivg++){
         if ((*ivg).totmots==0){
            (*ivg).discarded=1;
         }
      }
   }
   for (unsigned int j=0;j<nbchrom;j++){
      for (ivginst ivg=groupedinst[j].begin();ivg!=groupedinst[j].end();ivg++){
         if (!(*ivg).discarded){
            potregs.push_back(*ivg);
         }
      }
   }

   cout << "output results" << endl;
   for (unsigned int i=1;i<nbmots_for_score+1;i++){
      cout << "Motif " << i << endl;
      outputresults(i);
   }

}
   
   void
scangen_args_init()
{
   if (!strcmp(scangen_args.species_arg,"droso")){
      species="droso";
      nbspecies=12;
      conca=0.3; 
      nbchrom=6;
   } else if (!strcmp(scangen_args.species_arg,"eutherian")){
      species="eutherian";
      nbspecies=12;
      conca=0.263;
     nbchrom=21; 
   }
   concc=0.5-conca;
   conct=conca;
   concg=concc;
   
   scanwidth=scangen_args.scanwidth_arg;
   scanstep=scangen_args.scanstep_arg;
   annotextent=scangen_args.annotextent_arg;

   nbmots_for_score=scangen_args.nbmots_arg;
   
   
   neighbext=scangen_args.neighbext_arg;

}

int
cmd_scangen(int argc, char **argv)
{

   if ( scangen_cmdline_parser(argc,argv, & scangen_args)!=0)
      exit(1);

   scangen_args_init();
      
   cout << "Loading Motifs" << endl;
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


   //cout << "Loading phenotypes" << endl;
   //loadannots();
   
   cout << "Scanning genome for conserved motif instances" << endl;
   scanmots();

   cout << "Defining instances" << endl;
   //initgroupedinst();
   compgroupedinst();

   ifstream potregs;
   if (species=="droso") potregs.open("/home/santolin/these/files/droso/align/all/align-files.dat");
   else if (species=="eutherian") potregs.open("/home/santolin/these/files/mus/epo/align-files.dat");
   vstring regs;
   back_insert_iterator<vstring> dest(regs);
   copy(iisstring(potregs),iisstring(),dest);
   potregs.close();

   cout << "Scanning seqs" << endl;
   scanseqs(regs);


}
