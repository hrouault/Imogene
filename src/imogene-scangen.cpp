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
 */
#include<cmath>
#include<iostream>
#include<iomanip>
#include<vector>
#include<string>
#include<fstream>
#include<limits>
#include<sstream>
#include<algorithm>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include <numeric>

using namespace std;

#include "const.hpp"
#include "vectortypes.hpp"
#include "random.hpp"
#include "sequence.hpp"
#include "motif.hpp"
#include "tree.hpp"
#include "montecarlo.hpp"
#include "scangen.hpp"
#include ""

gengetopt_args_info args_info;
vmot motsdef;
vginst potregs;
vvginst groupedinst;
vstring phenos;
vstring gbacks;
vchrom chromints;

unsigned int sizepos,sizeneg;
unsigned int cutoff_for_combination=3;

// Print backrgound regions from intergenic coordinates (in vcds)
// TODO some cleanup.
   void
printbackreg(vcoord & vcds,string folder)
{
   ofstream outf("backreg.dat");
   unsigned int numback(0);
   unsigned int i(0);
   unsigned int backsize;
   if (args_info.size_given) backsize=args_info.size_arg;
   else backsize=2000;
   while (numback<10000){
      if (i>=vcds.size()) i=0;
      int range;
      range=vcds[i].stop-backsize-vcds[i].start;
      //cout << vcds[i];
      if (range<0){
         i++;
         //     cout << "range < 0" << endl;
         continue;
      }
      unsigned int newstart;
      unsigned int newpos;
      newpos= gsl_rng_uniform_int (gslran,range+1);
      newstart=vcds[i].start+newpos;
      Coordinate cdtmp;
      cdtmp.start=newstart;
      cdtmp.stop=newstart+backsize-1;
      cdtmp.chrom=vcds[i].chrom;
      cdtmp.name=vcds[i].name;

      Sequence & s=seq;
      for (int k=0;k<nbspecies;k++){
         if (s.species[k]){
            if (k==0) outf << ">" << numtospecies(k) << " " <<
               "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
            else  outf << ">" << numtospecies(k) << endl;
            outf << s.seqsrealigned[k] << endl;
         }; 
      }
      outf.close();
      numback++;
   }
   Sequence seq=coordtoseq(cdtmp);
   int nbfr=0;
   bool bc=0;
   if (species==1){
      if (seq.species[5]) nbfr++; 
      if (seq.species[6] || seq.species[7]) nbfr++;
      if (seq.species[8]) nbfr++;
      if (seq.species[9] || seq.species[10] || seq.species[11]) nbfr++;
      if (nbfr>1) bc=1;
   }
   else if (species==2){
      if (seq.species[2] || seq.species[3] || seq.species[4] || seq.species[5]) nbfr++;
      if (seq.species[6] || seq.species[7]) nbfr++;
      if (seq.species[8]) nbfr++;
      if (seq.species[9]) nbfr++;
      if (nbfr>1) bc=1;
   }
   if (seq.species[0] && bc && seq.nbtb>backsize/2){
      cout << "->" << numback << ". " <<cdtmp;
      stringstream os;
      os << numback;
      os >> seq.name;
      string filename=folder;
      filename+=seq.name;
      filename.append(".fa");
      outf.open(filename.c_str());
      Sequence & s=seq;
      for (int k=0;k<nbspecies;k++){
         if (s.species[k]){
            if (k==0) outf << ">" << numtospecies(k) << " " <<
               "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
            else  outf << ">" << numtospecies(k) << endl;
            outf << s.seqsrealigned[k] << endl;
         }; 
      }
      outf.close();
      numback++;
   }
   else {
      //cout << "sequence not conserved or too much masked" << endl;
   }
   cout.flush();
   i++;
}
outf.close();
}

// Print backrgound regions from intergenic coordinates (in vcds), given coordinates (in regints)
// TODO some cleanup
void
printbackregwcoords(vcoord & vcds,string folder)
{

   //vcds is a vector of random shuffled intergenic regions
   //regints are our interest sequences
   //
   ofstream outf;
   unsigned int i(0);
   unsigned int numback(1);

   unsigned int repeat=10;
   if (args_info.numrepeat_given) repeat=args_info.numrepeat_arg;

   for (int j=1;j<=repeat;j++){
      for (ivseq ivs=regints.begin();ivs!=regints.end();ivs++){

         //size of our interest sequence
         unsigned int backsize=ivs->stop-ivs->start+1;
         unsigned int state=0;

         //if we find a corresponding background region, then state=1;
         while (state==0){

            //if we have searched all intergenic sequences we loop on the first
            if (i==vcds.size()-1) i=0;

            int range;
            range=vcds[i].stop-backsize-vcds[i].start+1;

            //intergenic size has to be superior to interest sequence
            if (range<1){
               i++;
               continue;
            }

            //randomly choose a region of a same size than interest in the intergenic region
            unsigned int newstart;
            unsigned int newpos;
            newpos= gsl_rng_uniform_int (gslran,range);
            newstart=vcds[i].start+newpos;
            Coordinate cdtmp;
            cdtmp.start=newstart;
            cdtmp.stop=newstart+backsize-1;
            cdtmp.chrom=vcds[i].chrom;
            cdtmp.name=vcds[i].name;

            //extract the alignment
            Sequence seq=coordtoseq(cdtmp);
            int nbfr=0;
            bool bc=0;

            if (species==1){
               if (seq.species[5]) nbfr++; 
               if (seq.species[6] || seq.species[7]) nbfr++;
               if (seq.species[8]) nbfr++;
               if (seq.species[9] || seq.species[10] || seq.species[11]) nbfr++;
               if (nbfr>1) bc=1;
            }
            else if (species==2){
               if (seq.species[2] || seq.species[3] || seq.species[4] || seq.species[5]) nbfr++;
               if (seq.species[6] || seq.species[7]) nbfr++;
               if (seq.species[8]) nbfr++;
               if (seq.species[9]) nbfr++;
               if (nbfr>1) bc=1;
            }


            //we want mus to be present, the sequence to be conserved, and at least half of bases to be unmasked
            if (seq.species[0] && bc && seq.nbtb>backsize/2){
               cout << "->" << numback << ". " <<cdtmp;
               cout.flush();
               stringstream os;
               os << numback;
               os >> seq.name;
               string filename=folder;
               filename+=seq.name;
               filename.append(".fa");
               outf.open(filename.c_str());
               Sequence & s=seq;
               for (int k=0;k<nbspecies;k++){
                  if (s.species[k]){
                     if (k==0) outf << ">" << numtospecies(k) << " " <<
                        "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
                     else  outf << ">" << numtospecies(k) << endl;
                     outf << s.seqsrealigned[k] << endl;
                  }; 
               }
               outf.close();
               numback++;
               state=1;
            }
            i++;
         }

      }
   }

   outf.close();
}


   void
scanseqforinstancesnmask(Sequence &seq,vmot & mots)
{
   // We mask a tmp sequence
   Sequence seqtomask=seq;
   for (ivmot im=mots.begin();im!=min(mots.end(),mots.begin()+nbmots_for_score);im++){
      im->findinstancesnmask(seqtomask);
   }
   seq.instances=seqtomask.instances;
   return;
}

   void
scanseqsforinstancesnmask(vseq & align,vmot & mots)
{     
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
      scanseqforinstancesnmask(*ivs,mots);
   }
   return;
}

   void
scanseqforinstances(Sequence &seq,vmot & mots)
{
   seq.instances.clear();
   for (ivmot im=mots.begin();im!=min(mots.end(),mots.begin()+nbmots_for_score);im++){
      im->findinstances(seq);
   }
   return;
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
scanseqsforinstances(vseq & align,vmot & mots)
{     
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
      scanseqforinstances(*ivs,mots);
   }
   return;
}

   void
scanseqsforinstances(vseq & align,Motif & mot)
{     
   vmot mots;
   mots.push_back(mot);
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
      scanseqforinstances(*ivs,mots);
   }
   return;
}

   void
scanseq(Sequence &seq,vmot & mots)
{
   //cout << seq.name << " " << seq.finame << "\n";
   int nbmot=0;
   vd nbcorr;
   double nmcorr=0;
   unsigned int moti=0;
   for (ivmot im=mots.begin();im!=mots.begin()+nbmots_for_score;im++){//im!=mots.end();im++){}
      int nm=0;
      if (args_info.disp_svg_given){
         nm=(*im).nbmatchnmaskforsvg(seq,moti);
      } else {
         //nm=(*im).nbmatchnmask(seq,moti);
         nm=(*im).nbmatchwomask(seq,moti);
      }
      nbmot=nm;

      if (args_info.weightmots_given){
         nmcorr+=nm*log((*im).lambdatrain/(*im).lambda);
      }
      else{
         nmcorr+=nm;
      }
      nbcorr.push_back(nmcorr);
      //if (nbmot!=0){cout << moti << "->" << nbmot << "\n";};
      moti++;
      //cout << seq.name << " " << nbmot << " " << endl;
//            for (ivinst ivs=(*im).instances[0].begin();ivs!=(*im).instances[0].end();ivs++){
//               cout << ivs->coord << " ";
//            }
//            cout << endl;
   }
   //      for (ivd iv=nbcorr.begin();iv!=nbcorr.end();iv++){
   //        cout << *iv << " ";
   //    }
}

   void
scanseqs(vstring & regs)
{
   vcoord pcoords;
   if (args_info.masktrain_given){

      cout << "(Masking training set)" << endl;
      ifstream inf;
      inf.open(args_info.masktrain_arg);
      back_insert_iterator<vcoord> pdest(pcoords);
      copy(iiscoord(inf),iiscoord(),pdest);
      inf.close();

   }
   
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
         scanseq(seq,motsdef);
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
   if (args_info.all_mots_only_given){
         outputresults(nbmots_for_score);
   }
   else{
      for (unsigned int i=1;i<nbmots_for_score+1;i++){
         cout << "Motif " << i << endl;
         outputresults(i);
      }
   }

}

   void
scanseqs(ifstream & list,vcoord & coords)
{
   
   vstring regs;
   back_insert_iterator<vstring> dest(regs);
   copy(iisstring(list),iisstring(),dest);

   string pchrom("");
   for (ivstring is=regs.begin();is!=regs.end();is++){
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
      firstline >> seq.start;
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
      
      scanseq(seq,motsdef);

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

   cout << "Discarding empty Instances" << endl;
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
   if (args_info.all_mots_only_given){
         outputresults(nbmots_for_score);
   }
   else{
      for (unsigned int i=1;i<nbmots_for_score+1;i++){
         cout << "Motif " << i << endl;
         outputresults(i);
      }
   }
}
   
   void
scanmots()
{
   ifstream potregs;
   if (species==1) potregs.open("/home/santolin/these/files/droso/align/all/align-files.dat");
   else if (species==2) potregs.open("/home/santolin/these/files/mus/epo/align-files.dat");
   //else if (species==2) potregs.open("/home/santolin/these/files/transfac/matrices/align-test.dat");
   vstring regs;
   back_insert_iterator<vstring> dest(regs);
   copy(iisstring(potregs),iisstring(),dest);
   potregs.close();

   string pchrom("");
   system("if ! test -d scanmots;then mkdir scanmots;fi;");      
   unsigned int totlen(0),totlentb(0);
   //random_shuffle(regs.begin(),regs.end());
   for (ivstring is=regs.begin();is!=regs.end();is++){
      Sequence seq;
      ifstream fseq;
      fseq.open((*is).c_str());
      string fseqline;
      getline(fseq,fseqline);
      seq.name=fseqline;
//cout << seq.name << endl;
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
      firstline >> seq.start;
      firstline >> seq.stop;
      seq.finame=*is;
      string dum;
      seq.species=vint(nbspecies,0);
      vint dumi;
      seq.iseqs=vvint(nbspecies,dumi);
      seq.imaps=vvint(nbspecies,dumi);
      seq.imapsinv=vvint(nbspecies,dumi);
      while (!fseq.eof()){
//cout << fseqline << endl;
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

      scanseqforconsinstances(seq,motsdef);

      fseq.close();
      
      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
         //sort(ivm->refinstances_short.begin(),ivm->refinstances_short.end());
         ostringstream outfs;
         outfs << "scanmots/" << ivm->name << "_" << ivm->motscorethr2 << ".dat";
         ofstream outf;
         if (totlen==0) outf.open(outfs.str().c_str());
         else outf.open(outfs.str().c_str(),ios::app);

         for (ivinst ivi=ivm->refinstances_short.begin();ivi!=ivm->refinstances_short.end();ivi++){
            outf << ivi->site << "\t";
            outf << chromfromint(ivi->chrom) << "\t";
            outf << ivi->coord << "\t";
            outf << ivi->sens << "\t";
            outf << ivi->score << endl;
         }
         outf.close();
         ivm->refinstances_short.clear();
      }

      totlen+=seq.nbtb+seq.nbN;
      totlentb+=seq.nbtb;
   }

   cout << "Total length of the alignement: " << totlen << " bp, including " << totlentb << " unmasked bp" << endl;

   return;
}

bool operator<(const Sequence & seqscr1,const Sequence & seqscr2)
{
   if (seqscr1.score == seqscr2.score)
   {
      return seqscr1.sign > seqscr2.sign; //for identical scores, we want the negatives first (for ROC consistency)
   }
   else
   {
      return seqscr1.score > seqscr2.score; 
   };
}

bool operator<(const Motif & mot1,const Motif & mot2)
{
   if (mot1.optauc == mot2.optauc)
   {
      return mot1.index < mot2.index;//the lower index first 
   }
   else
   {
      return mot1.optauc > mot2.optauc;//the best auc first
   };
}

   void
calcscore(Motif & im,Sequence & seqscr)//(vmot & mots)
{
   unsigned int ncons=im.nbmatchcons(seqscr);//seqscr.motis[im->index];
   //   unsigned int ncons=im.nbmatchmat(seqscr.seq);//seqscr.motis[im->index];
   //   unsigned int ncons=seqscr.motis[im.index];
   double lseq=seqscr.nbtb;
   double score=0;
   //cout << ncons << " ";
   //score= ncons;//*log(im.lambdatrain/im.lambda) + lseq*(im.lambda-im.lambdatrain);
   double lnj=0;
   for (int i=0;i<ncons;i++){
      lnj+=log(i+1);
   }
   seqscr.motis[im.index]=ncons;
   score=ncons;
   //score=ncons*log(im.lambdatrain/im.lambda);
   //score= lnj-ncons*log(im.lambda*lseq) + lseq*im.lambda;
   //   score= (double)ncons/lseq;//*log(im.lambdatrain/im.lambda) + lseq*(im.lambda-im.lambdatrain);
   seqscr.score=score;
   return;
   //cout << seqscr.score << endl;
}

   void
calcscore(vmot & vm,Sequence & seqscr)//(vmot & mots)
{

   Sequence seqtmp=seqscr;
   double lseq=seqscr.nbtb;
   double score=0;
   for (ivmot im=vm.begin();im!=min(vm.begin()+nbmots_for_score,vm.end());im++){
      width=im->bsinit.size();
      unsigned int ncons;
      if (args_info.maskforscore_given){
         ncons=im->nbmatchconsnmask(seqtmp);//seqscr.motis[im->index];
      } else{
         ncons=im->nbmatchcons(seqtmp);//seqscr.motis[im->index];
      }
      //   score+=ncons;
      //      double lnj=0;
      //      for (int i=0;i<ncons;i++){
      //         lnj+=log(i+1);
      //      }
      //   //   score+=lnj-ncons*log(im->lambda*lseq) + lseq*im->lambda;
      //      if (args_info.byaffinity_given){
      //         score+=im->scorematchcons(seqtmp);
      //      } else {
      //         score+=ncons;
      //      }
      if (args_info.weightmots_given){

         score+= ncons*log(im->lambdatrain/im->lambda);// + lseq*(im.lambda-im.lambdatrain);
      }
      else {
         score+= ncons;//*log(im->lambdatrain/im->lambda);// + lseq*(im.lambda-im.lambdatrain);
      }
      seqscr.motis[im->index]=ncons;
   }
   seqscr.score=score;
   return;
   //cout << seqscr.score << endl;
}

   void
calcscore(vmot & vm,vseq & vs)//(vmot & mots)
{
   for (ivseq ivs=vs.begin();ivs!=vs.end();ivs++){
      calcscore(vm,*ivs);
   }
   return;
   //cout << seqscr.score << endl;
}

// allncons is a vint where each int is a number of conserved motifs
   double
calcscore(vmot & vm,vint & allncons,int lseq)
{

   double score=0;
   for (ivmot im=vm.begin();im!=min(vm.begin()+nbmots_for_score,vm.end());im++){
      unsigned int ncons=allncons[im->index];
      double lnj=0;
      for (int i=0;i<ncons;i++){
         lnj+=log(i+1);
      }
      //score+=lnj-ncons*log(im->lambda*lseq) + lseq*im->lambda;
      //    score+=allncons[im->index];
      score+= allncons[im->index];//*log(im->lambdatrain/im->lambda);// + lseq*(im.lambda-im.lambdatrain);
   }
   return score;
}


svg::svg()
{
   xsize=800;
   ysize=600;
   xoffset=140;
   yoffset=65;
   pos=0;
};

   void
svginit(ofstream & svgfile, svg s)
{
   svgfile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
   svgfile << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n";
   svgfile << "<svg version=\"1.0\" id=\"Calque_1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" x=\"0px\" y=\"0px\"\n";
   svgfile << "	 width=\""<< s.xsize <<"px\" height=\""<< s.ysize <<"px\" viewBox=\"0 0 "<< s.xsize <<" "<< s.ysize <<
      "\" enable-background=\"new 0 0 "<< s.xsize <<" "<< s.ysize <<"\" xml:space=\"preserve\">\n";

   svgfile << "<line fill=\"none\" stroke=\"" << "black" << "\" stroke-width=\"" << 3 << "\" x1=\"" << 0 << "\" y1=\"" << 0 << "\" x2=\"" << 4 << "\" y2=\"" << 0 << "\"/>\n";
}

   void
svgdisplay(ofstream & svgfile,Sequence & seq, svg & s)
{
   double ytext=s.yoffset+300*s.pos;
   double xbegin=s.xoffset;
   double xend=xbegin+0.4*seq.imaps[0].size();

   svgfile << "<text transform=\"matrix(1 0 0 1 55.5 " << ytext << ")\" font-size=\"12\">" << seq.finame << "</text>\n";

   for (unsigned int i=0;i<nbspecies;i++){
      if (seq.species[i]){

         double yline=s.yoffset+20+300*s.pos+20*i;

         string dro;
         dro = seq.speciesname[i];
         svgfile << "<text transform=\"matrix(1 0 0 1 40 " << yline << ")\" font-size=\"12\">" << dro << "</text>\n";
         svgfile << "<line fill=\"none\" stroke=\"#000000\" stroke-width=\"0.5\" x1=\"" << xbegin << "\" y1=\"" << yline << "\" x2=\"" << xend << "\" y2=\"" << yline << "\"/>\n";

         vvint coords;
         vint curcoord;
         int seqorali=0;
         unsigned int p=0;
         for (istring is=seq.seqsrealigned[i].begin();is!=seq.seqsrealigned[i].end();is++){
            if (*is=='-'){
               if (seqorali==1){
                  curcoord.push_back(p);
                  coords.push_back(curcoord);
                  curcoord.clear();
                  seqorali=0;
               }
            } else {
               if (seqorali==0){
                  curcoord.push_back(p);
                  seqorali=1;
               }
               if (seqorali==1 && is==seq.seqsrealigned[i].end()-1){
                  curcoord.push_back(p);
                  coords.push_back(curcoord);
                  curcoord.clear();
               }
            }
            p++;
         }

         for (ivvint ivv=coords.begin();ivv!=coords.end();ivv++){
            svgfile << "<line fill=\"none\" stroke=\"#000000\" stroke-width=\"3\" x1=\"" << xbegin+0.4*(*ivv)[0] << 
               "\" y1=\"" << yline << "\" x2=\"" << xbegin+0.4*(*ivv)[1] << "\" y2=\"" << yline << "\"/>\n";
         }
      }
   }

   for (ivinstseq ivi=seq.instances.begin();ivi!=seq.instances.end();ivi++){
      int moti=(*ivi).motindex;
      if (moti<8){
         string color;
         if (moti==0){
            color="red";
         } else if (moti==1){
            color="blue";
         } else if (moti==2){
            color="green";
         } else if (moti==3){
            color="purple";
         } else if (moti==4){
            color="gray";
         } else if (moti==5){
            color="orange";
         } else if (moti==6){
            color="brown";
         } else if (moti==7){
            color="gold";
         }
         double xmot=xbegin+0.4*(*ivi).pos;
         double yline=s.yoffset+20+300*s.pos+20*(*ivi).species;
         string width;
         if ((*ivi).score>scorethr2) width="3";
         else width="1";
         svgfile << "<line fill=\"none\" stroke=\"" << color << "\" stroke-width=\"" << width << "\" x1=\"" << xmot << "\" y1=\"" << yline-5 << "\" x2=\"" << xmot << "\" y2=\"" << yline+5 << "\"/>\n";
         svgfile << "<text transform=\"matrix(1 0 0 1 " << xmot-2 << " " << yline-8 << ")\" font-size=\"8\">" << fixed << setprecision(1) << (*ivi).score << "</text>\n";
      }
   }
   s.pos++;
}


   void
svgclose(ofstream & svgfile)
{
   svgfile << "</svg>\n";
}

   void
scanseqsforsvg(vseq & align)
{

   system("if ! test -d display;then mkdir display;fi;");      
   system("if ! test -d display/svg;then mkdir display/svg;fi;");      

   for (ivseq is=align.begin();is!=align.end();is++){

      int xsize(0);
      int ysize(0);
      svg s;
      string filename("display/svg/");
      filename+=(*is).name;
      filename+=".svg";
      ofstream svgfile(filename.c_str());

      scanseq(*is,motsdef);

      //we set the size for the svg file
      xsize=s.xoffset+(int)(0.4*(*is).imaps[0].size());
      if (xsize>s.xsize) s.xsize=xsize;
      s.xsize+=s.xoffset;
      ysize=s.yoffset+340; 
      s.ysize=ysize+s.yoffset;

      svginit(svgfile,s);
      svgdisplay(svgfile,*is,s);
      svgclose(svgfile);

      //cout << s.xsize << " " << s.ysize << endl; 
   }
}

   void
outputresults(unsigned int nbmots_score)
{
   cout << "Shuffling" << endl;
   random_shuffle(potregs.begin(),potregs.end());
   for (ivginst ivg=potregs.begin();ivg!=potregs.end();ivg++){
      if (args_info.weightmots_given){
         (*ivg).compscoreweight(motsdef,nbmots_score);
      }
      else{
         (*ivg).compscore(motsdef,nbmots_score);
      }
   }
   cout << "Sorting" << endl;
   sort(potregs.begin(),potregs.end());
   cout << "Discarding" << endl;
   vstring gnames;
   vginst potregs_def;

   // Keep best scoring enhancer per gene
   // OR
   // Keep all enhancers per gene, removing overlapping low score ones
   if (args_info.discard_on_gene_names_given){
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

   if (args_info.phenotype_given){
      stringstream filehist;
      filehist << "hist" << nbmots_score << ".dat";
      ofstream histo(filehist.str().c_str());
      displayhist(potregs_def,histo);
      histo.close();
      if (args_info.print_histo_sets_given){
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


vseq regtests;
vseq regints;


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
seqanalysis(Sequence & currseq,ofstream & streamfile)
{
   unsigned int i=0;
//      for (int j=0;j<nbspecies;j++){
//      cout << currseq.iseqs[j] << endl;
//      }
   for (vint::iterator istr=currseq.iseqs[0].begin();istr!=currseq.iseqs[0].end()-width+1;istr++){
      //cout << "\r" << i+1 << "/" << currseq.iseqs[0].size()-width+1 ; 
      vint bs(istr,istr+width);
     //cout << i << " " << bs << endl;
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
            
            if (args_info.MCMCtest_given){
               cout << "TEST!" << endl;
               currmot.compprec_test();
               exit(1);
            }
            
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
         currmot.display(streamfile);
         //		for (ivma ima=currmot.seqs.begin();ima!=currmot.seqs.end();ima++){
         //  cout << (*ima).alignseq[0] << endl;
         //		}
         //cout << currmot.matprec << endl;
      }
      i++;
   }
   cout << endl;
}
   
   void
refseqanalysis(Sequence & currseq,ofstream & streamfile)
{
   unsigned int i=0;
   double scoresave=scorethr2;
   //   for (int j=0;j<nbspecies;j++){
   //   cout << currseq.iseqs[j] << endl;
   //   }
   for (vint::iterator istr=currseq.iseqs[0].begin();istr!=currseq.iseqs[0].end()-width+1;istr++){
      //cout << "\r" << i+1 << "/" << currseq.iseqs[0].size()-width+1 ; 
      vint bs(istr,istr+width);
      //cout << i << " " << bs << endl;
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
            if (nb>2) scorethr2=scoresave;
            else scorethr2=scorethr;
            
            currmot.matinit(scorethr2);
            if (currmot.nbmot<1) break;
            
            currmot.comprefinstances(regints,0);
            currmot.comprefmot();
            max=distcv(currmot.matprec,pmat);
            pmat=currmot.matprec;
            //        cout << max << endl;
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
      cout << currmot.nbmot << " " ; 
      cout.flush();

      currmot.comprefinstances(regints,0);
      currmot.matinit(scorethr2);
      if (currmot.nbmot>2){
         currmot.pvaluecomp();
         currmot.display(streamfile);
         //		for (ivma ima=currmot.seqs.begin();ima!=currmot.seqs.end();ima++){
         //  cout << (*ima).alignseq[0] << endl;
         //		}
         //cout << currmot.matprec << endl;
      }
      i++;
   }
   cout << endl;
}

vTSS TSSall;

   void
loadannots()
{
   ifstream annots;
   annots.open(args_info.phenotype_arg);

   back_insert_iterator<vstring> dest(phenos);
   copy(iisstring(annots),iisstring(),dest);

   annots.close();

   ifstream glist;
   //!!!! FIRST scangens were done with genelist-wosensory
   //   if (species==1) glist.open("/home/santolin/these/files/droso/plaza/phenos/neg-ovo-pheno.dat");
   //   else if (species==2) glist.open("/home/santolin/these/files/mus/affymetrix/e10/genelist-wo-pos-down-uniq.dat");
   if (args_info.phenoback_given){
      glist.open(args_info.phenoback_arg);
   }
   else {
      if (species==1) glist.open("/home/rouault/these/sequence/genomes/genelist.dat");
      else if (species==2) glist.open("/home/santolin/these/files/mus/biomart/genelist-protein-coding+miRNA.dat");
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

   vcoord
loadcoords()
{
   ifstream inf;
   inf.open(args_info.coord_file_arg);
   vcoord coords;
   back_insert_iterator<vcoord> dest(coords);
   copy(iiscoord(inf),iiscoord(),dest);
   return coords;
}

   void
loadchroms()
{
   ifstream chromsf;
   if (species==1){
      chromsf.open("/home/rouault/these/sequence/genomes/melano-only/dmel-all-chromosome-r4.3.fasta");
   } else if (species==2){
      chromsf.open("/home/santolin/these/files/mus/genome/fasta/all/all.chromosome_no_MT.fasta");
      //chromsf.open("/home/santolin/these/files/mus/genome/fasta/all/chr4.fa");
      //chromsf.open("/home/santolin/these/files/mus/genome/fasta/repeat_masked/all.chromosome_no_MT.fa");
   }
   back_insert_iterator<vchrom> dest(chromints);
   copy(iischrom(chromsf),iischrom(),dest);
   sort(chromints.begin(),chromints.end());

   for (int i=0;i<chromints.size();i++){
      cout << "chrom " << chromints[i].name << " ";
      cout << "of length " << chromints[i].seq.size() << " : ";
      cout << chromints[i].seq.substr(0,50) << " ... ";
      cout << chromints[i].seq.substr(chromints[i].seq.size()-50);
      cout << "\n";
   }


   chromsf.close();
}

   vint 
loadlengthchrom()
{

   vint lchr(nbchrom,0);

   if (species==1){

      lchr[0]=22410834;//chr2L
      lchr[1]=20769785;//chr2R
      lchr[2]=23774897;//chr3L
      lchr[3]=27908053;//chr3R
      lchr[4]=1284640;//chr4
      lchr[5]=22227390;//chrX


   } else if (species==2){

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
initgroupedinst()
{
   cout << "Loading length chroms..." << endl;
   lengthchrom=loadlengthchrom();

   vginst dumvginst;
   groupedinst=vvginst(nbchrom,dumvginst);
   for (unsigned int i=0;i<nbchrom;i++){
      for (unsigned int j=0;j<lengthchrom[i]/scanstep+1;j++){
         groupedinst[i].push_back(GroupInstance(j*scanstep,j*scanstep+scanwidth,i));
      }
   }
   ifstream annots;
   if (species==1){
      annots.open("/home/rouault/these/sequence/genomes/regres-wellform-all.dat");
   } else if (species==2){
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
getsites(ifstream & file, Sequence & bds)
{
   string dum;
   string dumtest;
   int i=0;
   cout << "seq " << endl;
   getline(file,dum);
   while(!file.eof()){
      bds.sitenames.push_back(dum.substr(1));
      getline(file,dum);
      if (i>0 && dum.length()!=dumtest.length()) {
         cout << dum << endl;
         cout << "Error in site length" << endl;
         exit(2);
      }
      dumtest=dum;
      i++;
      bds.seqs.push_back(dum);
      cout << dum << endl;
      bds.iseqs.push_back(stringtoint(dum));
      getline(file,dum);
   }
}

   void
disptex(vseq & seqs)
{
   system("if ! test -d display;then mkdir display;fi;");      
   //   system("if ! test -d display/tex;then mkdir display/tex;fi;");      

   double scoreinit=scorethr2;
   string filename("display/");
   filename+="results.tex";
   ofstream outf(filename.c_str());
   disptexinit(outf);
   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){ 
      cout << "Scanning " << (*ivs).name << endl;
      dispseqwmots(*ivs,motsdef,outf,scoreinit);
   }
   disptexclose(outf);

   outf.close();
}

   void
disptexwgaps(vseq & align)
{
   system("if ! test -d display;then mkdir display;fi;");      

   cout << "Scanning sequences for instances..." << endl;
   scanseqsforinstancesnmask(align,motsdef);

   cout << "Defining conserved instances..." << endl;
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
      ivs->instances2instancescons();
      //cout << ivs->name << "\n" << ivs->instancescons;
   }

   string folder("display/");
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){ 
      stringstream file;
      file << folder;
      file << ivs->name << "_";
      file << chromfromint(ivs->chrom) << "_";
      file << ivs->start << "_";
      file << ivs->stop;
      file << ".tex";
      ofstream outf(file.str().c_str());
      disptexinit(outf);
      cout << "Scanning " << (*ivs).name << endl;
      dispseqwmotswgaps(*ivs,outf);
      disptexclose(outf);
      outf.close();
   }

}


//Loads align-file (fasta) or coord-file (name/chrom/start/stop)
   vseq
loadseqs()
{
   vseq seqs;

   if (args_info.align_file_given){

      ifstream inf;
      inf.open(args_info.align_file_arg);
      seqs=loadsequencesconserv(inf);
      inf.close();
   }
   else if (args_info.coord_file_given){

      ifstream inf;
      inf.open(args_info.coord_file_arg);

      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align-masked.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align-masked.dat");

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

//Loads align-file (fasta) or coord-file (name/chrom/start/stop)
   vseq
loadseqscons()
{
   vseq seqs;

   if (args_info.align_file_given){

      ifstream inf;
      inf.open(args_info.align_file_arg);
      seqs=loadsequencesconservonly(inf);
      inf.close();
   }
   else if (args_info.coord_file_given){

      ifstream inf;
      inf.open(args_info.coord_file_arg);

      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align-masked.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align-masked.dat");

      alignscoord=loadcoordconserv(align);

      align.close();

      vcoord coords;
      back_insert_iterator<vcoord> dest(coords);
      copy(iiscoord(inf),iiscoord(),dest);
      for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
         Sequence seqtoimport=coordtoseq(*ivc);
         //         cout << (*ivc).name << " " << (*ivc).start << endl;
         //         cout << seqtoimport.name << " " << seqtoimport.start << endl;
         vint spe;
         int i=0;
         for (ivint iv=seqtoimport.species.begin();iv!=seqtoimport.species.end();iv++){
            if (*iv==1) spe.push_back(i);
            i++;
         }
         if (iscons(spe) && seqtoimport.species[0] && seqtoimport.nbtb>0){
            seqs.push_back(seqtoimport);
         }
         if (!iscons(spe)) cout << ivc->name << ": sequence not conserved\n";
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

//Loads align-file (fasta) or coord-file (name/chrom/start/stop)
   vseq
loadseqsnotmasked()
{
   vseq seqs;

   if (args_info.align_file_given){

      ifstream inf;
      inf.open(args_info.align_file_arg);
      seqs=loadsequencesconserv(inf);
      inf.close();
   }
   else if (args_info.coord_file_given){

      ifstream inf;
      inf.open(args_info.coord_file_arg);

      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align.dat");

      alignscoord=loadcoordconserv(align);

      align.close();

      vcoord coords;
      back_insert_iterator<vcoord> dest(coords);
      copy(iiscoord(inf),iiscoord(),dest);
      for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
         //         cout << (*ivc).name << " " << (*ivc).start << endl;
         Sequence seqtoimport=coordtoseq(*ivc);
         //         cout << seqtoimport.name << " " << seqtoimport.start << endl;
         //         for (int i=0;i<nbspecies;i++){
         //         cout << seqtoimport.iseqs[i] << endl; 
         //         }
         //         exit(9);
         vint spe;
         int i=0;
         for (ivint iv=seqtoimport.species.begin();iv!=seqtoimport.species.end();iv++){
            if (*iv==1) spe.push_back(i);
            i++;
         }
         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
            if (iscons(spe)){
               seqtoimport.cons=1;
            }
            else {
               //cout << ivc->name << ": sequence not conserved\n";
            }
            //      cout << "\r" <<  ivc->name;
            //      cout.flush();
            seqs.push_back(seqtoimport);
         }
      }
      inf.close();
   }
   else {
      cout << "Error in loadseqs: please give a coord/alignment file. Exiting..." << endl;
      exit(1);
   }
   return seqs;
}

//Loads align-file (fasta) or coord-file (name/chrom/start/stop) 
// keep only cons alignments!!
   vseq
loadseqsnotmaskedcons()
{
   vseq seqs;

   if (args_info.align_file_given){

      ifstream inf;
      inf.open(args_info.align_file_arg);
      seqs=loadsequencesconservonly(inf);
      inf.close();
   }
   else if (args_info.coord_file_given){

      ifstream inf;
      inf.open(args_info.coord_file_arg);

      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align.dat");

      alignscoord=loadcoordconserv(align);

      align.close();

      vcoord coords;
      back_insert_iterator<vcoord> dest(coords);
      copy(iiscoord(inf),iiscoord(),dest);
      for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
         Sequence seqtoimport=coordtoseq(*ivc);
         //         cout << (*ivc).name << " " << (*ivc).start << endl;
         //         cout << seqtoimport.name << " " << seqtoimport.start << endl;
         vint spe;
         int i=0;
         for (ivint iv=seqtoimport.species.begin();iv!=seqtoimport.species.end();iv++){
            if (*iv==1) spe.push_back(i);
            i++;
         }
         if (iscons(spe) && seqtoimport.species[0] && seqtoimport.nbtb>0){
            seqs.push_back(seqtoimport);
         }
         if (!iscons(spe)) cout << ivc->name << ": sequence not conserved\n";
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
weightmots()
{
   vseq seqs;
   cout << "Loading pos/neg sequences..." << endl;
   loadseqsforscore(seqs);

   for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      cout << ivm->name << endl;
      ivm->nbmot=0;
      width=ivm->motwidth;
      scorethr2=ivm->motscorethr2;
      ivm->calclambdaposneg(seqs);
      double scoresave=ivm->motscorethr2;
      while (ivm->lambda==0 || ivm->lambdatrain==0){
         ivm->motscorethr2=ivm->motscorethr2-0.5;
         ivm->calclambdaposneg(seqs);
         if (ivm->motscorethr2<0){
            ivm->motscorethr2=scoresave;
            ivm->calclambdaposneg(seqs);
            if (ivm->lambda==0 && ivm->lambdatrain>0){ivm->lambda=exp(1.)*ivm->lambdatrain;}
            else if (ivm->lambda>0 && ivm->lambdatrain==0){ivm->lambdatrain=exp(-1.)*ivm->lambda;}
            else { ivm->lambdatrain=1.e-4;ivm->lambda=ivm->lambdatrain/exp(1.);}
            break;
         }
      }
      cout << "Motif " << ivm->index+1 << " score " << ivm->motscorethr2 <<  " lb " << ivm->lambda << " lt " << ivm->lambdatrain << " log(lt/lb) " << log(ivm->lambdatrain/ivm->lambda) << endl;
      ivm->motscorethr2=scoresave;
   }
   return;
}

   void
args_init()
{
   width=args_info.width_arg;
   scorethr2=args_info.threshold_arg;
   species=args_info.species_arg;
   if (species==1){
      conca=0.3; 
      concc=0.2;
      nbspecies=12;
      cout << "Species: droso" << endl;
   }
   else if (species==2){
      conca=0.263; 
      concc=0.5-conca;
      nbspecies=10;
      cout << "Species: mus" << endl;
   }
   conct=conca;
   concg=concc;
   evolutionary_model=args_info.evolutionary_model_arg;

   scorethrcons=scorethr2-1.0;
   scorethr=scorethr2-1.0;

   neighbext=args_info.neighbext_arg;
   nbmots_for_score=args_info.nbmots_arg;
   scanwidth=args_info.scanwidth_arg;

   if (species==1) nbchrom=6;
   else if (species==2) nbchrom=21;

}

   void
args_init_scangen()
{
   scanwidth=args_info.scanwidth_arg;
   scanstep=args_info.scanstep_arg;

   nbmots_for_score=args_info.nbmots_arg;
}

   int
main(int argc, char** argv)
{
   if ( cmdline_parser(argc,argv, & args_info)!=0) 
      exit(1);

   args_init();

   rnginit();

   //                        compalpha();
   //                        exit(9);
   //   loadmots(args_info.motifs_arg,motsdef);
   //   inittreedist();

   if (args_info.scangen_given){

      compalpha();
      args_init_scangen();
      cout << "Loading Motifs" << endl;
      if (args_info.mots_w_names_given){
         loadmotswnames(args_info.mots_w_names_arg,motsdef);
      } 
      else{
         loadmots(args_info.motifs_arg,motsdef); 
      }
      cout << "Loaded " << motsdef.size() << " motifs." << endl;
      cout << "Nb mots for score: " << nbmots_for_score  << endl;

      if (args_info.opt_thr_file_given){
         ifstream opt;
         opt.open(args_info.opt_thr_file_arg);
         string dum;
         getline(opt,dum);
         int index=0;
         while (getline(opt,dum)){
            istringstream line(dum);
            line >> dum;
            double thr;
            line >> thr;
            motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
            motsdef[index].motscorethr=motsdef[index].motwidth*(thr-1)/10;
            motsdef[index].motscorethrcons=motsdef[index].motwidth*(thr-1)/10;
            index++;
         }
         opt.close();
      }
      else if (args_info.mean_info_given){
         cout << "Setting threshold to mean info..." << endl;
         for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
            ivm->setscorethr2meaninfo();
            cout << ivm->name << " mean thr = " << ivm->motscorethr2 << endl;
         }
      }
      else {
         cout << "Setting threshold to given value..." << endl;
         cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;
      }
      if (args_info.weightmots_given){
         cout << "Weighting Motifs" << endl;
         weightmots();
      }
      //      cout << "Setting threshold to mean info..." << endl;
      //      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      //         ivm->setscorethr2meaninfo();
      //         ivm->id.append("_");
      //         ivm->id.append(ivm->name);
      //         ivm->name=ivm->id;
      //      }

      if (args_info.phenotype_given){

         //cout << "Loading chroms" << "\n";
         //loadchroms();

         cout << "Loading phenotypes" << endl;
         loadannots();

         cout << "Loading instances" << endl;
         initgroupedinst();

         ifstream potregs;
         //if (species==1) potregs.open("/home/rouault/these/sequence/genomes/seqs.dat");
         //if (species==1) potregs.open("/home/santolin/these/files/droso/align/all/align-masked-files.dat");
         if (species==1) potregs.open("/home/santolin/these/files/droso/align/all/align-files.dat");
         else if (species==2) potregs.open("/home/santolin/these/files/mus/epo/align-files.dat");
         vstring regs;
         back_insert_iterator<vstring> dest(regs);
         copy(iisstring(potregs),iisstring(),dest);
         potregs.close();

         cout << "Scanning seqs" << endl;
         scanseqs(regs);

      } 
      else if (args_info.coord_file_given){

         cout << "Loading chromosomes" << "\n";
         loadchroms();

         cout << "Loading coordinates" << "\n";
         vcoord coords;
         coords=loadcoords();

         cout << "Loading instances" << endl;
         initgroupedinst();

         cout << "Loading aligments" << endl;
         ifstream align;
         if (species==1) align.open("/home/santolin/these/files/droso/align/all/align.dat");
         else if (species==2) align.open("/home/santolin/these/files/mus/epo/align.dat");
         alignscoord=loadcoordconserv(align);
         align.close();

         vstring regs;
         for (ivcoord ivc=alignscoord.begin();ivc!=alignscoord.end();ivc++){
            for (ivcoord ivc1=coords.begin();ivc1!=coords.end();ivc1++){
               if (ivc->chrom != ivc1->chrom || ivc->stop <= ivc1->start || ivc->start >= ivc1->stop) continue;
               else regs.push_back(ivc->name);
            }
         }

         for (ivstring iv=regs.begin();iv!=regs.end();iv++){
            cout << *iv << endl;
         }


         cout << "Scanning seqs" << endl;
         scanseqs(regs);
      }
   } 
   else if (args_info.display_given){


      cout << "Loading (unmasked) alignments " << endl;
      vseq align;
      align=loadseqsnotmasked();
      cout << "Nb sequences to scan: " << align.size() << endl;

      if (args_info.jaspardb_given){

         cout << "Loading Jaspar Motifs..." << endl;
         loadjaspardb(motsdef);
         cout << "Setting threshold to mean info..." << endl;
         for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
            ivm->setscorethr2meaninfo();
            ivm->id.append("_");
            ivm->id.append(ivm->name);
            ivm->name=ivm->id;
         }
         nbmots_for_score=motsdef.size();
         cout << "Loaded " << motsdef.size() << " motifs." << endl;

         cout << "Scanning sequences for instances..." << endl;
         scanseqsforinstances(align,motsdef);

         cout << "Defining conserved instances..." << endl;

         system("if ! test -d TFBS-pos;then mkdir TFBS-pos;fi;");      
         string folder="TFBS-pos/";
         ofstream outf;
         cout << "Writing TFBS-pos files..." << endl;
         for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
            ivs->instances2instancescons();
            if (ivs->instancescons.size()){
               stringstream file;
               file << folder;
               file << ivs->name << "_";
               file << chromfromint(ivs->chrom) << "_";
               file << ivs->start << "_";
               file << ivs->stop << ".dat";
               outf.open(file.str().c_str());
               outf << ivs->instancescons;
               outf.close();
            }
         }
         stringstream file;
         file << folder;
         file << "all-motifs.dat";
         outf.open(file.str().c_str());
         for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
            bool show(0);
            for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
               if (ivs->instancescons.size()){
                  for (ivvinstseq ivv=ivs->instancescons.begin();ivv!=ivs->instancescons.end();ivv++){
                     Instanceseq ist=ivv->at(0);
                     if (ist.motindex!=ivm->index) continue;
                     if (!show){
                        outf << ivm->name << "\t";
                        show=1;
                     }
                     outf << "chr" << ivs->chrom << ":";
                     int pos=ivs->start+ist.truepos;
                     if (ist.sens==1) outf << pos << "-" << pos+ivm->motwidth-1;
                     else if (ist.sens==-1) outf << pos-ivm->motwidth+1 << "-" << pos;
                     outf << " (" <<ivs->name << ",";
                     outf <<  ist.seq << ",";
                     outf << ist.score << ")";
                     //if (iv!=ivs->instancescons[0].end()-1) 
                     outf <<",\t";
                  }
               }
            }
            if (show==1) outf << endl;
         }
         outf.close();
      }
      else{

         cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;
         cout << "Loading Motifs" << endl;
         if (args_info.mots_w_names_given){
            loadmotswnames(args_info.mots_w_names_arg,motsdef);
         } else{
            loadmots(args_info.motifs_arg,motsdef); 
         }

         if (nbmots_for_score<motsdef.size()) motsdef.erase(motsdef.begin()+nbmots_for_score,motsdef.end());
         cout << "Loaded " << motsdef.size() << " motifs." << endl;
         for (ivmot iv=motsdef.begin();iv!=motsdef.end();iv++){
            iv->motscorethrcons=iv->motscorethr2;
         }
         if (args_info.opt_thr_file_given){
            ifstream opt;
            opt.open(args_info.opt_thr_file_arg);
            string dum;
            getline(opt,dum);
            int index=0;
            while (getline(opt,dum)){
               istringstream line(dum);
               line >> dum;
               double thr;
               line >> thr;
               motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
               motsdef[index].motscorethr=motsdef[index].motwidth*thr/10;
               motsdef[index].motscorethrcons=motsdef[index].motwidth*thr/10;
               index++;
            }
            opt.close();
         }


         if (args_info.disp_tex_given){

            ofstream outf;
            cout << "Creating fasta/tex files... " << endl;
            if (args_info.wgaps_given){
               disptexwgaps(align);
            } else{
               disptex(align);
            }
         }
         if (args_info.disp_svg_given){
            scanseqsforsvg(align);
         }
      }
   }
   else if (args_info.dispscore_given){
      cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;

      system("if ! test -d score;then mkdir score;fi;");      

      cout << "Loading Motifs..." << endl;
      if (args_info.mots_w_names_given) loadmotswnames(args_info.mots_w_names_arg,motsdef);
      else loadmots(args_info.motifs_arg,motsdef);
      cout << "Loaded " << motsdef.size() << " motifs." << endl;
      cout << "Nb mots for score: " << nbmots_for_score  << endl;

      if (args_info.opt_thr_file_given){
         ifstream opt;
         opt.open(args_info.opt_thr_file_arg);
         string dum;
         getline(opt,dum);
         int index=0;
         while (getline(opt,dum)){
            istringstream line(dum);
            line >> dum;
            double thr;
            line >> thr;
            motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
            motsdef[index].motscorethr=motsdef[index].motwidth*thr/10;
            motsdef[index].motscorethrcons=motsdef[index].motwidth*thr/10;
            index++;
         }
         opt.close();
      }
      if (args_info.weightmots_given){
         cout << "Weighting Motifs" << endl;
         weightmots();
      }

      vseq vscore;
      cout << "Loading pos/neg sequences..." << endl;
      loadseqsforscore(vscore);
      for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
         if (ivs->name.find('_')!=string::npos)
         {
            (*ivs).name.insert((*ivs).name.find('_'),"\\");
         }
      }

      system("if ! test -d score;then mkdir score;fi;");      
      system("if ! test -d score/dispscore;then mkdir score/dispscore;fi;");      
      calcscore(motsdef,vscore);
      sort(vscore.begin(),vscore.end());
      dispscore(vscore,motsdef);
      //for (ivmot iv=motsdef.begin();iv!=motsdef.end();iv++){
      //   dispscore(vscore,*iv);
      //}


   }
   else if (args_info.motgen_given){

      scorethr=scorethr2-2.*(double)width/10;//(double)width/10*6.5;
      scorethrcons=scorethr2-1.*(double)width/10;//(double)width/10*6.5;
      compalpha();
      cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;
      cout << "alpha: " << alpha << endl;
      ofstream motmeldb;
      motmeldb.open("motmeldb.txt");

      cout << "Loading background set..." << endl;
      ifstream backreg;
      //      if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000-align.dat"); //regs2000.fa");
      //      else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs2000-align.dat");//regs2000.fa");
      //      regtests=loadsequencesconserv(backreg);//loadsequences(backreg);
      if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000.fa");
      else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs2000-masked.fa");
      regtests=loadsequences(backreg);
      backreg.close();
      cout << "Background set size : " << regtests.size() << endl;

      cout << "Loading training set..." << endl;
      regints=loadseqs();
      cout << "Training set size : " << regints.size() << endl;
      inittreedist();

      for (vseq::iterator iseq=regints.begin();iseq!=regints.end();iseq++){
         cout << (*iseq).name << endl;
         if (args_info.ref_given){
            refseqanalysis(*iseq,motmeldb);
         }
         else{
            seqanalysis(*iseq,motmeldb);
         }
         cout << endl;
      }
      motmeldb.close();
   }
   else if (args_info.backreg_given){

      cout << "Loading genomic alignments..." << endl;

      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align.dat");//-masked.dat");
      alignscoord=loadcoordconserv(align);
      align.close();

      cout << "Loading intergenic coordinates..." << endl;
      vcoord vcds;

      ifstream intgen;
      if (species==1) intgen.open("/home/santolin/these/files/droso/genome/noncoding/nonCDS.dat");
      else if (species==2) intgen.open("/home/santolin/these/files/mus/all/files/noncoding.dat");
      vcds=loadcoordconserv(intgen);
      intgen.close();

      cout << "Random shuffling intergenic sequences..." << endl;

      random_shuffle(vcds.begin(),vcds.end());

      if (args_info.coord_file_given){
         //cout << "Loading available and conserved sequences from coordinate file..." << endl;
         //regints=loadseqsnotmaskedcons();
         cout << "Loading available sequences from coordinate file..." << endl;
         regints=loadseqsnotmasked();
         cout << "Loaded " << regints.size() << " sequences" << endl;
      }

      cout << "Writing background files..." << endl;
      system("if ! test -d backreg;then mkdir backreg;fi;");      
      string folder="backreg/";//"/home/santolin/these/files/mus/backreg/align/reg_";
      if (args_info.coord_file_given){
         printbackregwcoords(vcds,folder);

      } else{
         printbackreg(vcds,folder);
      }
   }
else {
   cout << "No mode given" << endl;
}


gsl_rng_free(gslran);
cmdline_parser_free(&args_info);
cout << "exit normally" << endl;
return 0;
}

