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

#include <cmath>
#include <cstdlib>
#include <algorithm>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multimin.h>
#include <math.h>
#include <numeric>
#include <time.h>
#include <cstring>
#include <errno.h>

#include "const.hpp"
#include "vectortypes.hpp"
#include "random.hpp"
#include "motif.hpp"
#include "sequence.hpp"
#include "genmot.hpp"
#include "tree.hpp"

using namespace std;

Motif::Motif()
{
   // *** only way to do if we want to define distwith in the cpp 
   distmot=new int[distwidth];
   for (unsigned int i=0;i<distwidth;i++){
      distmot[i]=0;
   }
   name="";
   nbmatch=0;
   nbmatchback=0;
   lambda=0;
   lambdatrain=0;
   pvalue=0;
   scorepoiss=0;
   nbmot=0;
   ntrain=0;
   check=true;
   vinst dumvinst;
   instances=vvinst(nbchrom,dumvinst);

   motscorethr2=scorethr2;
   motwidth=width;
   motscorethr=scorethr2-2*(double)motwidth/10;
   motscorethrcons=scorethr2-(double)motwidth/10;

   tottest=0;
}


   void
Motif::matinitforscanmots(Sequence & seq)
{
   unsigned int len=seq.iseqs[0].size();
   unsigned int i=0;
   for (vint::const_iterator istr=seq.iseqs[0].begin();istr!=seq.iseqs[0].end()-motwidth+1;istr++){
      if (scoref(istr,matprec)>motscorethr2){
         vint::const_iterator endci=seq.iseqs[0].end()-motwidth+1;
         unsigned int sh=shift(istr,matprec,endci,motwidth);
         Motalign ma(i+sh,seq,*this, 1);
         if (ma.iscons()){
            Instance refinst(seq.chrom,seq.start+i+sh,1,index,scoref(istr,matprec),vinttostring(ma.alignseq[0]));
            refinstances_short.push_back(refinst);
            if (i+sh<len-2*motwidth){
               istr+=motwidth-1+sh;
               i+=motwidth-1+sh;
            }
            else {
               break;
            }
         }
      }
      else if (scoref(istr,matprecrevcomp)>motscorethr2){ 
         vint::const_iterator endci=seq.iseqs[0].end()-motwidth+1;
         unsigned int sh=shift(istr,matprecrevcomp,endci,motwidth);
         Motalign ma(i+sh,seq,*this,-1);
         if (ma.iscons()){
            Instance refinst(seq.chrom,seq.start+i+sh,-1,index,scoref(istr,matprecrevcomp),vinttostring(ma.alignseq[0]));
            refinstances_short.push_back(refinst);
            if (i+sh<len-2*motwidth){
               istr+=motwidth-1+sh;
               i+=motwidth-1+sh;
            }
            else {
               break;
            }
         }
      }
      i++;
      if (istr>seq.iseqs[0].end()-motwidth) break;
   }
}
   

   void 
freqtolog(vvd & mat)
{
   int j=0;
   for (ivvd col=mat.begin();col!=mat.end();col++){
      int i=0;
      for (ivd line=(*col).begin();line!=(*col).end();line++){
         if (i==0||i==1){
            if ((*line)==0) (*line)=-1000;
            else (*line)=log((*line)/conca);
         }
         else {
            if ((*line)==0) (*line)=-1000;
            else (*line)=log((*line)/concc);	
         }
         //if ((*line)<-1000) (*line)=-1000;
         i++;
      }
      j++;
   } 
}

   void
countfreq(vvd & mat)
{
   int j=0;
   for (ivvd col=mat.begin();col!=mat.end();col++){

      int ntot(0);
      for (ivd line=(*col).begin();line!=(*col).end();line++){
         ntot+=(*line);
      }

      int i(0);
      for (ivd line=(*col).begin();line!=(*col).end();line++){
         if (i==0||i==1){
            (*line)=((*line)+alpha)/(ntot+2*alpha+2*beta);
         }
         else {
            (*line)=((*line)+beta)/(ntot+2*alpha+2*beta);
         }
         i++;
      }
      j++;
   } 
}

   void
countbases(Motif & mot,Sequence & bds)
{
   for (ivvint line=bds.iseqs.begin();line!=bds.iseqs.end();line++){
      int j=0;
      for (ivint i=(*line).begin();i!=(*line).end();i++){
         mot.matrice[j][*i]+=1;
         //cout << mot.matrice[j][*i] << " ";
         j++;
      }
      //cout << endl;
   }
}

   void
Motif::setscorethr2meaninfo()
{
   double maxinfo(0);
   double meaninfo(0);
   int j(0);
   for (ivvd ivv=matprec.begin();ivv!=matprec.end();ivv++){
      int i(0);
      double maxcol(-10);
      for (ivd iv=ivv->begin();iv!=ivv->end();iv++){
         if (*iv>maxcol) maxcol=*iv;
         meaninfo+=matfreq[j][i]*(*iv);
         i++;
      }
      j++;
      maxinfo+=maxcol;
   }
   motscorethr2=min(0.95*maxinfo,meaninfo);
   motscorethr=motwidth*(motscorethr2-1)/10;
   motscorethrcons=motwidth*(motscorethr2-1)/10;
}


GroupInstance::GroupInstance()
{
   nbmots=vint(nbmots_for_score,0);
   discarded=0;
   totmots=0;
}

GroupInstance::GroupInstance(int sta,int sto,int chr)
{
   start=sta;
   stop=sto;
   chrom=chr;
   nbmots=vint(nbmots_for_score,0);
   discarded=0;
   totmots=0;
}

   void
Motif::updatebacksites(Sequence & seq)
{
   unsigned int nbbacktemp=0;

   tottest+=seq.nbtb;
   if (motwidth>seq.nbtb){
      nbbacktemp=0;
   }
   else nbbacktemp=nbmatchcons(seq);
   if (nbbacktemp<distwidth) distmot[nbbacktemp]++;
   
   nbmatchback+=nbbacktemp;
}

   void
Motif::pvaluecomp()
{
   // Density of conserved binding sites in the background
   lambda=nbmatchback/(double)tottest;

   // Chi2 calculation
   calcscorepoiss();
   vseq::iterator iseq;
   pvalue=0.0;


   int nbtot=0;
   int nbbtrain=0;

   for (iseq=regints.begin();iseq!=regints.end();iseq++){
      // Here we use nbmot (cons) and not nmot (not cons)
      if (lambda>0) pvalue+=log(gsl_ran_poisson_pdf((*iseq).nbmot,lambda*(*iseq).nbtb));
      else pvalue+=0;
      nbbtrain+=(*iseq).nbtb;
      nbtot+=(*iseq).nbmot;
   }

   // Density of conserved binding sites in the training set
   lambdatrain=(double)nbtot/(double)nbbtrain;
   ntrain=nbtot;
}

   void
Motif::calcmeanpoiss()
{
   double sum=0;
   meanpoiss=0;
   for (unsigned int i=0;i<distwidth;i++){
      sum+=distmot[i];
      meanpoiss+=i*distmot[i];
   }
   meanpoiss/=sum;
}

   void
Motif::calcscorepoiss()
{
   calcmeanpoiss();
   scorepoiss=0.0;
   for (unsigned int j=0;j<distwidth;j++){
      if (distmot[j]>0 && meanpoiss>0){
         scorepoiss += gsl_pow_2((double)distmot[j]-regtests.size()*gsl_ran_poisson_pdf(j,meanpoiss))/(regtests.size()*gsl_ran_poisson_pdf(j,meanpoiss));
      }
   }
}

   void
Motif::corrprec()
{
   for (ivvd iv=matprec.begin();iv!=matprec.end();iv++){
      vd & col=*iv;
      double fra=conca*exp(col[0])+alpha;
      double frt=conct*exp(col[1])+alpha;
      double frc=concc*exp(col[2])+beta;
      double frg=concg*exp(col[3])+beta;
      double sum=fra+frt+frc+frg;
      fra/=sum;
      frt/=sum;
      frc/=sum;
      frg/=sum;
      col[0]=log(fra/conca);
      col[1]=log(frt/conct);
      col[2]=log(frc/concc);
      col[3]=log(frg/concg);
   }
}

   void
Motif::display(ostream & streamfile)
{
   streamfile << bsinit;
   streamfile << "\t";
   streamfile << pvalue << "\t";
   streamfile << scorepoiss << "\t";
   streamfile << nbmot << "\t";
   streamfile << ntrain << "\t";
   streamfile << setprecision(3);
   streamfile << lambdatrain << "\t";
   streamfile << lambda << "\t";
   streamfile << matprec;

   streamfile << setprecision(5);
   for (unsigned int i=0;i<distwidth;i++){
      streamfile << distmot[i];
      if (i!=distwidth-1) streamfile << ",";
   }
   streamfile << "\t";

   streamfile << endl;
}

Instance::Instance(int chr,int pos, int sen, int moti)
{
   motindex=moti;
   chrom=chr;
   coord=pos;
   sens=sen;
}

Instance::Instance(int chr,int pos, int sen, int moti,double sco,string s)
{
   motindex=moti;
   chrom=chr;
   coord=pos;
   sens=sen;
   score=sco;
   site=s;
}

ostream & operator <<(ostream &os,const Instance & inst)
{
   os << inst.motindex << "\t";
   os << inst.site << "\t";
   os << chromfromint(inst.chrom) << "\t";
   os << inst.coord << "\t";
   os << inst.sens << "\t";
   os << inst.score << endl;
   return os;
}

bool operator<(const Instance & inst1,const Instance & inst2)
{
   if (inst1.chrom != inst2.chrom) return inst1.chrom < inst2.chrom;
   else return inst1.coord < inst2.coord;
}

bool operator<(const GroupInstance & ginst1,const GroupInstance & ginst2)
{
   return ginst1.score > ginst2.score;
}

ostream & operator <<(ostream &os,const GroupInstance & ginst)
{
   os << chromfromint(ginst.chrom) << " ";
   os << ginst.start << " ";
   os << ginst.stop << " ";
   os << ginst.score << " ";
   os << "\n";
   for (civinst iv=ginst.instances.begin();iv!=ginst.instances.end();iv++){
      os << *iv;
   }
   return os;
}

bool operator<(const Combination & comb1,const Combination & comb2)
{
   //if score are equal, return highest number of positives first
   if (comb1.score==comb2.score) return comb1.TP > comb2.TP;
   else return comb1.score > comb2.score;
}

   ostream &
operator <<(ostream &os,const Combination & comb)
{
   for (civint ivi=comb.motis.begin();ivi!=comb.motis.end();ivi++){
      os << *ivi << "\t";
   }
   os << comb.TPR << "\t";
   os << comb.TP << "\t";
   os << comb.FPR << "\t";
   os << comb.FP << "\t";
   os << comb.score;

   return os;
}

   ostream &
operator <<(ostream &os,const vcombi & vcomb)
{
   for (civcombi ivc=vcomb.begin();ivc!=vcomb.end();ivc++){
      for (civint ivi=ivc->motis.begin();ivi!=ivc->motis.end();ivi++){
         os << *ivi << "\t";
      }
      os << ivc->TPR << "\t";
      os << ivc->TP << "\t";
      os << ivc->FPR << "\t";
      os << ivc->FP << "\t";
      os << ivc->score << "\n";
   }

   return os;
}

   void
GroupInstance::compbestannot()
{
   unsigned int dist=annotextent;
   for (ivTSS ivt=TSSs.begin();ivt!=TSSs.end();ivt++){
      if (abs((int)((*ivt).coord)-(int)start+(int)scanwidth/2)<(int)dist){
         dist=abs((int)(*ivt).coord-(int)start+(int)scanwidth/2);
         besttss=*ivt;
         //cout << besttss.gene << endl;
      }
   }
}

   int
Motif::nbmatchcons (Sequence & seq)
{
   width=motwidth;
   int nmat=0;
   unsigned int i=0;
   unsigned int len=seq.iseqs[0].size();
   for (civint istr=seq.iseqs[0].begin();istr!=seq.iseqs[0].end()-width+1;istr++){
      double score=scoref(istr,matprec);

      if (score>motscorethr2){
         Motalign ma(i,seq,*this, 1);
         if (ma.iscons()){
            nmat++;
            if (i<len-2*width){
               istr+=width-1;
               i+=width-1;
            }
            else {
               break;
            }
         }
      }
      else {
         score=scoref(istr,matprecrevcomp);
         if (score>motscorethr2){
            Motalign ma(i,seq,*this, -1);
            if (ma.iscons()){
               nmat++;
               if (i<len-2*width){
                  istr+=width-1;
                  i+=width-1;
               }
               else {
                  break;
               }
            }
         }
      }
      i++;
   }
   return nmat;
}


   void
Motif::findinstancesnmask (Sequence & seq)
{ 
   double thr;
   for (unsigned int spe=0;spe<nbspecies;spe++){
      if (spe==0) thr=motscorethr2;
      else thr=motscorethrcons;
      if (seq.species[spe]){
         unsigned int len=seq.iseqs[spe].size();
         unsigned int i=0;
         for (civint istr=seq.iseqs[spe].begin();istr!=seq.iseqs[spe].end()-motwidth+1;istr++){
            double score=scoref(istr,matprec);
            if (score>thr){
               //pos relative to start and including gaps (for inter-species comparison)
               Instanceseq inst(index,1,seq.imapsinv[spe][i],i,score,spe,name);
               inst.iseq=vint(istr,istr+motwidth);
               inst.seq=vinttostring(inst.iseq);
               seq.instances.push_back(inst);

               for (unsigned int j=i;j<i+motwidth;j++){
                  seq.iseqs[spe][j]=4;
               }
               if (i<len-2*motwidth){
                  istr+=motwidth-1;
                  i+=motwidth-1;
               }
               else {
                  break;
               }
            }
            else {
               score=scoref(istr,matprecrevcomp);
               if (score>thr){
                  Instanceseq inst(index,-1,seq.imapsinv[spe][i],i,score,spe,name);
                  inst.iseq=vint(istr,istr+motwidth);
                  inst.seq=vinttostring(inst.iseq);
                  seq.instances.push_back(inst);

                  for (unsigned int j=i;j<i+motwidth;j++){
                     seq.iseqs[spe][j]=4;
                  }
                  if (i<len-2*motwidth){
                     istr+=motwidth-1;
                     i+=motwidth-1;
                  }
                  else {
                     break;
                  }
               }
            }
            i++;
         }
      }
   }
   return;
}

   void
Motif::findinstances (Sequence & seq)
{ 
   double thr;
   for (unsigned int spe=0;spe<nbspecies;spe++){
      if (spe==0) thr=motscorethr2;
      else thr=motscorethrcons;
      if (seq.species[spe]){
         unsigned int len=seq.iseqs[spe].size();
         unsigned int i=0;
         for (civint istr=seq.iseqs[spe].begin();istr!=seq.iseqs[spe].end()-motwidth+1;istr++){
            double score=scoref(istr,matprec);
            if (score>thr){
               //pos relative to start and including gaps (for inter-species comparison)
               Instanceseq inst(index,1,seq.imapsinv[spe][i],i,score,spe,name);
               inst.iseq=vint(istr,istr+motwidth);
               inst.seq=vinttostring(inst.iseq);
               seq.instances.push_back(inst);

               if (i<len-2*motwidth){
                  istr+=motwidth-1;
                  i+=motwidth-1;
               }
               else {
                  break;
               }
            }
            else {
               score=scoref(istr,matprecrevcomp);
               if (score>thr){
                  Instanceseq inst(index,-1,seq.imapsinv[spe][i],i,score,spe,name);
                  inst.iseq=vint(istr,istr+motwidth);
                  inst.seq=vinttostring(inst.iseq);
                  seq.instances.push_back(inst);

                  if (i<len-2*motwidth){
                     istr+=motwidth-1;
                     i+=motwidth-1;
                  }
                  else {
                     break;
                  }
               }
            }
            i++;
         }
      }
   }
   return;
}

   void
Motif::findinstances (vseq & vs)
{ 
   for (ivseq ivs=vs.begin();ivs!=vs.end();ivs++){
      ivs->instances.clear();
      findinstances(*ivs);
   }
   return;
}
   int
Motif::dispmots (Sequence & seq, int motindex)
{
   int nmat=0;
   unsigned int i=0;
   unsigned int len=seq.iseqs[0].size();
   for (civint istr=seq.iseqs[0].begin();istr!=seq.iseqs[0].end()-width+1;istr++){
      double score=scoref(istr,matprec);
      if (score>motscorethr2){
         Motalign ma(i,seq,*this, 1);
         nmat++;
         cout << "motif " << motindex << " at position : " << i;
         if (ma.iscons()) cout << " (conserved)\\\\\n";
         else cout << "\\\\\n";
         if (i<len-2*width){
            istr+=width-1;
            i+=width-1;
         }
         else {
            break;
         }
      }
      else {
         score=scoref(istr,matprecrevcomp);
         if (score>motscorethr2){
            Motalign ma(i,seq,*this, -1);
            nmat++;
            cout << "motif " << motindex << " (revcom) at position : " << i;
            if (ma.iscons()) cout << " (conserved)\\\\\n";
            else cout << "\\\\\n";
            if (i<len-2*width){
               istr+=width-1;
               i+=width-1;
            }
            else {
               break;
            }
         }
      }
      i++;
   }
   return nmat;
}

   int
Motif::statemot (Sequence & seq,int pos, int num, double & scoremot)
{
   double scthr(motscorethr2);
   if (num>0) scthr=motscorethrcons;
   civint istr=seq.iseqs[num].begin()+pos;
   double sc=scoref(istr,matprec);
   vint temp;
   for (civint iv=istr;iv!=istr+width;iv++){
      temp.push_back(*iv);
   }
   if (sc>scthr){
      scoremot=sc;
      if (num==0){
         Motalign ma(pos,seq,*this, 1);
         if (ma.iscons()) return 2;
      }
      return 1;
   }
   else { 
      sc=scoref(istr,matprecrevcomp);
      if (sc>scthr){
         scoremot=sc;
         if (num==0){
            Motalign ma(pos,seq,*this, -1);
            if (ma.iscons()) return 2;
         }
         return 1;
      }
   }
   return 0;
}
   
   void
Motif::matinit(double scth)
{
   seqs.clear();
   for (ivseq iseq=regints.begin();iseq!=regints.end();iseq++){
      (*iseq).nmot=0;
      (*iseq).nbmot=0;
   }
   //   cout << "begin init\n";
   for (ivseq iseq=regints.begin();iseq!=regints.end();iseq++){
      Sequence & seq=*iseq;
      unsigned int len=seq.iseqs[0].size();
      unsigned int i=0;
      for (vint::const_iterator istr=seq.iseqs[0].begin();istr!=seq.iseqs[0].end()-width+1;istr++){
         //         double score1,score2;
         //         if (scoref(istr,matprec)>scth){
         //            vint::const_iterator endci=seq.iseqs[0].end()-width+1;
         //            unsigned int sh=shift(istr,matprec,endci,width);
         //            score1=scoref(istr+sh,matprec);
         //         }
         //         else if (scoref(istr,matprecrevcomp)>scth){ 
         //            vint::const_iterator endci=seq.iseqs[0].end()-width+1;
         //            unsigned int sh=shift(istr,matprecrevcomp,endci,width);
         //            score2=scoref(istr+sh,matprecrevcomp);
         //         }

         if (scoref(istr,matprec)>scth){
            vint::const_iterator endci=seq.iseqs[0].end()-width+1;
            unsigned int sh=shift(istr,matprec,endci,width);
            Motalign ma(i+sh,seq,*this, 1);
            //            ivint idro=ma.matches.begin();
            //            for (ivvint iv=ma.alignseq.begin();iv!=ma.alignseq.end();iv++){
            //               cout << *idro << " " << *iv << endl;
            //               idro++;
            //            }
            //            cout << endl;
            (*iseq).nmot++;
            if (ma.iscons()){
               //               ma.print();
               iseq->nbmot++;
               seqs.push_back(ma);
               if (i+sh<len-2*width){
                  istr+=width-1+sh;
                  i+=width-1+sh;
               }
               else {
                  break;
               }
            }
         }
         else if (scoref(istr,matprecrevcomp)>scth){ 
            vint::const_iterator endci=seq.iseqs[0].end()-width+1;
            unsigned int sh=shift(istr,matprecrevcomp,endci,width);
            Motalign ma(i+sh,seq,*this,-1);
            //            for (ivvint iv=ma.alignseq.begin();iv!=ma.alignseq.end();iv++){
            //               cout << *iv << endl;
            //            }
            //            cout << endl;
            (*iseq).nmot++;
            if (ma.iscons()){
               iseq->nbmot++;
               seqs.push_back(ma);
               //               ma.print();
               if (i+sh<len-2*width){
                  istr+=width-1+sh;
                  i+=width-1+sh;
               }
               else {
                  break;
               }
            }
         }
         i++;
         if (istr>seq.iseqs[0].end()-width) break;
      }
   }
   //   cout << "end init\n";
   nbmot=seqs.size();
}

   bool
Motalign::iscons()
{
   if (species=="droso"){
      int nbfr=0;
      if (matches[5]) nbfr++; 
      if (matches[6] || matches[7]) nbfr++;
      if (matches[8]) nbfr++;
      if (matches[9] || matches[10] || matches[11]) nbfr++;
      if (nbfr>1) return true;
   }
   else if (species=="eutherian"){
      int nbfr=0;
      if (matches[2] || matches[3] || matches[4] || matches[5] || matches[6] || matches[7]) nbfr++; // primates
      if (matches[8]) nbfr++; // horse
      if (matches[9]) nbfr++; // dog
      if (matches[10]) nbfr++; // wild boar
      if (matches[11]) nbfr++; // cow
      if (nbfr>1) return true;
   }
   return false;	
}

   void
Motalign::mask()
{
   for (unsigned int i=1;i<nbspecies;i++){
      if (matches[i]){
         for (unsigned int j=0;j<width;j++){
            *(matchespos[i]+j)=4;
         }
      }
   }
}

   void
Motif::compprec()
{
   vd col;
   vvd matopti;

   for (unsigned int i=0;i<width;i++){
      col=colopti(i,this);
      matopti.push_back(col);
   }
   matprec=matopti;
   matprecrevcomp=reversecomp(matprec);
}

   void
Motif::compprec_MCMC()
{
   vd col;
   vvd matmean;

   for (unsigned int i=0;i<motwidth;i++){
      col=colmean(i,this);
      matmean.push_back(col);
   }
   matprec=matmean;
   matprecrevcomp=reversecomp(matprec);
}

   vvd
mattofreq(vvd & mat)
{
   vd dumd(4,0.0);
   vvd mfreq(width,dumd);
   int i=0;
   for (ivvd iv=mat.begin();iv!=mat.end();iv++){
      vd & col=*iv;
      mfreq[i][0]=conca*exp(col[0]);//A
      mfreq[i][1]=conct*exp(col[1]);//T
      mfreq[i][2]=concc*exp(col[2]);//C
      mfreq[i][3]=concg*exp(col[3]);//G
      // AVOID PROBLEMS OF >1 TOTAL FREQUENCY!
      double sum=mfreq[i][0]+mfreq[i][1]+mfreq[i][2]+mfreq[i][3];
      mfreq[i][0]/=sum;
      mfreq[i][1]/=sum;
      mfreq[i][2]/=sum;
      mfreq[i][3]/=sum;
      i++;
   }
   return mfreq;
}

   vd
colopti(unsigned int pos,Motif * mot)
{
   const gsl_multimin_fminimizer_type *T = 
      gsl_multimin_fminimizer_nmsimplex;
   gsl_multimin_fminimizer *s = NULL;
   gsl_vector *ss, *x;
   gsl_multimin_function minex_func;

   /* Starting point */
   x = gsl_vector_alloc (3);
   gsl_vector_set (x, 0, .25);
   gsl_vector_set (x, 1, .25);
   gsl_vector_set (x, 2, .25);

   /* Set initial step sizes to 0.05 */
   ss = gsl_vector_alloc (3);
   gsl_vector_set_all (ss, 0.1);

   /* Initialize method and iterate */
   minex_func.n = 3;
   minex_func.f = &loglikely;
   void * par[2]={(void *)mot,&pos};
   minex_func.params = par;

   s = gsl_multimin_fminimizer_alloc (T, 3);
   gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

   int iter=0;
   int status=0;
   double size;
   do
   {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status) 
         break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-4);

      if (status == GSL_SUCCESS)
      {
         //printf ("converged to minimum at\n");
      }

      //            printf ("%5d %10.4e %10.4e %10.4e f() = %7.4f size = %.5f\n", 
      //                  iter,
      //                  gsl_vector_get (s->x, 0), 
      //                  gsl_vector_get (s->x, 1), 
      //                  gsl_vector_get (s->x, 2), 
      //                  s->fval, size);
   }
   while (status == GSL_CONTINUE && iter < 1000);
   vd res;
   double w0=gsl_vector_get(s->x,0);
   double w1=gsl_vector_get(s->x,1);
   double w2=gsl_vector_get(s->x,2);
   double w3=1.0-w0-w1-w2;
   res.push_back(log(w0/conca));
   res.push_back(log(w1/conct));
   res.push_back(log(w2/concc));
   res.push_back(log(w3/concg));

   gsl_vector_free(x);
   gsl_vector_free(ss);
   gsl_multimin_fminimizer_free (s);
   return res;
}

   vd
colmean(unsigned int pos,Motif * mot)
{

   double al(alpha),be;
   if (al<0.05) al=0.1;
   be=concc/conca*al;
   //   clock_t start,finish;
   //   double dif;
   //   start=clock();
   unsigned int Ndir(0);
   //cout << mot->seqs.size() << endl;
   for (ivma ivm=mot->seqs.begin();ivm!=mot->seqs.end();ivm++){
      Ndir+=accumulate(ivm->matches.begin(),ivm->matches.end(),0)/2; 
   }
   //cout << "Ndir=" << Ndir << endl;

   //unsigned int nsample(20);
   unsigned int nsample(10);
   //cout << "nsample=" << nsample << endl;

   //unsigned int nmean(100);
   //cout << "nmean=" << nmean << endl;
   double dkl=1.;
   double dcutoff(1e-4);
   //cout << "dcutoff=" << dcutoff << endl;

   // INIT
   vd winit(4,0.);
   winit[0]=conca*exp(mot->matprec[pos][0]);
   winit[1]=conct*exp(mot->matprec[pos][1]);
   winit[2]=concc*exp(mot->matprec[pos][2]);
   winit[3]=concg*exp(mot->matprec[pos][3]);
   double sum=winit[0]+winit[1]+winit[2]+winit[3];
   winit[0]/=sum;
   winit[1]/=sum;
   winit[2]/=sum;
   winit[3]/=sum;

   double fprev,fafter;
   void * par[2]={(void *)mot,&pos};

   vd wprev(4,0.),wafter(4,0.),wmean(4,0.),wmeanforstat(4,0.);

   // Starting point
   double alpinit[4]={al+Ndir*winit[0],al+Ndir*winit[1],be+Ndir*winit[2],be+Ndir*winit[3]};
   double thetainit[]={0.,0.,0.,0.};
   gsl_ran_dirichlet(gslran,4,alpinit,thetainit);
   wprev[0]=thetainit[0];
   wprev[1]=thetainit[1];
   wprev[2]=thetainit[2];
   wprev[3]=1-thetainit[0]-thetainit[1]-thetainit[2];

   fprev=loglikelyhood(wprev,par);

   // =========
   // MAIN LOOP
   // =========

   // counter counts accepted trials 
   unsigned int counter(1),iterok(1),iter(0);
   double logalph,prob; // for Metropolis Hastings
   // BURN-IN
   unsigned int bicount(0),bicutoff(nsample);

   //while (counter<=nmean) {
   while (dkl>=dcutoff) {

      double alp[4]={al+Ndir*wprev[0],al+Ndir*wprev[1],be+Ndir*wprev[2],be+Ndir*wprev[3]};
      double theta[]={0.,0.,0.,0.};
      gsl_ran_dirichlet(gslran,4,alp,theta);
      wafter[0]=theta[0];
      wafter[1]=theta[1];
      wafter[2]=theta[2];
      wafter[3]=1-theta[0]-theta[1]-theta[2];
         
      // check if new probabilities are consistent, else reject
      int wtest(0);
      for (int i=0;i<4;i++) if (wafter[i]<0 || wafter[i]>1) wtest=1;
      if (wtest) continue;

      fafter=loglikelyhood(wafter,par);

      // Metropolis-Hastings
      // Move to new position with probability alpha, else reject
      double theta2[]={wprev[0],wprev[1],wprev[2],wprev[3]};
      double alp2[4]={al+Ndir*wafter[0],al+Ndir*wafter[1],be+Ndir*wafter[2],be+Ndir*wafter[3]};
      logalph=min(0.,fafter+log(gsl_ran_dirichlet_pdf(4,alp2,theta2))-(fprev+log(gsl_ran_dirichlet_pdf(4,alp,theta))));

      prob=gsl_rng_uniform (gslran);

      if (log(prob)<=logalph){

         wprev=wafter;
         fprev=fafter;

         if (bicount<bicutoff){
            bicount++;
            continue;
         }

         if (iterok%nsample==0){
            vd wprevmean=wmeanforstat;
            for (unsigned int ib=0;ib<4;ib++){
               //wmeanforstat[ib]+=wafter[ib]/nmean;
               wmeanforstat[ib]=(counter-1.)/counter*wprevmean[ib]+wafter[ib]/counter;
            }

            if (wmeanforstat[0]<=0 || wmeanforstat[1]<=0 || wmeanforstat[2]<=0 || wmeanforstat[3]<=0){
               for (unsigned int i=0;i<4;i++){
                  if (wmeanforstat[i]<0) cout << "Warning in colmean, pos=" << pos << ", base=" <<
                     inttostring(i) << ", wmeanforstat=" << wmeanforstat[i] << endl;
               }
               wmeanforstat=wprevmean;
               iter++;
               continue;
            }

            dkl=0.;
            for (int ib=0;ib<4;ib++){
               dkl+=wmeanforstat[ib]*log(wmeanforstat[ib]/wprevmean[ib])/log(2);
            }

            counter++;
         }
               
         
         for (unsigned int ib=0;ib<4;ib++) wmean[ib]=(iterok-1.)/iterok*wmean[ib]+wafter[ib]/iterok;
        
         iterok++;
      }
      iter++;
   } 

   //      finish=clock();
   //      dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
   //   printf ("---> time(colmean): %.1f milliseconds.\n", dif );
   //   cout << "iter=" << iter << " ";
   //   cout << "iterOK=" << iterok << " ";
   //   cout << "counter=" << counter << endl;

   // HERE WE AVOID A PRECISION PROBLEM FOR w=0
   // VALUES ARE FLUCTUATING AND CAN BE NEGATIVE
   for (unsigned int i=0;i<4;i++){
      if (wmean[i]<0) cout << "Warning in colmean, pos=" << pos << ", base=" <<
         inttostring(i) << ", wmean=" << wmean[i] << endl;
   }
   sum=0;
   for (unsigned int i=0;i<4;i++){
      wmean[i]=fabs(wmean[i]);
      sum+=wmean[i];
   }
   for (unsigned int i=0;i<4;i++) wmean[i]/=sum;

   vd res;
   res.push_back(log(wmean[0]/conca));
   res.push_back(log(wmean[1]/conct));
   res.push_back(log(wmean[2]/concc));
   res.push_back(log(wmean[3]/concg));


   return res;
}

   void
Motalign::print()
{
   int i=0;
   for (ivint imat=matches.begin();imat!=matches.end();imat++){
      if (*imat){
         cout << vinttostring(alignseq[i]) << "\n";
      }
      i++;
   }
   cout << "\n";
}

Motalign::Motalign()
{
}

Motalign::Motalign(unsigned int pos, Sequence & seq, Motif & mot,int sens)
{
   unsigned int motwidth=mot.motwidth;
   seq_start=seq.iseqs[0].begin();
   seq_stop=seq.iseqs[0].end();
   civvint imap=seq.imaps.begin()+1;
   int seqnum=1;
   matches=vint(nbspecies,0);
   vint seqdum=vint(motwidth,4);
   ivint posdum;
   matchespos.push_back(seq.iseqs[0].begin()+pos);
   vint seqmel=vint(seq.iseqs[0].begin()+pos,seq.iseqs[0].begin()+pos+motwidth);
   //   cout << seqmel << endl;
   if (sens==1){
      alignseq.push_back(seqmel);
      strand=1;
      //                  cout << vinttostring(seqmel) << "\n";
   } else {
      alignseq.push_back(reversecomp(seqmel));
      strand=-1;
      //      vint v=reversecomp(seqmel);
      //                  cout << vinttostring(v) << "\n";
   }
   matches[0]=1;
   for (ivvint is=seq.iseqs.begin()+1;is!=seq.iseqs.end();is++){
      if (seq.species[seqnum]){
         int truepos=(*imap)[seq.imapsinv[0][pos]];
         int start=(int)truepos-neighbext;
         if (start<0) start=0;
         unsigned int stop=truepos+neighbext;
         if (stop>(*is).size()-motwidth) stop=(*is).size()-motwidth;
         int hasmatch=0;
         for (unsigned int i=start;i<stop;i++){
            civint startsite=(*is).begin()+i;
            if (sens==1){
               if (scoref(startsite,mot.matprec)>mot.motscorethrcons){
                  vint::const_iterator endci=(*is).end()-motwidth+1;
                  unsigned int sh=shift((*is).begin()+start,mot.matprec,endci,2*neighbext);
                  //                  int deca=(int)sh-neighbext;
                  //                  if (seqnum==1) cout << deca << "\n";
                  vint seqt=vint((*is).begin()+start+sh,(*is).begin()+start+sh+motwidth);
                  //                  cout << vinttostring(seqt) << "\n";
                  alignseq.push_back(seqt);
                  matchespos.push_back((*is).begin()+start+sh);
                  hasmatch=1;
                  //                  cout << i-truepos << "\n";
                  break;
               }
            } else if (sens==-1){
               if (scoref(startsite,mot.matprecrevcomp)>mot.motscorethrcons){
                  vint::const_iterator endci=(*is).end()-motwidth+1;
                  unsigned int sh=shift((*is).begin()+start,mot.matprecrevcomp,endci,2*neighbext);
                  //                  int deca=-((int)sh-neighbext);
                  //                  if (seqnum==1) cout << deca << "\n";
                  vint mot=vint((*is).begin()+start+sh,(*is).begin()+start+sh+motwidth);
                  vint v=reversecomp(mot);
                  //                  cout << vinttostring(v) << "\n";
                  alignseq.push_back(reversecomp(mot));
                  matchespos.push_back((*is).begin()+start+sh);
                  hasmatch=1;
                  //                  cout << i-truepos << "\n";
                  break;
               }
            }
         }
         if (hasmatch){
            matches[seqnum]=1;
         } else {
            alignseq.push_back(seqdum);
            matchespos.push_back(posdum);
            matches[seqnum]=0;
         }
      }
      else {
         alignseq.push_back(seqdum);
         matchespos.push_back(posdum);
         matches[seqnum]=0;
      }
      seqnum++;
      imap++;
   }
   //   cout << "\n";
}


   void
loadmots ( const char * filename, vmot & mots )
{
   if (!filename){
      cout << "Please give a motifs file. Exiting..." << endl;
      exit(1);
   }

   ifstream fmotifs;
   fmotifs.open(filename);
   if (fmotifs.fail()){
      cerr << "Cannot open motif file: " << strerror(errno) << endl;
      exit(-1);
   }

   string dum;
   fmotifs >> dum;

   unsigned int i(0);//!!to comply withh cpp convention
   while (! fmotifs.eof()){
      Motif mot1;
      mot1.index=i;
      stringstream name;
      name << "Mot_";
      name << i+1;
      name >> mot1.name;
      mot1.bsinit=dum;
      width=mot1.bsinit.size();
      mot1.motwidth=mot1.bsinit.size();
      fmotifs >> mot1.pvalue;
      fmotifs >> mot1.scorepoiss;
      fmotifs >> mot1.nbmot;
      fmotifs >> mot1.ntrain;
      fmotifs >> mot1.lambdatrain;
      fmotifs >> mot1.lambda;
      vd dumd(4,0.0);
      vvd dummat(mot1.motwidth,dumd);
      mot1.matprec=dummat;
      fmotifs >> mot1.matprec;
      for (ivvd ivmat=mot1.matprec.begin();ivmat!=mot1.matprec.end();ivmat++){
         for (ivd ivrow=ivmat->begin();ivrow!=ivmat->end();ivrow++){
            if (*ivrow<-6.5) *ivrow=-6.5;
         }
      }
      mot1.matfreq=mattofreq(mot1.matprec);
      mot1.matprecrevcomp=reversecomp(mot1.matprec);
      fmotifs >> mot1.distmot;
      fmotifs >> dum;
      mot1.motscorethr2=mot1.motwidth*scorethr2/10;
      mot1.motscorethr=mot1.motwidth*(scorethr2-1)/10;
      mot1.motscorethrcons=mot1.motwidth*(scorethr2-1)/10;
      mots.push_back(mot1);
      i++;
   }
   fmotifs.close();

   if (i==0){//No motifs
      cout << "No motifs" << endl;
      exit(1);
   } else if (i<nbmots_for_score-1) {
      nbmots_for_score=i;
      //cout << "Changed nbmots_for_score to value " << i << endl;
   }
}		/* -----  end of function loadmots  ----- */


   void
GroupInstance::compscore(vmot & lmots,unsigned int nbmots_score)
{
   score=0;
   unsigned int imot=0;
   for (ivmot ivm=lmots.begin();ivm!=lmots.begin()+nbmots_score;ivm++){
      //      ivmot ivm=lmots.begin()+nbmots_score-1;
      score+=nbmots[imot]*log((*ivm).lambdatrain/(*ivm).lambda);
      imot++;
   }
}


   int
GroupInstance::distance(const GroupInstance & gi)
{
   if (chrom!=gi.chrom) return 10000000;
   else return abs(start-gi.start);
}

   void
GroupInstance::isdiscarded()
{
   if (besttss.gene=="") discarded=1;
//   for (ivTSS ivt=TSSs.begin();ivt!=TSSs.end();ivt++){
//      if ((*ivt).gene=="phyl" ||
//            (*ivt).gene=="spdo" ||
//            (*ivt).gene=="vvl" ||
//            (*ivt).gene=="neur" ||
//            (*ivt).gene=="chn" ||
//            (*ivt).gene=="PFE" ||
//            (*ivt).gene=="CG32150" ||
//            (*ivt).gene=="mira" ||
//            (*ivt).gene=="CG9363" ||
//            (*ivt).gene=="cpo" ||
//            (*ivt).gene=="sens" ||
//            (*ivt).gene=="CG32392" ||
//            (*ivt).gene=="sv" ||
//            (*ivt).gene=="insv") {
//         discarded=1;
//         break;
//      }
//   }
} 

   void
displayhist(vginst & vgi,ostream & ostr)
{
   ostr << "Pos\tName\n";
   unsigned int pos=1;
   for (ivginst ivg=vgi.begin();ivg!=vgi.end();ivg++){
      if (!(*ivg).discarded){
         if ((*ivg).goodpheno){
            ostr << pos << "\t" << (*ivg).besttss.gene << "\n";
         }
         //if ((*ivg).besttss.gene.substr(0,2)!="CG"){
            pos++;
         //}
      }
   }
}

   void
displayhist_set(vginst & vgi, vstring geneset,ostream & ostr)
{
   ostr << "Score\tName\n";
   for (ivstring ivs=geneset.begin();ivs!=geneset.end();ivs++){
      double scoregene=0;
      for (ivginst ivg=vgi.begin();ivg!=vgi.end();ivg++){
         if ((*ivg).besttss.gene==*ivs){
            if ((*ivg).score>scoregene){
               scoregene=(*ivg).score;
            }
         }
      }
      ostr << scoregene << "\t" << *ivs << "\n";
   }
}

double ic;

   double
funcroot (double x, void *params)
{

   double f=2*conca*(gsl_sf_psi(x+1)-gsl_sf_psi(x/conca+1)-log(conca))+2*concc*(gsl_sf_psi(concc*x/conca+1)-gsl_sf_psi(x/conca+1)-log(concc))-ic;

   return f;
}

   double
funcroot_deriv (double x, void *params)
{

   double df=2*conca*(gsl_sf_psi_1(x+1)-gsl_sf_psi_1(x/conca+1)/conca)+2*concc*(gsl_sf_psi_1(concc*x/conca+1)*concc/conca-gsl_sf_psi_1(x/conca+1)/conca);

   return df;
}

   void
funcroot_fdf (double x, void *params, double *y, double *dy)
{

   *y=funcroot(x, NULL);
   *dy=funcroot_deriv(x,NULL);
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  compalpha
 *  Description:  Compute the base priors
 * =====================================================================================
 */
   int
compalpha()
{
   int status;
   int iter = 0, max_iter = 100;
   const gsl_root_fdfsolver_type *T;
   gsl_root_fdfsolver *s;
   double x0, x = 0.1;
   gsl_function_fdf FDF;

   ic=scorethr2/width;

   FDF.f = &funcroot;
   FDF.df = &funcroot_deriv;
   FDF.fdf = &funcroot_fdf;
   FDF.params = NULL;

   T = gsl_root_fdfsolver_newton;
   s = gsl_root_fdfsolver_alloc (T);
   gsl_root_fdfsolver_set (s, &FDF, x);

   do
   {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, 1e-3);
   }
   while (status == GSL_CONTINUE && iter < max_iter);

   alpha=x;
   beta=concc/conca*alpha;


   gsl_root_fdfsolver_free (s);
   return status;
}

   void 
displaymat(vvd & mat)
{
   cout.precision(4);
   for (int i=0;i<4;i++){
      for (unsigned int j=0;j<mat.size();j++){
         cout << mat[j][i] << "\t";
      }
      cout << "\n";
   }
   cout << "\n";
   cout.precision(6);
   return;
}

   void 
matfreqdisp(vvd& matrice)
{
   vd dum(4,0.0);
   vvd mat(width,dum);
   int j=0;
   for (vvd::iterator imat=matrice.begin();imat!=matrice.end();imat++){
      vd & col=*imat;
      double col0=0.3*exp(col[0]);//A
      double col2=0.2*exp(col[2]);//C
      double col3=0.2*exp(col[3]);//G
      double col1=0.3*exp(col[1]);//T
      //convert in A/C/G/T format.
      mat[j][0]=col0;//floor(100*col0);//A
      mat[j][1]=col1;//floor(100*col2);//T
      mat[j][2]=col2;//floor(100*col3);//C
      mat[j][3]=col3;//floor(100*col1);//G
      j++;
   }
   cout << "A ";
   for (vvd::const_iterator imat=mat.begin();imat!=mat.end();imat++){
      cout << (*imat)[0] << " ";
   }
   cout << "\nT ";
   for (vvd::const_iterator imat=mat.begin();imat!=mat.end();imat++){
      cout << (*imat)[1] << " ";
   }
   cout << "\nC ";
   for (vvd::const_iterator imat=mat.begin();imat!=mat.end();imat++){
      cout << (*imat)[2] << " ";
   }
   cout << "\nG ";
   for (vvd::const_iterator imat=mat.begin();imat!=mat.end();imat++){
      cout << (*imat)[3] << " ";
   }
   cout << "\n";
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


