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
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <gsl/gsl_randist.h>

using namespace std;

#include "test_cmdline.h"
#include "test.hpp"
#include "const.hpp"
#include "random.hpp"
#include "motif.hpp"
#include "tree.hpp"
#include "distinfo.hpp"
#include "sequence.hpp"
#include "display.hpp"
#include "genmot.hpp"

test_args_info test_args;
unsigned int numhamm;

   void
evolvesite(Motif & mot)
{

   //number of sites for evolution stat:
   unsigned int numforstat(200);

   ofstream outf("evol_site.dat");

   //distances to ref
   vd distances(nbspecies,0);
   if (species=="droso"){
      distances[0]=0;
      distances[1]=0.0968;
      distances[2]=0.1008;
      distances[3]=0.2334;
      distances[4]=0.2229;
      distances[5]=1.5141;
      distances[6]=1.2973;
      distances[7]=1.3057;
      distances[8]=1.5334;
      distances[9]=1.4741;
      distances[10]=1.5203;
      distances[11]=1.4661;
   }
   else if (species=="eutherian"){
      distances[0]=0; //mus
      distances[1]=0.1761; // rat
      distances[2]=0.5004; //pan
      distances[3]=0.5004; //hom
      distances[4]=0.5003; //gor
      distances[5]=0.5001; //pon
      distances[6]=0.5051; //mac
      distances[7]=0.5117; //cal
      distances[8]=0.5301; //equ
      distances[9]=0.573; //can
      distances[10]=0.5092; //sus
      distances[11]=0.5991; //bos
   }

   evolutionary_model=2;
   inittreedist();
   vd dum(4,0.);
   vvd vdum(4,dum);
   vvd wdum(mot.motwidth,dum);

   double initscore;

   for (unsigned int n=1;n<nbspecies;n++){
      cout << numtospecies(n) << endl;
      for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){

         if (!ivm->matches[n]) continue;

         initscore=scoref(ivm->alignseq[0],mot.matprec);

         vint site(ivm->alignseq[0]);
         vint initsite(site);

         // Initialize evolution matrix to Id
         vvd pm(vdum);
         pm[0][0]=1.;
         pm[1][1]=1.;
         pm[2][2]=1.;
         pm[3][3]=1.;

         // initial probability from ref site
         vvd Pdisp_f(wdum);
         vvd Pdisp_h(wdum);
         vvd Initprob(wdum);
         for (unsigned int i=0;i<mot.motwidth;i++) 
            Initprob[i][ivm->alignseq[0][i]]=1;
         Pdisp_f=Initprob;
         Pdisp_h=Initprob;

         // HALPERN
         for (unsigned int i=0;i<mot.motwidth;i++){
            Pdisp_h[i]=evolvedist_halpern(Initprob[i],mot.matfreq[i],distances[n]);
         }

         // FELSEN
         for (unsigned int i=0;i<mot.motwidth;i++){
            Pdisp_f[i]=evolvedist_felsen(Initprob[i],mot.matfreq[i],distances[n]);
         }
//         cout << "Felsen, distance=" << distances[n] << endl;
//         displaymat(Pdisp_f);
//         cout << "Halpern, distance=" << distances[n] << endl;
//         displaymat(Pdisp_h);


         // SCORES COMPUTATION
         //
         // First, the "true" TFBS
         //

         outf << numtospecies(n) << " ";
         outf << "Data" << " ";
         outf << scoref(ivm->alignseq[n],mot.matprec)-initscore << " ";
         outf << scorefhamming(ivm->alignseq[n],initsite) << endl;


         // draw numforstat sites for evolved sites scores computation
         //cout << "data " << vinttostring(ivm->alignseq[n]) << endl;
         for (unsigned int ks=0;ks<numforstat;ks++){

            vint dumsite(mot.motwidth,4);

            // HALPERN
            vint tempsite(dumsite);
            //while (scorefhamming(tempsite,initsite)>numhamm){
            while (scoref(tempsite,mot.matprec)<mot.motscorethrcons){
               for (unsigned int j=0;j<mot.motwidth;j++){
                  double pdist[4]={Pdisp_h[j][0],Pdisp_h[j][1],Pdisp_h[j][2],Pdisp_h[j][3]};
                  gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
                  unsigned int base=gsl_ran_discrete (gslran,gsldist);
                  gsl_ran_discrete_free (gsldist);
                  site[j]=base;
               }
               tempsite=site;
            }
            outf << numtospecies(n) << " ";
            outf << "Halpern" << " ";
            outf << scoref(site,mot.matprec)-initscore << " ";
            outf << scorefhamming(site,initsite) << endl;

//            if (ks==1)
//            cout << vinttostring(ivm->alignseq[n]) << " " << vinttostring(site) << " (Halpern)" << endl;

            //cout << "halpern " << vinttostring(site) << endl;

            // FELSEN
            tempsite=dumsite;
            //while (scorefhamming(tempsite,initsite) > numhamm){
            while (scoref(tempsite,mot.matprec)<mot.motscorethrcons){
               for (unsigned int j=0;j<mot.motwidth;j++){
                  double pdist[4]={Pdisp_f[j][0],Pdisp_f[j][1],Pdisp_f[j][2],Pdisp_f[j][3]};
                  gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
                  unsigned int base=gsl_ran_discrete (gslran,gsldist);
                  gsl_ran_discrete_free (gsldist);
                  site[j]=base;
               }
               tempsite=site;
            }
            outf << numtospecies(n) << " ";
            outf << "Felsen" << " ";
            outf << scoref(site,mot.matprec)-initscore << " ";
            outf << scorefhamming(site,initsite) << endl;
            
//            if (ks==1)
//            cout << vinttostring(ivm->alignseq[n]) << " " << vinttostring(site) << " (Felsen)" << endl;

         }
         //         getchar();
      }
   }
   
   outf.close();
   outf.open("evol_site_models.dat");
   

   for (double d = 0;d < 2; d += 0.05){
      for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){

         initscore=scoref(ivm->alignseq[0],mot.matprec);

         vint site(ivm->alignseq[0]);
         vint initsite(site);

         // Initialize evolution matrix to Id
         vvd pm(vdum);
         pm[0][0]=1.;
         pm[1][1]=1.;
         pm[2][2]=1.;
         pm[3][3]=1.;

         // initial probability from ref site
         vvd Pdisp_f(wdum);
         vvd Pdisp_h(wdum);
         vvd Initprob(wdum);
         for (unsigned int i=0;i<mot.motwidth;i++) 
            Initprob[i][ivm->alignseq[0][i]]=1;
         Pdisp_f=Initprob;
         Pdisp_h=Initprob;

         // HALPERN
         for (unsigned int i=0;i<mot.motwidth;i++){
            Pdisp_h[i]=evolvedist_halpern(Initprob[i],mot.matfreq[i],d);
         }

         // FELSEN
         for (unsigned int i=0;i<mot.motwidth;i++){
            Pdisp_f[i]=evolvedist_felsen(Initprob[i],mot.matfreq[i],d);
         }
//         cout << "Felsen, distance=" << distances[n] << endl;
//         displaymat(Pdisp_f);
//         cout << "Halpern, distance=" << distances[n] << endl;
//         displaymat(Pdisp_h);


         // SCORES COMPUTATION
         //

         // draw numforstat sites for evolved sites scores computation
         //cout << "data " << vinttostring(ivm->alignseq[n]) << endl;
         for (unsigned int ks=0;ks<numforstat;ks++){

            vint dumsite(mot.motwidth,4);

            // HALPERN
            vint tempsite(dumsite);
            //while (scorefhamming(tempsite,initsite)>numhamm){
            while (scoref(tempsite,mot.matprec)<mot.motscorethrcons){
               for (unsigned int j=0;j<mot.motwidth;j++){
                  double pdist[4]={Pdisp_h[j][0],Pdisp_h[j][1],Pdisp_h[j][2],Pdisp_h[j][3]};
                  gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
                  unsigned int base=gsl_ran_discrete (gslran,gsldist);
                  gsl_ran_discrete_free (gsldist);
                  site[j]=base;
               }
               tempsite=site;
            }
            outf << d << " ";
            outf << "Halpern" << " ";
            outf << scoref(site,mot.matprec)-initscore << " ";
            outf << scorefhamming(site,initsite) << endl;

//            if (ks==1)
//            cout << vinttostring(ivm->alignseq[n]) << " " << vinttostring(site) << " (Halpern)" << endl;

            //cout << "halpern " << vinttostring(site) << endl;

            // FELSEN
            tempsite=dumsite;
            //while (scorefhamming(tempsite,initsite) > numhamm){
            while (scoref(tempsite,mot.matprec)<mot.motscorethrcons){
               for (unsigned int j=0;j<mot.motwidth;j++){
                  double pdist[4]={Pdisp_f[j][0],Pdisp_f[j][1],Pdisp_f[j][2],Pdisp_f[j][3]};
                  gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
                  unsigned int base=gsl_ran_discrete (gslran,gsldist);
                  gsl_ran_discrete_free (gsldist);
                  site[j]=base;
               }
               tempsite=site;
            }
            outf << d << " ";
            outf << "Felsen" << " ";
            outf << scoref(site,mot.matprec)-initscore << " ";
            outf << scorefhamming(site,initsite) << endl;
            
//            if (ks==1)
//            cout << vinttostring(ivm->alignseq[n]) << " " << vinttostring(site) << " (Felsen)" << endl;

         }
         //         getchar();
      }
   }

   outf.close();
}

void
evolvecv(Motif & mot)
{
      
   ofstream outf;
   outf.open("evol_cv.dat");

   Motif refmot = mot;

   unsigned int nshuffle=test_args.nshuffle_arg;
   unsigned int nshufflecons=nshuffle;
   int iter(0);
   int logiter(0);


   // TEST FELSEN CONVERGENCE
   cout << "Testing Felsen convergence..." << endl;
   evolutionary_model=1;
   inittreedist();
   for (unsigned int k=0;k<nshufflecons;k++){
      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      vma tmpseqs;
      logiter=0;

      for (int i=0;i<refmot.seqs.size();i++){
         tmpseqs.push_back(refmot.seqs[i]);

         int numinst=tmpseqs.size(); 
         iter=(int)pow(10,0.1*logiter);
         if (numinst<iter) continue;
         logiter++;

         cout << "# of sites for Felsen: " <<   tmpseqs.size() <<  " (iter=" << k+1 << "/" << nshufflecons << ")" << endl;
         mot.seqs=tmpseqs;
         // MAXIMUM LIKELIHOOD
         mot.compprec();
         mot.matfreq=mattofreq(mot.matprec);
         double distkl=0.;
         for (int pos=0;pos<mot.motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=mot.matfreq[pos][ib]*log(mot.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            }
         }
         outf << "Felsen_colopti " << numinst << " " << distkl << endl;

      }
   }
   
   // TEST HALPERN CONVERGENCE
   cout << "Testing Halpern convergence..." << endl;
   evolutionary_model=2;
   inittreedist();
   for (unsigned int k=0;k<nshufflecons;k++){

      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      vma tmpseqs;
      logiter=0;

      for (int i=0;i<refmot.seqs.size();i++){
         tmpseqs.push_back(refmot.seqs[i]);

         int numinst=tmpseqs.size(); 
         iter=(int)pow(10,0.1*logiter);
         if (numinst<iter) continue;
         logiter++;

         cout << "# of sites for Halpern: " <<   tmpseqs.size() <<  " (iter=" << k+1 << "/" << nshufflecons << ")" << endl;
         mot.seqs=tmpseqs;
         // MAXIMUM LIKELIHOOD
         mot.compprec();
         mot.matfreq=mattofreq(mot.matprec);
         double distkl=0.;
         for (int pos=0;pos<mot.motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=mot.matfreq[pos][ib]*log(mot.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            }
         }
         outf << "Halpern_colopti " << numinst << " " << distkl << endl;
      }
   }
   
   // TEST INDEPENDENT SPECIES CONVERGENCE
   cout << "Testing Indenpendent species convergence..." << endl;
   for (unsigned int k=0;k<nshufflecons;k++){
      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      vma tmpseqs;
      logiter=0;

      for (int i=0;i<refmot.seqs.size();i++){
         tmpseqs.push_back(refmot.seqs[i]);

         int numinst=tmpseqs.size(); 
         iter=(int)pow(10,0.1*logiter);
         if (numinst<iter) continue;
         logiter++;

         cout << "# of sites for Independent species: " <<   tmpseqs.size() <<  " (iter=" << k+1 << "/" << nshufflecons << ")" << endl;
         mot.seqs=tmpseqs;
         mot.compprec_inde();
         mot.matfreq=mattofreq(mot.matprec);
               
         double distkl=0.;
         for (int pos=0;pos<mot.motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=mot.matfreq[pos][ib]*log(mot.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            }
         }
         outf << "Inde " << numinst << " " << distkl << endl;

      }
   }
   
   // TEST CONVERGENCE ON CONSERVED SITES ON EACH SPECIES
   cout << "Testing conserved sites convergence on each species..." << endl;
   for (unsigned int n=0;n<nbspecies;n++)
   {
      for (unsigned int k=0;k<nshuffle;k++){
         // RANDOM SHUFFLING
         random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
         mot.seqs.clear();
         logiter=0;
         for (int i=0;i<refmot.seqs.size();i++){
            if (refmot.seqs[i].matches[n]){
               mot.seqs.push_back(refmot.seqs[i]);
               int numinst=mot.seqs.size(); 
               iter=(int)pow(10,0.1*logiter);
               if (numinst<iter) continue;
               logiter++;

               cout << "# of sites for ConsRef: " <<   mot.seqs.size() << 
                  " (species=" << n+1 << "/" << nbspecies << ", iter=" << k+1 << "/" << nshuffle << ")" << endl;
               mot = comprefmot(mot, n);
               double distkl=0.;
               for (int pos=0;pos<mot.motwidth;pos++){
                  for (int ib=0;ib<4;ib++){
                     distkl+=mot.matfreq[pos][ib]*log(mot.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
                  }
               }
               outf << "ConsRef_" << numtospecies(n) << " " << numinst << " " << distkl << endl;
            }
         }
      }
   }
   
   // TEST CONVERGENCE ON ~Neffx REF CONSERVED SITES
   cout << "Testing Neff x conserved sites convergence on ref species..." << endl;
   unsigned int Neff;
   if (species=="droso") Neff=4;
   else if (species=="eutherian") Neff=2;
   int n = 0;
   for (unsigned int k=0;k<nshufflecons;k++){
      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      mot.seqs.clear();
      logiter=0;
      for (int i=0;i<refmot.seqs.size();i++){
         if (refmot.seqs[i].matches[n]){
            for (unsigned int nsite=0;nsite<Neff;nsite++){
               mot.seqs.push_back(refmot.seqs[i]);
            }
            int numinst=mot.seqs.size() / Neff; 
            iter=(int)pow(10,0.1*logiter);
            if (numinst<iter) continue;
            logiter++;

            cout << "# of sites for Control: " <<   mot.seqs.size() / Neff << "x" << Neff <<
               ", iter=" << k+1 << "/" << nshufflecons << ")" << endl;
            mot = comprefmot(mot, n);
            double distkl=0.;
            for (int pos=0;pos<mot.motwidth;pos++){
               for (int ib=0;ib<4;ib++){
                  distkl+=mot.matfreq[pos][ib]*log(mot.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
               }
            }
            outf << "Control " << numinst << " " << distkl << endl;
         }
      }
   }

   outf.close();


   return;
}

void
evolvecvcrossval(Motif & mot)
{
      
   ofstream outf;
   outf.open("evol_cv.dat");

   Motif refmot = mot;

   unsigned int nshuffle=test_args.nshuffle_arg;
   unsigned int nshufflecons=nshuffle;
   unsigned int Nhalf = floor(refmot.seqs.size() / 2);
   unsigned int N = refmot.seqs.size();
   int iter(0);
   int logiter(0);


   // TEST FELSEN CONVERGENCE
   cout << "Testing Felsen convergence..." << endl;
   evolutionary_model=1;
   inittreedist();
   for (unsigned int k=0;k<nshufflecons;k++){
      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      vma trainseqs, testseqs;
      vvd trainmat, testmat;
      logiter=0;
      
      for (int i = Nhalf; i < N; i++){
         testseqs.push_back(refmot.seqs[i]);
         mot.seqs=testseqs;
         mot = comprefmot(mot, 0);
         testmat = mattofreq(mot.matprec);
      }

      for (int i = 0; i < Nhalf; i++){
         trainseqs.push_back(refmot.seqs[i]);

         int numinst=trainseqs.size(); 
         iter=(int)pow(10,0.1*logiter);
         if (numinst<iter) continue;
         logiter++;

         cout << "# of sites for Felsen: " <<   trainseqs.size() <<  " (iter=" << k+1 << "/" << nshufflecons << ")" << endl;
         // MAXIMUM LIKELIHOOD
         mot.seqs=trainseqs;
         mot.compprec();
         trainmat = mattofreq(mot.matprec);
               
         double distkl=0.;
         for (int pos=0;pos<mot.motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl += trainmat[pos][ib]*log(trainmat[pos][ib] / testmat[pos][ib]);
            }
         }
         outf << "Felsen_colopti " << numinst << " " << distkl << endl;

      }
   }
   
   // TEST HALPERN CONVERGENCE
   cout << "Testing Halpern convergence..." << endl;
   inittreedist();
   for (unsigned int k=0;k<nshufflecons;k++){
      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      vma trainseqs, testseqs;
      vvd trainmat, testmat;
      logiter=0;
      
      for (int i = Nhalf; i < N; i++){
         testseqs.push_back(refmot.seqs[i]);
         mot.seqs=testseqs;
         mot = comprefmot(mot, 0);
         testmat = mattofreq(mot.matprec);
      }

      for (int i = 0; i < Nhalf; i++){
         trainseqs.push_back(refmot.seqs[i]);

         int numinst=trainseqs.size(); 
         iter=(int)pow(10,0.1*logiter);
         if (numinst<iter) continue;
         logiter++;

         cout << "# of sites for Halpern: " <<   trainseqs.size() <<  " (iter=" << k+1 << "/" << nshufflecons << ")" << endl;
         // MAXIMUM LIKELIHOOD
         mot.seqs=trainseqs;
         mot.compprec();
         trainmat = mattofreq(mot.matprec);
               
         double distkl=0.;
         for (int pos=0;pos<mot.motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl += trainmat[pos][ib]*log(trainmat[pos][ib] / testmat[pos][ib]);
            }
         }
         outf << "Halpern_colopti " << numinst << " " << distkl << endl;

      }
   }
   
   // TEST INDEPENDENT SPECIES CONVERGENCE
   cout << "Testing Indenpendent species convergence..." << endl;
   for (unsigned int k=0;k<nshufflecons;k++){
      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      vma trainseqs, testseqs;
      vvd trainmat, testmat;
      logiter=0;

      for (int i = Nhalf; i < N; i++){
         testseqs.push_back(refmot.seqs[i]);
         mot.seqs=testseqs;
         mot = comprefmot(mot, 0);
         testmat = mattofreq(mot.matprec);
      }

      for (int i = 0; i < Nhalf; i++){
         trainseqs.push_back(refmot.seqs[i]);

         int numinst=trainseqs.size(); 
         iter=(int)pow(10,0.1*logiter);
         if (numinst<iter) continue;
         logiter++;

         cout << "# of sites for Independent species: " <<   trainseqs.size() <<  " (iter=" << k+1 << "/" << nshufflecons << ")" << endl;
         
         mot.seqs=trainseqs;
         mot.compprec_inde();
         trainmat = mattofreq(mot.matprec);
               
         double distkl=0.;
         for (int pos=0;pos<mot.motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl += trainmat[pos][ib]*log(trainmat[pos][ib] / testmat[pos][ib]);
            }
         }
         outf << "Inde " << numinst << " " << distkl << endl;

      }
   }
   
   
   // TEST CONVERGENCE ON CONSERVED SITES ON EACH SPECIES
   cout << "Testing conserved sites convergence on each species..." << endl;
   for (unsigned int n=0;n<nbspecies;n++)
   {
      for (unsigned int k=0;k<nshuffle;k++){
         // RANDOM SHUFFLING
         random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
         vma trainseqs, testseqs;
         vvd trainmat, testmat;
         
         logiter=0;
         for (int i = Nhalf; i < N; i++){
            testseqs.push_back(refmot.seqs[i]);
            mot.seqs=testseqs;
            mot = comprefmot(mot, 0);
            testmat = mattofreq(mot.matprec);
         }

         for (int i = 0; i < Nhalf; i++){
            if (refmot.seqs[i].matches[n]){
               trainseqs.push_back(refmot.seqs[i]);
               
               int numinst=trainseqs.size(); 
               iter=(int)pow(10,0.1*logiter);
               if (numinst<iter) continue;
               logiter++;

               cout << "# of sites for ConsRef: " <<   trainseqs.size() << 
                  " (species=" << n+1 << "/" << nbspecies << ", iter=" << k+1 << "/" << nshuffle << ")" << endl;
         
               mot.seqs=trainseqs;
               mot = comprefmot(mot, n);
               trainmat = mattofreq(mot.matprec);
               
               double distkl=0.;
               for (int pos=0;pos<mot.motwidth;pos++){
                  for (int ib=0;ib<4;ib++){
                     distkl += trainmat[pos][ib]*log(trainmat[pos][ib] / testmat[pos][ib]);
                  }
               }
               outf << "ConsRef_" << numtospecies(n) << " " << numinst << " " << distkl << endl;
            }
         }
      }
   }
   
   // TEST CONVERGENCE ON ~Neffx REF CONSERVED SITES
   cout << "Testing Neff x conserved sites convergence on ref species..." << endl;
   unsigned int Neff;
   if (species=="droso") Neff=4;
   else if (species=="eutherian") Neff=2;
   int n = 0;
   for (unsigned int k=0;k<nshufflecons;k++){
      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      vma trainseqs, testseqs;
      vvd trainmat, testmat;
      
      logiter=0;
      for (int i = Nhalf; i < N; i++){
         testseqs.push_back(refmot.seqs[i]);
         mot.seqs=testseqs;
         mot = comprefmot(mot, 0);
         testmat = mattofreq(mot.matprec);
      }

      for (int i = 0; i < Nhalf; i++){
         if (refmot.seqs[i].matches[n]){
            for (unsigned int nsite=0;nsite<Neff;nsite++){
               trainseqs.push_back(refmot.seqs[i]);
            }
            
            int numinst = trainseqs.size() / Neff; 
            iter=(int)pow(10,0.1*logiter);
            if (numinst<iter) continue;
            logiter++;

            cout << "# of sites for Control: " <<   trainseqs.size() / Neff << "x" << Neff <<
               ", iter=" << k+1 << "/" << nshufflecons << ")" << endl;
            mot.seqs=trainseqs;
            mot = comprefmot(mot, n);
            trainmat = mattofreq(mot.matprec);
            
            double distkl=0.;
            for (int pos=0;pos<mot.motwidth;pos++){
               for (int ib=0;ib<4;ib++){
                  distkl += trainmat[pos][ib]*log(trainmat[pos][ib] / testmat[pos][ib]);
               }
            }
            outf << "Control " << numinst << " " << distkl << endl;
         }
      }
   }

   outf.close();


   return;
}

   void
refinemotif(Motif & mot, vseq & align)
{
   mot.setscorethr2meaninfo();
   //mot.matinithamming(mot.motscorethr2, numhamm);
   mot.matinit(mot.motscorethr2);

   // refining motif on ref sequences
   scorethr2 = mot.motscorethr2;
   width = mot.motwidth;
   if (test_args.nops_given){
      alpha = 1e-2;
      beta = concc / conca * alpha;
   } else {
      compalpha();
   }
   mot = comprefmot(mot, 0);

   // deleting unessential flanking bases
   mot.cutflanking();
   
   mot.setscorethr2meaninfo();
   mot.matinit(mot.motscorethr2);
   scorethr2 = mot.motscorethr2;
   width = mot.motwidth;
   if (test_args.nops_given){
      alpha = 1e-2;
      beta = concc / conca * alpha;
   } else {
      compalpha();
   }
   mot = comprefmot(mot, 0);

   //mot.matinithamming(mot.motscorethr2, numhamm);

   return;
}

   void
test_args_init()
{
   if (!strcmp(test_args.species_arg, "droso")) {
      species = "droso";
      nbspecies = 12;
      conca = 0.3;
   } else if (!strcmp(test_args.species_arg, "eutherian")) {
      species = "eutherian";
      nbspecies = 12;
      conca = 0.263;
   }
   concc = 0.5 - conca;
   conct = conca;
   concg = concc;
   width = test_args.width_arg;
   scorethr2 = test_args.threshold_arg * log(2);
   scorethr = scorethr2 * (1 - 2.0 / width);
   scorethrcons = scorethr2 * (1 - 1.0 / width);
   evolutionary_model = test_args.evolutionary_model_arg;
   numhamm = test_args.numhamm_arg;
   if (test_args.progress_given) progress = true;
}

string test_datapath;

/**
 * ===  FUNCTION  =============================================================
 *         Name:  cmd_test
 *  Description:  Run tests
 * ============================================================================
 */
   int
cmd_test(int argc, char ** argv)
{
   if (test_cmdline_parser(argc, argv, & test_args) != 0)
      exit(EXIT_FAILURE);
   test_args_init();
   rnginit();
   const char * imo_test_datapath = getenv("IMOGENE_DATA");
   if (imo_test_datapath == NULL) {
      test_datapath = DATA_PATH;
   } else {
      test_datapath = imo_test_datapath;
   }

   cout << "Loading alignments " << endl;
   regints = loadseqs(test_args.align_arg);
   cout << "Nb sequences to scan: " << regints.size() << endl;
       
   cout << "Loading Motifs" << endl;
   vmot mots;
   if (test_args.motifs_ATCG_given){
      loadmotsATCG(test_args.motifs_ATCG_arg, mots);
   } else if (test_args.motifs_given){
      loadmots(test_args.motifs_arg, mots);
   } else {
      exit(EXIT_FAILURE);
   }
   
   Motif mot;
   vmot logomots;
   
   mot = mots[0];
   logomots.push_back(mot);
   
   refinemotif(mot, regints);
   logomots.push_back(mot);
   
   // display
   ofstream motmeldb("motifs.txt");
   mot.display(motmeldb);
   motmeldb.close();
   dispweblogo(logomots);

   cout << "Output conserved sites..." << endl;
   ofstream sites("sites.txt");
   for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){
       int i = 0;
       for (ivint imat = ivm->matches.begin(); imat != ivm->matches.end(); imat++) {
           if (*imat) {
               sites << vinttostring(ivm->alignseq[i]) << "\n";
           }
           i++;
       }
      sites << "\n";
   }
   sites.close();

   if (test_args.evolve_site_given)
   {
      cout << "Testing site conservation..." << endl;
      evolvesite(mot);
   }
   else if (test_args.evolve_cv_given)
   {
      cout << "Testing convergence..." << endl;
      if (test_args.crossval_given) evolvecvcrossval(mot);
      else evolvecv(mot);
   }

   gsl_rng_free(gslran);
   test_cmdline_parser_free(&test_args);
   return EXIT_SUCCESS;
}		/* -----  end of function test  ----- */


