
#include <cmath>
#include <cstdlib>
#include<algorithm>

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


#include "const.hpp"
#include "vectortypes.hpp"
#include "random.hpp"
#include "motif.hpp"
#include "sequence.hpp"
#include "scangen.hpp"
#include "tree.hpp"
#include "montecarlo.hpp"

using namespace std;


Motif::Motif()
{
   for (unsigned int i=0;i<distwidth;i++){
      distmot[i]=0;
   }
   name="";
   id="";
   nbmatch=0;
   nbmatchback=0;
   lambda=0;
   lambdatrain=0;
   pvalue=0;
   scorepoiss=0;
   meanval=0;
   nbmot=0;
   ntrain=0;
   check=true;
   vinst dumvinst;
   instances=vvinst(nbchrom,dumvinst);

   motscorethr2=scorethr2;
   motwidth=width;
   motscorethr=scorethr2-2*(double)motwidth/10;
   motscorethrcons=scorethr2-(double)motwidth/10;

   optauc=0.5;
}


TFBS::TFBS()
{
   num=0;
   numrand=0;
   probref=0;
   probini=0;
   prob=0;
   probmix=0;
}

// recursive function
// enumerate all possible TFBS above a motscorethr2 threshold
// n begins at motwidth and decreases to 0
// site is the temporary TFBS
// sites are all accepted sites
void
Motif::enumeratewthres(int n,vint& site,vtfbs& sites)
{
   if (n==0)
   {
      if (scoref(site,matprec)>motscorethr2){
         TFBS bs;
         bs.site=site;
         bs.score=scoref(site,matprec);
         bs.prob=compprob(site,matfreq);
         sites.push_back(bs);
      }
   }
   else
   {
         for(int j=0;j<4;j++)
         {
            vint sitetemp(site);
            if (n==motwidth)
            {
               sitetemp.clear();
               sitetemp.push_back(j);
            }
            else
            {
               sitetemp.push_back(j);
            }
            enumeratewthres(n-1,sitetemp,sites);
         }
   }
   return;
}

// recursive function
// enumerate all possible TFBS
// n begins at motwidth and decreases to 0
// site is the temporary TFBS
// sites are all accepted sites
double totprob(0),totini(0),totref(0),totmix(0);
void enumerateallsitesncomp(int n,vint& site,vtfbs& sites,Motif & motini,Motif & refmot,vmot & motmix,vd & kprobs,ofstream & outf)
{
   if (n==0)
   {
      TFBS bs;
      bs.site=site;
      for (ivvint ivv=refmot.refinstances.iseqs.begin();ivv!=refmot.refinstances.iseqs.end();ivv++){
         if (*ivv==site) bs.num++;
      }
      bs.prob=(double)bs.num/refmot.refinstances.seqs.size();
      
      // REFERENCE
      bs.score=scoref(site,refmot.matprec);
      bs.probref=compprob(site,refmot.matfreq);
      // INITIAL
      bs.probini=compprob(site,motini.matfreq);

      bs.probmix=0;
      for (unsigned int j=0;j<kprobs.size();j++){
         double sc=1;
         unsigned int pos=0;
         for (ivint iv=site.begin();iv!=site.end();iv++){
            sc*=motmix[j].matfreq[pos][*iv];
            pos++;
         }
         bs.probmix+=kprobs[j]*sc;
      }
      sites.push_back(bs);
      
      outf << vinttostring(site) << "\t";
      outf << bs.prob << "\t";
      outf << bs.probini << "\t";
      outf << bs.probref << "\t";
      outf << bs.probmix << "\n";
      totprob+=bs.prob;
      totini+=bs.probini;
      totref+=bs.probref;
      totmix+=bs.probmix;
   }
   else
   {
      for(int j=0;j<4;j++)
      {
         vint sitetemp(site);
         // "insert" for lexicographic order + letter in alphabetical order
         unsigned int base;
         if (j==0) base=0;
         if (j==1) base=2;
         if (j==2) base=3;
         if (j==3) base=1;
         if (n==refmot.motwidth)
         {
            sitetemp.clear();
            //sitetemp.push_back(j);
            sitetemp.insert(sitetemp.begin(),base);
         }
         else
         {
            //sitetemp.push_back(j);
            sitetemp.insert(sitetemp.begin(),base);
         }
         enumerateallsitesncomp(n-1,sitetemp,sites,motini,refmot,motmix,kprobs,outf);

      }
   }
   return;
}

void enumeratesitesforkmeans(int n,vint& site,vtfbs& sites,Motif & motini,Motif & refmot,vvmot & vmotmix,vvd & vkprobs,ofstream & outf)
{
   if (n==0)
   {
      TFBS bs;
      bs.site=site;
      for (ivvint ivv=refmot.refinstances.iseqs.begin();ivv!=refmot.refinstances.iseqs.end();ivv++){
         if (*ivv==site) bs.num++;
      }
      bs.prob=(double)bs.num/refmot.refinstances.seqs.size();
      
      // REFERENCE
      bs.score=scoref(site,refmot.matprec);
      bs.probref=compprob(site,refmot.matfreq);
      // INITIAL
      bs.probini=compprob(site,motini.matfreq);
      
      outf << vinttostring(site) << "\t";
      outf << bs.prob << "\t";
      outf << bs.probini << "\t";
      outf << bs.probref << "\t";

      for (unsigned int k=0;k<vkprobs.size();k++){
         bs.probmix=0;
         for (unsigned int j=0;j<vkprobs[k].size();j++){
            double sc=1;
            unsigned int pos=0;
            for (ivint iv=site.begin();iv!=site.end();iv++){
               sc*=vmotmix[k][j].matfreq[pos][*iv];
               pos++;
            }
            bs.probmix+=vkprobs[k][j]*sc;
         }
         outf << bs.probmix << "\t";
      }
      
      outf << "\n";
      
      sites.push_back(bs);
   }
   else
   {
      for(int j=0;j<4;j++)
      {
         vint sitetemp(site);
         // "insert" for lexicographic order + letter in alphabetical order
         unsigned int base;
         if (j==0) base=0;
         if (j==1) base=2;
         if (j==2) base=3;
         if (j==3) base=1;
         if (n==refmot.motwidth)
         {
            sitetemp.clear();
            //sitetemp.push_back(j);
            sitetemp.insert(sitetemp.begin(),base);
         }
         else
         {
            //sitetemp.push_back(j);
            sitetemp.insert(sitetemp.begin(),base);
         }
         enumeratesitesforkmeans(n-1,sitetemp,sites,motini,refmot,vmotmix,vkprobs,outf);

      }
   }
   return;
}

   vvvd 
Motif::correlations()
{
   vd dum(16,0.); // 2-bases correlation
   vvd m(motwidth,dum);
   vvvd Mcorr(motwidth,m);
   vvvd sxy(Mcorr),sx(Mcorr),sy(Mcorr);

   vd ddum(4,0.);
   vvd mean(motwidth,ddum);

   // CONS
   //matinit(mot.motscorethr2);
   // ALL
   comprefinstances(regints,0);
   comprefmot();

   cout << "There are " << refinstances.iseqs.size() << " sites for computation."<< endl;

   // First compute the means
   unsigned int basex,basey,linpos;
   double xi,yi;
   for (ivvint iv=refinstances.iseqs.begin();iv!=refinstances.iseqs.end();iv++){
      vint site=*iv;
      //cout << site << endl;
      for (unsigned int posx=0;posx<motwidth;posx++){
         basex=site[posx];
         mean[posx][basex]+=1./refinstances.iseqs.size();
      }
   }
   displaymat(mean);

   for (ivvint iv=refinstances.iseqs.begin();iv!=refinstances.iseqs.end();iv++){
      vint site=*iv;
      for (unsigned int posx=0;posx<motwidth;posx++){
         basex=site[posx];
         for (unsigned int posy=posx;posy<motwidth;posy++){
            basey=site[posy];
            for (unsigned int bx=0;bx<4;bx++){
               for (unsigned int by=0;by<4;by++){
                  linpos=4*bx+by;
                  if (bx==basex) xi=1;
                  else xi=0;
                  if (by==basey) yi=1;
                  else yi=0;
                  sxy[posy][posx][linpos]+=xi*yi;
                  //sxy[posy][posx][linpos]+=(xi-mean[posx][bx])*(yi-mean[posy][by]);
                  //sx[posy][posx][linpos]+=pow(xi-mean[posx][bx],2);
                  //sy[posy][posx][linpos]+=pow(yi-mean[posy][by],2);
               }
            }
         }
      }
   }

   for (unsigned int posx=0;posx<motwidth;posx++){
      for (unsigned int posy=posx;posy<motwidth;posy++){
         for (unsigned int bx=0;bx<4;bx++){
            for (unsigned int by=0;by<4;by++){
               linpos=4*bx+by;
               //double SX=sx[posy][posx][linpos];
               //double SY=sy[posy][posx][linpos];
               double SXY=sxy[posy][posx][linpos];
               double M;

               //if (SX==0 || SY==0){
               //M=0;
               //} else {
               //M=SXY/sqrt(SX*SY);
               M=SXY/refinstances.iseqs.size()-mean[posx][bx]*mean[posy][by];
               if (abs(M)<1e-6) M=0;
               //}
               //Mcorr[posy][posx][linpos]=sxy[posy][posx][linpos]/refinstances.iseqs.size()-mean[posx][bx]*mean[posy][by];
               if (abs(M)>=0.05 && posx!=posy) cout << posx << " " << posy << " " << inttostring(bx) << " " << inttostring(by) << " " << M << endl;
               //Mcorr[posy][posx][linpos]=M;
            }
         }
         //cout << posx << " " << posy << " " <<  Mcorr[posy][posx] << endl;
      }
   }

   return Mcorr;
}

   void 
Motif::testpseudocount(ofstream & outf)
{
   unsigned int numloop=10;

   findinstances(regints);

   matfreq=mattofreq(matprec);
   matinit(motscorethr2);
   vma initseqs=seqs;

   Motif motinit=*this;

   // COMPUTE ALLREF PWM
   cout << "Computing reference motif on reference species instances..." << endl;
   Motif refmot=motinit;
   double dkl=1.;
   cout << "Convergence: ";
   //   refmot.comprefinstances(regints,0);
   cout << "Using conserved instances for PWM" << endl;
   refmot.comprefinstancescons(0);
   int maxloop(50),countloop(1);
   while(dkl>1e-6 && countloop<maxloop)
   {
      vvd prevmat=refmot.matfreq;
      refmot.comprefmot();
      //      refmot.comprefinstances(regints,0);
      dkl=0.;
      for (int pos=0;pos<motwidth;pos++){
         for (int ib=0;ib<4;ib++){
            dkl+=refmot.matfreq[pos][ib]*log(refmot.matfreq[pos][ib]/prevmat[pos][ib]);
         }
      }
      prevmat=refmot.matprec;
      refmot.comprefinstancescons(0);
      cout << refmot.refinstances.iseqs.size() << " ";
      cout.flush();
      countloop++;
   }
   cout << endl;
   if (countloop==maxloop) cout << "Did not converge. Taking last PWM as reference." << endl;
   else cout << "Converged! " << endl;
   refmot.comprefinstances(regints,0);
   cout << "# of sites for AllRef: " <<   refmot.refinstances.iseqs.size() << endl;
   //   refmot.matinit(motscorethr2);

   // DEFINING CONSERVED PWM AS REFERENCE
   refmot.matfreq=mattofreq(refmot.matprec);
   cout << "# of conserved instances: " <<   refmot.seqs.size() << endl;

   int iter(0);
   int logiter(0);

   // TEST CONVERGENCE ON CONSERVED SITES
   cout << "Testing conserved sites convergence on each species..." << endl;
   for (unsigned int k=0;k<numloop;k++){
      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      refinstances.seqs.clear();
      refinstances.iseqs.clear();
      logiter=0;
      for (int i=0;i<refmot.seqs.size();i++){
         if (refmot.seqs[i].matches[0]){
            refinstances.iseqs.push_back(refmot.seqs[i].alignseq[0]);
            refinstances.seqs.push_back(vinttostring(refmot.seqs[i].alignseq[0]));

            int numinst=refinstances.seqs.size(); 
            iter=(int)pow(10,0.1*logiter);
            if (numinst<iter) continue;
            logiter++;

            cout << "# of sites: " <<   refinstances.iseqs.size() << 
               " (iter=" << k+1 << "/" << numloop << ")" << endl;
            comprefmot();
            double distkl=0.;
            for (int pos=0;pos<motwidth;pos++){
               for (int ib=0;ib<4;ib++){
                  distkl+=matfreq[pos][ib]*log(matfreq[pos][ib]/refmot.matfreq[pos][ib]);
               }
            }
            outf << "ConsRef " << alpha << " " << numinst << " " << distkl  << endl;
         }
      }
   }

   // TEST FELSEN CONVERGENCE
   //   /*
   cout << "Testing Felsen convergence..." << endl;
   evolutionary_model=1;
   inittreedist();
   for (unsigned int k=0;k<numloop;k++){
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

         cout << "# of sites for Felsen: " <<   tmpseqs.size() <<  " (iter=" << k+1 << "/" << numloop << ")" << endl;
         seqs=tmpseqs;

         // MEAN POSTERIOR
         //compprec_MCMC();
         compprec();
         matfreq=mattofreq(matprec);
         double distkl=0.;
         for (int pos=0;pos<motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=matfreq[pos][ib]*log(matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            }
         }
         outf << "Felsen_colmean " << alpha << " " << numinst << " " << distkl << endl;
      }
   }


   return;
}

   void 
Motif::testpwmcv(ofstream & outf,ofstream & motmeldb)
{

   findinstances(regints);

   matfreq=mattofreq(matprec);
   matinit(motscorethr2);
   vma initseqs=seqs;

   Motif motinit=*this;

   // COMPUTE ALLREF PWM
   cout << "Computing reference motif on reference species instances..." << endl;
   Motif refmot=motinit;
   double dkl=1.;
   cout << "Convergence: ";
   //cout << "Using conserved instances for PWM" << endl;
   cout << "Using all instances for PWM" << endl;
   refmot.comprefinstances(regints,0);
   //refmot.comprefinstancescons(0);
   while(dkl>1e-6)
   {
      vvd prevmat=refmot.matfreq;
      refmot.comprefmot();
      dkl=0.;
      for (int pos=0;pos<motwidth;pos++){
         for (int ib=0;ib<4;ib++){
            dkl+=refmot.matfreq[pos][ib]*log(refmot.matfreq[pos][ib]/prevmat[pos][ib]);
         }
      }
      prevmat=refmot.matprec;
      //refmot.comprefinstancescons(0);
      refmot.comprefinstances(regints,0);
      cout << refmot.refinstances.iseqs.size() << " ";
      cout.flush();
   }
   cout << endl;
   cout << "Converged! " << endl;
   refmot.comprefinstances(regints,0);
   cout << "# of sites for AllRef: " <<   refmot.refinstances.iseqs.size() << endl;
   //   refmot.matinit(motscorethr2);

   cout << "DEFINING ALL SITES PWM AS REFERENCE" << endl;
   //cout << "DEFINING CONSERVED PWM AS REFERENCE" << endl;
   refmot.matfreq=mattofreq(refmot.matprec);
   cout << "# of conserved instances: " <<   refmot.seqs.size() << endl;
   refmot.display(motmeldb);
   
   refmot.comprefinstancescons(0);
   refmot.comprefmot();
   refmot.matfreq=mattofreq(refmot.matprec);
   refmot.display(motmeldb);
   //exit(9);
   
   refmot.comprefinstances(regints,0);
   refmot.comprefmot();
   
   unsigned int nshuffle=40;
   unsigned int nshufflecons=40;
   int iter(0);
   int logiter(0);
   
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
         seqs=tmpseqs;
         // MAXIMUM LIKELIHOOD
         compprec();
         matfreq=mattofreq(matprec);
         double distkl=0.;
         for (int pos=0;pos<motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=matfreq[pos][ib]*log(matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            }
         }
         outf << "Halpern_colopti " << numinst << " " << distkl << endl;

         // MEAN POSTERIOR
         compprec_MCMC();
         matfreq=mattofreq(matprec);
         distkl=0.;
         for (int pos=0;pos<motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=matfreq[pos][ib]*log(matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            }
         }
         outf << "Halpern_colmean " << numinst << " " << distkl << endl;
      }
   }
   display(motmeldb);

   // TEST FELSEN CONVERGENCE
   //   /*
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
         seqs=tmpseqs;
         // MAXIMUM LIKELIHOOD
         compprec();
         matfreq=mattofreq(matprec);
         double distkl=0.;
         for (int pos=0;pos<motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=matfreq[pos][ib]*log(matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            }
         }
         outf << "Felsen_colopti " << numinst << " " << distkl << endl;

         // MEAN POSTERIOR
         compprec_MCMC();
         matfreq=mattofreq(matprec);
         distkl=0.;
         for (int pos=0;pos<motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=matfreq[pos][ib]*log(matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            }
         }
         outf << "Felsen_colmean " << numinst << " " << distkl << endl;
      }
   }
   display(motmeldb);
   // */
   
   // TEST FELSEN W/ MODIFIED DISTANCE CONVERGENCE
   //   /*
   cout << "Testing Felsen convergence w/ modified distance..." << endl;
   evolutionary_model=1;
   inittreedist();
   for (ivnoe iv=treedist.begin();iv!=treedist.end();iv++){
      //iv->prox1/=10.;
      //iv->prox2/=10.;
      iv->prox1/=2.5;
      iv->prox2/=2.5;
   }
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
         seqs=tmpseqs;
         // MAXIMUM LIKELIHOOD
         compprec();
         matfreq=mattofreq(matprec);
         double distkl=0.;
         for (int pos=0;pos<motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=matfreq[pos][ib]*log(matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            }
         }
         outf << "Felsen_newdist_colopti " << numinst << " " << distkl << endl;

         // MEAN POSTERIOR
         compprec_MCMC();
         matfreq=mattofreq(matprec);
         distkl=0.;
         for (int pos=0;pos<motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=matfreq[pos][ib]*log(matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            }
         }
         outf << "Felsen_newdist_colmean " << numinst << " " << distkl << endl;
      }
   }
   display(motmeldb);
   // */

   // TEST ALLREF CONVERGENCE
   Sequence initRefInstances=refmot.refinstances;
   // /*
   cout << "Testing all sites convergence on reference species..." << endl;
   for (unsigned int k=0;k<nshuffle;k++)
   {
      // RANDOM SHUFFLING
      random_shuffle(initRefInstances.iseqs.begin(),initRefInstances.iseqs.end());
      refinstances.iseqs.clear();
      refinstances.seqs.clear();
      logiter=0;
      for (ivvint iv=initRefInstances.iseqs.begin();iv!=initRefInstances.iseqs.begin()+refmot.seqs.size();iv++){
         refinstances.iseqs.push_back(*iv);
         refinstances.seqs.push_back(vinttostring(*iv));

         // 0.1 step in log scale
         int numinst=refinstances.seqs.size(); 
         iter=(int)pow(10,0.1*logiter);
         if (numinst<iter) continue;
         logiter++;

         cout << "# of sites for AllRef: " <<   refinstances.iseqs.size() << 
            " (iter=" << k+1 << "/" << nshuffle << ")" << endl;
         comprefmot();
         double distkl=0.;
         for (int pos=0;pos<motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=matfreq[pos][ib]*log(matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            }
         }
         outf << "AllRef " << numinst << " " << distkl << endl;
      }
   }
   display(motmeldb);
   //*/

   // TEST CONVERGENCE ON CONSERVED SITES ON EACH SPECIES
   ///*
   cout << "Testing conserved sites convergence on each species..." << endl;
   for (unsigned int n=0;n<nbspecies;n++)
   {
      for (unsigned int k=0;k<nshuffle;k++){
         // RANDOM SHUFFLING
         random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
         refinstances.seqs.clear();
         refinstances.iseqs.clear();
         logiter=0;
         for (int i=0;i<refmot.seqs.size();i++){
            if (refmot.seqs[i].matches[n]){
               refinstances.iseqs.push_back(refmot.seqs[i].alignseq[n]);
               refinstances.seqs.push_back(vinttostring(refmot.seqs[i].alignseq[n]));

               int numinst=refinstances.seqs.size(); 
               iter=(int)pow(10,0.1*logiter);
               if (numinst<iter) continue;
               logiter++;

               cout << "# of sites for ConsRef: " <<   refinstances.iseqs.size() << 
                  " (species=" << n+1 << "/" << nbspecies << ", iter=" << k+1 << "/" << nshuffle << ")" << endl;
               comprefmot();
               double distkl=0.;
               for (int pos=0;pos<motwidth;pos++){
                  for (int ib=0;ib<4;ib++){
                     distkl+=matfreq[pos][ib]*log(matfreq[pos][ib]/refmot.matfreq[pos][ib]);
                  }
               }
               outf << "ConsRef_" << numtospecies(n) << " " << numinst << " " << distkl << endl;
            }
         }
         display(motmeldb);
      }
   }
   //*/
   //


   // TEST CONVERGENCE ON ~Neffx REF CONSERVED SITES (CONTROL: IS CONSERVATION ONLY COUNTING FASTER?)
   ///*
   cout << "Testing conserved sites convergence on each species..." << endl;
   unsigned int Neff;
   if (species==1) Neff=4;
   else if (species==2) Neff=2;
   for (unsigned int k=0;k<nshufflecons;k++){
      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      refinstances.seqs.clear();
      refinstances.iseqs.clear();
      logiter=0;
      for (int i=0;i<refmot.seqs.size();i++){
         if (refmot.seqs[i].matches[0]){

            for (unsigned int nsite=0;nsite<Neff;nsite++){
               refinstances.iseqs.push_back(refmot.seqs[i].alignseq[0]);
               refinstances.seqs.push_back(vinttostring(refmot.seqs[i].alignseq[0]));
            }

            int numinst=refinstances.seqs.size()/Neff; 
            iter=(int)pow(10,0.1*logiter);
            if (numinst<iter) continue;
            logiter++;

            cout << "# of sites for Control: " <<   refinstances.iseqs.size()/Neff << "x" << Neff << 
               ", iter=" << k+1 << "/" << nshufflecons << ")" << endl;
            comprefmot();
            double distkl=0.;
            for (int pos=0;pos<motwidth;pos++){
               for (int ib=0;ib<4;ib++){
                  distkl+=matfreq[pos][ib]*log(matfreq[pos][ib]/refmot.matfreq[pos][ib]);
               }
            }
            outf << "Control " << numinst << " " << distkl << endl;
         }
      }
      display(motmeldb);
   }
   //*/



   return;
}

void
Motif::testpwmcvtimer(ofstream & outf)
{
   clock_t start,finish;
   double dif;

   findinstances(regints);

   matfreq=mattofreq(matprec);
   matinit(motscorethr2);
   vma initseqs=seqs;

   Motif motinit=*this;

   // COMPUTE ALLREF PWM
   cout << "Computing reference motif on reference species instances..." << endl;
   Motif refmot=motinit;
   double dkl=1.;
   cout << "Convergence: ";
   //   refmot.comprefinstances(regints,0);
   cout << "Using conserved instances for PWM" << endl;
   refmot.comprefinstancescons(0);
   while(dkl>1e-6)
   {
      vvd prevmat=refmot.matfreq;
      refmot.comprefmot();
      //      refmot.comprefinstances(regints,0);
      dkl=0.;
      for (int pos=0;pos<motwidth;pos++){
         for (int ib=0;ib<4;ib++){
            dkl+=refmot.matfreq[pos][ib]*log(refmot.matfreq[pos][ib]/prevmat[pos][ib]);
         }
      }
      prevmat=refmot.matprec;
      refmot.comprefinstancescons(0);
      cout << refmot.refinstances.iseqs.size() << " ";
      cout.flush();
   }
   cout << endl;
   cout << "Converged! " << endl;
   refmot.comprefinstances(regints,0);
   cout << "# of sites for AllRef: " <<   refmot.refinstances.iseqs.size() << endl;
   //   refmot.matinit(motscorethr2);

   // DEFINING CONSERVED PWM AS REFERENCE
   refmot.matfreq=mattofreq(refmot.matprec);
   cout << "# of conserved instances: " <<   refmot.seqs.size() << endl;
   
   unsigned int nshuffle=5;
   unsigned int nshufflecons=15;

   int iter(0);
   // TEST ALLREF CONVERGENCE
   Sequence initRefInstances=refmot.refinstances;
   // /*
   cout << "Testing all sites convergence on reference species..." << endl;
   for (unsigned int k=0;k<nshuffle;k++)
   {
      // RANDOM SHUFFLING
      random_shuffle(initRefInstances.iseqs.begin(),initRefInstances.iseqs.end());
      refinstances.iseqs.clear();
      refinstances.seqs.clear();
      for (ivvint iv=initRefInstances.iseqs.begin();iv!=initRefInstances.iseqs.begin()+refmot.seqs.size();iv++){
         refinstances.iseqs.push_back(*iv);
         refinstances.seqs.push_back(vinttostring(*iv));

         // 0.1 step in log scale
         int numinst=refinstances.seqs.size(); 
         if (numinst%10!=0 && numinst>1) continue;

         cout << "# of sites for AllRef: " <<   refinstances.iseqs.size() << 
            " (iter=" << k+1 << "/" << nshuffle << ")" << endl;
         start=clock();
         vd dum(4,0.0);
         matrice=vvd(width,dum);
         countbases(*this,refinstances);
         countfreq(matrice);
         matfreq=matrice;
         freqtolog(matrice);
         matprec=matrice;
         finish=clock();
         dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
         outf << "AllRef " << numinst << " " << dif << endl;
      }
   }
   //*/

   // TEST HALPERN CONVERGENCE
   cout << "Testing Halpern convergence..." << endl;
   evolutionary_model=2;
   inittreedist();
   for (unsigned int k=0;k<nshufflecons;k++){

      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      vma tmpseqs;

      for (int i=0;i<refmot.seqs.size();i++){
         tmpseqs.push_back(refmot.seqs[i]);

         int numinst=tmpseqs.size(); 
         if (numinst%10!=0 && numinst>1) continue;

         cout << "# of sites for Halpern: " <<   tmpseqs.size() <<  " (iter=" << k+1 << "/" << nshufflecons << ")" << endl;
         seqs=tmpseqs;
         // MAXIMUM LIKELIHOOD
         start=clock();
         compprec();
         finish=clock();
         matfreq=mattofreq(matprec);
         dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
         outf << "Halpern_colopti " << numinst << " " << dif << endl;

         // MEAN POSTERIOR
         start=clock();
         compprec_MCMC();
         matfreq=mattofreq(matprec);
         finish=clock();
         matfreq=mattofreq(matprec);
         dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
         outf << "Halpern_colmean " << numinst << " " << dif << endl;
      }
   }

   // TEST FELSEN CONVERGENCE
   //   /*
   cout << "Testing Felsen convergence..." << endl;
   evolutionary_model=1;
   inittreedist();
   for (unsigned int k=0;k<nshufflecons;k++){
      // RANDOM SHUFFLING
      random_shuffle(refmot.seqs.begin(),refmot.seqs.end());
      vma tmpseqs;

      for (int i=0;i<refmot.seqs.size();i++){
         tmpseqs.push_back(refmot.seqs[i]);

         int numinst=tmpseqs.size(); 
         if (numinst%10!=0 && numinst>1) continue;

         cout << "# of sites for Felsen: " <<   tmpseqs.size() <<  " (iter=" << k+1 << "/" << nshufflecons << ")" << endl;
         seqs=tmpseqs;
         // MAXIMUM LIKELIHOOD
         start=clock();
         compprec();
         finish=clock();
         matfreq=mattofreq(matprec);
         dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
         outf << "Felsen_colopti " << numinst << " " << dif << endl;

         // MEAN POSTERIOR
         start=clock();
         compprec_MCMC();
         matfreq=mattofreq(matprec);
         finish=clock();
         matfreq=mattofreq(matprec);
         dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
         outf << "Felsen_colmean " << numinst << " " << dif << endl;
      }
   }
   // */

   return;
}

Motif comprefmot(Motif & motinit)
{

   Motif mot;
   mot=motinit;
   vd dum(4,0.0);
   mot.matrice=vvd(width,dum);
   countbases(mot,motinit.refinstances);
   //   displaymat(mot.matrice);
   //   cout << endl;
   countfreq(mot.matrice);
   mot.matfreq=mot.matrice;
   //   displaymat(mot.matfreq);
   //   getchar();
   freqtolog(mot.matrice);
   mot.matprec=mot.matrice;
   //Consensus sequence
   mot.bsinit="";
   for (int i=0;i<width;i++){
      string letter="A";
      double max(mot.matprec[i][0]);
      if (mot.matprec[i][1]>max) letter="T";
      else if (mot.matprec[i][2]>max) letter="C";
      else if (mot.matprec[i][3]>max) letter="G";
      mot.bsinit.append(letter);
   }

   //cout << mot.matprec << endl;
   mot.pvaluecomp();

   return mot;

}

// !! Initialize instances first with comprefinstances
   void 
Motif::comprefmot()
{
   vd dum(4,0.0);
   matrice=vvd(width,dum);
   countbases(*this,refinstances);
   countfreq(matrice);
   matfreq=matrice;
   freqtolog(matrice);
   matprec=matrice;
   //Consensus sequence
   bsinit="";
   for (int i=0;i<width;i++){
      string letter="A";
      double max(matprec[i][0]);
      if (matprec[i][1]>max) letter="T";
      else if (matprec[i][2]>max) letter="C";
      else if (matprec[i][3]>max) letter="G";
      bsinit.append(letter);
   }

   pvaluecomp();

   return;
}

   void
Motif::comprefinstancescons(unsigned int nspe)
{
   matinit(motscorethr2);
   refinstances.seqs.clear();
   refinstances.iseqs.clear();
   for (ivma ivm=seqs.begin();ivm!=seqs.end();ivm++){
      if (ivm->matches[nspe]){
         refinstances.iseqs.push_back(ivm->alignseq[nspe]);
         refinstances.seqs.push_back(vinttostring(ivm->alignseq[nspe]));
      }
   }
   return;
}

   void
Motif::comprefinstances(vseq & regs,unsigned int nspe)
{
   findinstances(regs);
   refinstances.iseqs.clear();
   refinstances.seqs.clear();
   for (ivseq ivs=regs.begin();ivs!=regs.end();ivs++){
      ivs->nmot=0;
      for (ivinstseq ivi=ivs->instances.begin();ivi!=ivs->instances.end();ivi++){
         if (ivi->species==nspe){
            ivs->nmot++;
            vint iv=ivi->iseq;
            if (ivi->sens==-1) iv=reversecomp(ivi->iseq);
            refinstances.iseqs.push_back(iv);
            refinstances.seqs.push_back(vinttostring(iv));
         }
      }
   }
   nbmot=refinstances.iseqs.size(); // NON conserved here
   return;

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
getmatrices(ifstream & file, Motif & mot)
{
   string dum;
   cout << "\n mat \n" << endl;
   getline(file,dum);
   while(!file.eof()){
      if (dum.find("rawmatfreq")!=string::npos){
         getline(file,dum);
         int element;
         for (int i=0;i<4;i++){
            int j=0;
            stringstream line(dum);
            int nfreq;
            while (j!=width) {
               line>>nfreq;
               mot.matrice[j][i]+=nfreq;
               cout << nfreq << " ";
               j++;
            }
            cout << "\n";
            getline(file,dum);
         }
      } else if (dum.find("realmatfreq")!=string::npos){
         getline(file,dum);
         int ntot;
         stringstream tot;
         tot.str(dum); 
         tot >> ntot;
         getline(file,dum);
         for (int i=0;i<4;i++){
            int j=0;
            stringstream line(dum);
            double nfreq;
            while (j!=width) {
               line>>nfreq;
               mot.matrice[j][i]+=floor(nfreq*((double)ntot/100)+0.5);//round to nearest value
               cout << floor(nfreq*((double)ntot/100)+0.5) << " "; 
               j++;
            }
            cout << "\n";
            getline(file,dum);
         }
      } 
      else
      {
         cout << "Bad typography! Realmatfreq or rawmatfreq?" << endl;
         exit(1);
      }
      cout << "\n";
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
Motif::calclambdaposneg(vseq & vscore)
{
   unsigned int nbneg=0;
   unsigned int nbpos=0;
   unsigned int nbtbpos=0;
   unsigned int nbtbneg=0;
   for (ivseq iseq=vscore.begin();iseq!=vscore.end();iseq++){
      if (iseq->sign==1){
         nbtbpos+=(*iseq).nbtb;
         int ncons;
         if (motwidth>nbtbpos){
            ncons=0;
         }
         else ncons=nbmatchcons(*iseq);
         nbpos+=ncons;
         nbmot+=ncons;
      } else if (iseq->sign==-1){
         nbtbneg+=(*iseq).nbtb;
         if (motwidth<=nbtbneg){
            nbneg+=nbmatchcons(*iseq);
         }
      }
   }
   lambdatrain=(double)nbpos/nbtbpos;
   lambda=(double)nbneg/nbtbneg;
   if (lambda!=0)  pvalue=gsl_ran_poisson_pdf(nbpos,lambda*nbtbpos);
   else pvalue=0;
}

   void
Motif::calclambdaback()
{
   unsigned int nbback=0;
   unsigned int tottest=0;
   for (ivseq iseq=regtests.begin();iseq!=regtests.end();iseq++){
      tottest+=(*iseq).nbtb;
      nbback+=nbmatchcons(*iseq);
   }
   lambda=nbback/(double)tottest;

   if (lambda==0.) lambda=1.e-10;
}

   void
Motif::calclambda()
{
   unsigned int nbback=0;
   unsigned int tottest=0;
   for (ivseq iseq=regtests.begin();iseq!=regtests.end();iseq++){
      tottest+=(*iseq).nbtb;
      nbback+=nbmatchcons(*iseq);
   }
   lambda=nbback/(double)tottest;

   unsigned int nbints=0;
   unsigned int totints=0;
   for (ivseq iseq=regints.begin();iseq!=regints.end();iseq++){
      nbints+=(*iseq).nbtb;
      totints+=nbmatchcons(*iseq);
   }
   lambdatrain=totints/(double)nbints;

   if (lambda==0.) lambda=1.e-10;
   else if (lambdatrain==0.) lambdatrain=1.e-10;
}

   void
Motif::lambdacomp()
{
   unsigned int nbbacktemp=0;
   unsigned int nbback=0;
   unsigned int tottest=0;
   vseq::const_iterator iseq;
   for (iseq=regtests.begin();iseq!=regtests.end();iseq++){
      tottest+=(*iseq).nbtb;
      nbbacktemp=nbmatchmat(*iseq);
      if (nbbacktemp<distwidth) distmot[nbbacktemp]++;
      nbback+=nbbacktemp;
      // cout << (*iseq).nbtb << " " << tottest << " " << nbbacktemp << " " << nbback << endl;
   }
   nbmatchback=nbback;
   lambda=nbback/(double)tottest;
}

   void
Motif::lambdacompcons()
{
   unsigned int nbbacktemp=0;
   unsigned int nbback=0;
   unsigned int tottest=0;
   vseq::iterator iseq;
   for (iseq=regtests.begin();iseq!=regtests.end();iseq++){
      tottest+=(*iseq).nbtb;
      if (motwidth>tottest){
         nbbacktemp=0;
      }
      else nbbacktemp=nbmatchcons(*iseq);
      if (nbbacktemp<distwidth) distmot[nbbacktemp]++;
      nbback+=nbbacktemp;
      // cout << (*iseq).nbtb << " " << tottest << " " << nbbacktemp << " " << nbback << endl;
   }
   nbmatchback=nbback;
   lambda=nbback/(double)tottest;
}

   void
Motif::pvaluecomp()
{
   lambdacomp();

   calcscorepoiss();
   vseq::iterator iseq;
   pvalue=0.0;


   int nbtot=0;
   int nbbtrain=0;

   //   for (iseq=regints.begin();iseq!=regints.begin()+nbvalidated;iseq++){
   for (iseq=regints.begin();iseq!=regints.end();iseq++){
      pvalue+=log(gsl_ran_poisson_pdf((*iseq).nmot,lambda*(*iseq).nbtb));
      //	   cout << (*iseq).nmot << " " << lambda << " " << (*iseq).nbtb << " " << pvalue << endl;
      nbbtrain+=(*iseq).nbtb;
      nbtot+=(*iseq).nmot;
   }
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
Motif::displaywname(ostream & streamfile)
{
   streamfile << name << "\t";
   streamfile << bsinit << "\t";
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

   int
Motif::nbmatchmat (const Sequence & seq)
{
   int nmat=0;
   unsigned int i=0;
   unsigned int len=seq.iseqs[0].size();
   for (vint::const_iterator istr=seq.iseqs[0].begin();istr!=seq.iseqs[0].end()-motwidth+1;istr++){
      double score=scoref(istr,matprec);
      if (score>motscorethr2){
         nmat++;
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
         if (score>motscorethr2){
            nmat++;
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
   return nmat;
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

// draws a site with matfreq probs    
   vint 
Motif::drawsite(double scorethr)
{
   vint site(motwidth,4);
   while (scoref(site,matprec)<scorethr){
      for (unsigned int i=0;i<motwidth;i++){
         vd w=matfreq[i];
         double pdist[4]={w[0],w[1],w[2],w[3]};
         gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
         unsigned int base=gsl_ran_discrete (gslran,gsldist);
         site[i]=base;
      }
   }
   return site;
}

// compare binding sites by their affinity
// sorts with highests scoring TFBSs first
bool operator<(const TFBS & bs1,const TFBS & bs2)
{
   //return bs1.score > bs2.score;
   return bs1.prob > bs2.prob;
}

bool operator<(const Instance & inst1,const Instance & inst2)
{
   if (inst1.chrom != inst2.chrom) return inst1.chrom < inst2.chrom;
   else return inst1.coord < inst2.coord;
}

   ostream &
operator <<(ostream &os,const TFBS & bs)
{
   vint site=bs.site;
   os << vinttostring(site) << "\t";
   //os.precision(4);
   os << bs.score;
   //   os.precision(2);
   //   os << bs.prob << "\t";
   //   os.precision(6);
   //   os << floor(bs.num+0.5);
   return os;
}

   ostream &
operator <<(ostream &os,const vtfbs & vbs)
{
   for (civtfbs civ=vbs.begin();civ!=vbs.end();civ++){
      os << *civ << "\n";
   }
   return os;
}

bool operator<(const GroupInstance & ginst1,const GroupInstance & ginst2)
{
   return ginst1.score > ginst2.score;
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
   unsigned int dist=100000;
   for (ivTSS ivt=TSSs.begin();ivt!=TSSs.end();ivt++){
      if (abs((int)((*ivt).coord)-(int)start+(int)scanwidth/2)<(int)dist){
         dist=abs((int)(*ivt).coord-(int)start+(int)scanwidth/2);
         besttss=*ivt;
         //cout << besttss.gene << endl;
      }
   }
}

   unsigned int
nbmatches(vinst vin, int start)
{
   unsigned nbmat=0;
   for (ivinst ivin=vin.begin();ivin!=vin.end();ivin++){
      if ((*ivin).coord>=start && (*ivin).coord<start+scanwidth+1) nbmat++;
      if ((*ivin).coord>start+scanwidth) return nbmat;
   }
}

   int
Motif::nbmatchnmask (Sequence & seq,unsigned int moti)
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
            for (unsigned int j=i;j<i+width;j++){
               seq.iseqs[0][j]=4;
            }
            ma.mask();
            instances[seq.chrom].push_back(Instance(seq.chrom,seq.start+i,1,moti));
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
               for (unsigned int j=i;j<i+width;j++){
                  seq.iseqs[0][j]=4;
               }
               ma.mask();
               instances[seq.chrom].push_back(Instance(seq.chrom,seq.start+i,-1,moti));
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

   int
Motif::nbmatchwomask (Sequence & seq,unsigned int moti)
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
            instances[seq.chrom].push_back(Instance(seq.chrom,seq.start+i,1,moti));
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
               instances[seq.chrom].push_back(Instance(seq.chrom,seq.start+i,-1,moti));
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

   int
Motif::nbmatchconsbest1kb (Sequence & seq)
{
   unsigned int i=0;
   unsigned int len=seq.iseqs[0].size();
   vd vpos;
   for (civint istr=seq.iseqs[0].begin();istr!=seq.iseqs[0].end()-width+1;istr++){
      double score=scoref(istr,matprec);
      if (score>motscorethr2){
         Motalign ma(i,seq,*this, 1);
         if (ma.iscons()){
            vpos.push_back(i);
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
               vpos.push_back(i);
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

   unsigned int nmat=0;
   if (len>scanwidth){
      for (int j=0;j<len-scanwidth;j++){
         unsigned int counter=0;
         for (ivd iv=vpos.begin();iv!=vpos.end();iv++){
            if (*iv>=j && *iv<j+scanwidth){
               counter++;
            }
         }
         if (counter>nmat) nmat=counter;
      }
   }   
   else{
      nmat=vpos.size();
   }

   return nmat;
}

   int
Motif::nbmatchconsnmask (Sequence & seq)
{
   int nmat=0;
   unsigned int i=0;
   unsigned int len=seq.iseqs[0].size();
   for (civint istr=seq.iseqs[0].begin();istr!=seq.iseqs[0].end()-width+1;istr++){
      double score=scoref(istr,matprec);

      if (score>motscorethr2){
         Motalign ma(i,seq,*this, 1);
         if (ma.iscons()){
            for (unsigned int j=i;j<i+width;j++){
               seq.iseqs[0][j]=4;
            }
            ma.mask();
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
               for (unsigned int j=i;j<i+width;j++){
                  seq.iseqs[0][j]=4;
               }
               ma.mask();
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

   double
Motif::scorematchcons (Sequence & seq)
{
   double scr=0;
   unsigned int i=0;
   unsigned int len=seq.iseqs[0].size();
   for (civint istr=seq.iseqs[0].begin();istr!=seq.iseqs[0].end()-width+1;istr++){
      double score=scoref(istr,matprec);

      if (score>motscorethr2){
         Motalign ma(i,seq,*this, 1);
         if (ma.iscons()){
            scr+=score;
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
               scr+=score;
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
   return scr;
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


   int
Motif::nbmatchnmaskforsvg (Sequence & seq,unsigned int moti)
{
   int nmat=0;
   for (unsigned int spe=0;spe<nbspecies;spe++){
      if (seq.species[spe]){
         unsigned int len=seq.iseqs[spe].size();
         unsigned int i=0;
         for (civint istr=seq.iseqs[spe].begin();istr!=seq.iseqs[spe].end()-width+1;istr++){
            double score=scoref(istr,matprec);
            if (score>7){
               for (unsigned int j=i;j<i+width;j++){
                  seq.iseqs[spe][j]=4;
               }
               Instanceseq inst(moti,1,seq.imapsinv[spe][i],i,score,spe,name);
               seq.instances.push_back(inst);
               nmat++;
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
               if (score>7){
                  for (unsigned int j=i;j<i+width;j++){
                     seq.iseqs[spe][j]=4;
                  }
                  Instanceseq inst(moti,1,seq.imapsinv[spe][i],i,score,spe,name);
                  seq.instances.push_back(inst);
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
            i++;
         }
      }
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
Motif::matinithamming(double scth,unsigned int numhamm)
{
   seqs.clear();
   for (ivseq iseq=regints.begin();iseq!=regints.end();iseq++){
      (*iseq).nmot=0;
   }
   //   cout << "begin init\n";
   for (ivseq iseq=regints.begin();iseq!=regints.end();iseq++){
      Sequence & seq=*iseq;
      unsigned int len=seq.iseqs[0].size();
      unsigned int i=0;
      for (vint::const_iterator istr=seq.iseqs[0].begin();istr!=seq.iseqs[0].end()-motwidth+1;istr++){
         if (scoref(istr,matprec)>scth){
            vint::const_iterator endci=seq.iseqs[0].end()-motwidth+1;
            unsigned int sh=shift(istr,matprec,endci,motwidth);
            Motalign ma=mahamming(i+sh,seq,*this, 1, numhamm);
            (*iseq).nmot++;
            if (ma.iscons()){
               seqs.push_back(ma);
               if (i+sh<len-2*motwidth){
                  istr+=motwidth-1+sh;
                  i+=motwidth-1+sh;
               }
               else {
                  break;
               }
            }
         }
         else if (scoref(istr,matprecrevcomp)>scth){ 
            vint::const_iterator endci=seq.iseqs[0].end()-motwidth+1;
            unsigned int sh=shift(istr,matprecrevcomp,endci,motwidth);
            Motalign ma=mahamming(i+sh,seq,*this,-1, numhamm);
            (*iseq).nmot++;
            if (ma.iscons()){
               seqs.push_back(ma);
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
   //   cout << "end init\n";
   nbmot=seqs.size();
}

   bool
Motalign::iscons()
{
   if (species==1){
      int nbfr=0;
      if (matches[5]) nbfr++; 
      if (matches[6] || matches[7]) nbfr++;
      if (matches[8]) nbfr++;
      if (matches[9] || matches[10] || matches[11]) nbfr++;
      if (nbfr>1) return true;
   }
   else if (species==2){
      int nbfr=0;
      //     if (matches[0] || matches[1]) nbfr++;
      if (matches[2] || matches[3] || matches[4] || matches[5]) nbfr++;
      if (matches[6] || matches[7]) nbfr++;
      if (matches[8]) nbfr++;
      if (matches[9]) nbfr++;
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

   void
Motif::compprec_test()
{
   time_t tstart,tend;
   clock_t start,finish;
   double dif;
   vd col1,col2;
   vvd matopti,matmean,matoptimc;
   //for (ivma ivm=seqs.begin();ivm!=seqs.end();ivm++) ivm->print();
   cout << "There are " << seqs.size() << " conserved instances" << endl;

   for (unsigned int i=0;i<width;i++){
      vd col;
      //i=3;
      cout << "pos " << i+1 << endl;
      // colmean_test_likelihood_w_dirichlet(i,this);
      //colmean_testCV_w_dirichlet(i,this);
      //colmean_test_autocorr_w_dirichlet(i,this);
      //colmean_test_window_for_mean(i,this);
      //exit(9);
      //col=colmean_test_rejection_w_dirichlet_metropolis(i,this);
      //      col=colmean_test(i,this);
      //         col=colmean_testCV(i,this);
      //      time (&tstart);
      //      col1=colmean_test_numiter(i,this);
      //      time (&tend);
      //      dif = difftime (end,tstart);
      //      cout << "Time for random walk: " << dif << " seconds" << endl;
      //      time (&tstart);
      //      time (&tend);
      //      dif = difftime (tend,tstart);
      //      cout << "Time for Dirichlet: " << dif << " seconds " <<  endl;
      //      double distKL=0.;
      //      for (int ib=0;ib<4;ib++){
      //         distKL+=col2[ib]*log(col2[ib]/col1[ib]);
      //      }
      //      cout << "KL distance d(2,1): " << distKL <<  endl;
      //      distKL=0;
      //      for (int ib=0;ib<4;ib++){
      //         distKL+=col1[ib]*log(col1[ib]/col2[ib]);
      //      }
      //      cout << "KL distance d(1,2): " << distKL <<  endl;
      //      col=colmean_test_rejection_w_dirichlet_metropolis(i,this);
      //      col=colmean_test_rejection_w_dirichlet(i,this);
      //  col=colmean_test_rejection_n_autocorr(i,this);
      //               col=colmean_t_test(i,this);
      //               continue;
      //           col=colmean_testCV_new(i,this);
      //      exit(9);

      //      vvd dcol;
      //      time (&tstart);
      //      dcol=colmean(i,this);
      //      time (&tend);
      //cout << "MEAN WITH WEIGHT " << col << endl;
      //      cout << "MEAN (MCMC) " << dcol[0];
      //      printf ("---> time: %.6lf seconds.\n", dif );
      //matmean.push_back(dcol[0]);
      //cout << "MAX (MC) " << dcol[1];
      //matoptimc.push_back(dcol[1]);

      //   time (&tstart);
      start=clock();
      col=colmean(i,this);
      finish=clock();
      dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
      cout << "MEAN (Dirichlet) " << col;
      printf ("---> time: %.6lf milliseconds.\n", dif );
      matmean.push_back(col);
      start=clock();
      col=colopti(i,this);
      finish=clock();
      dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
      cout << "MAX (SD) " << col;
      printf ("---> time: %.6lf milliseconds.\n", dif );
      matopti.push_back(col);
   }

   ofstream outf("motmeldb_test.txt");
   // INIT
   display(outf);
   // MEAN
   matprec=matmean;
   display(outf);
   // MAX
   matprec=matopti;
   display(outf);
   outf.close();
   exit(3);
   matprecrevcomp=reversecomp(matprec);
}

   vvd 
mattoenergy(vvd & mat)
{
   vd dumd(4,0.0);
   vvd matenergy(width,dumd);
   int i=0;
   for (ivvd iv=mat.begin();iv!=mat.end();iv++){
      double max=-100.;
      for (ivd col=(*iv).begin();col!=(*iv).end();col++){
         if (*col>max) max=*col;
      }
      int j=0;
      for (ivd col=(*iv).begin();col!=(*iv).end();col++){
         matenergy[i][j]=max-*col;
         j++;
      }
      i++;
   }
   return matenergy;

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

// for different window sizes, compute the kl distance between means
   vd
colmean_test_window_for_mean(unsigned int pos,Motif * mot)
{
   //   clock_t start,finish;
   //   double dif;
   //   start=clock();
   unsigned int Ndir(0);
   //cout << mot->seqs.size() << endl;
   for (ivma ivm=mot->seqs.begin();ivm!=mot->seqs.end();ivm++){
      Ndir+=accumulate(ivm->matches.begin(),ivm->matches.end(),0)/2; 
   }
   //cout << "Ndir=" << Ndir << endl;

   unsigned int nsample(30);
   //cout << "nsample=" << nsample << endl;

   //unsigned int nmean(100);
   //cout << "nmean=" << nmean << endl;
   double dkl=1.;
   double dcutoff(1e-5);
   //cout << "dcutoff=" << dcutoff << endl;
   //
   ofstream outf("testwindow.dat");

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

   vd vdum(4,0.);
   vd wprev(4,0.),wafter(4,0.),wmean(4,0.);

   for (unsigned int ws=1;ws<100;ws+=2){


      vvd wwindow(ws,vdum);

      // Starting point
      double alpinit[4]={alpha+Ndir*winit[0],alpha+Ndir*winit[1],beta+Ndir*winit[2],beta+Ndir*winit[3]};
      double thetainit[]={0.,0.,0.,0.};
      gsl_ran_dirichlet(gslran,4,alpinit,thetainit);
      wprev[0]=thetainit[0];
      wprev[1]=thetainit[1];
      wprev[2]=thetainit[2];
      wprev[3]=1-thetainit[0]-thetainit[1]-thetainit[2];

      wmean=wprev;
      fprev=loglikelyhood(wprev,par);

      // =========
      // MAIN LOOP
      // =========

      // counter counts accepted trials 
      unsigned int counter(1),iterok(1),iter(0);
      double logalph,prob; // for Metropolis Hastings
      // BURN-IN
      unsigned int bicount(0),bicutoff(nsample);

      // we average on 10 experiments
      while (iterok<=10*nsample*ws) {

         //double alp[4]={alpha+Ndir*wprev[0],alpha+Ndir*wprev[1],beta+Ndir*wprev[2],beta+Ndir*wprev[3]};
         double alp[4]={alpha+Ndir*winit[0],alpha+Ndir*winit[1],beta+Ndir*winit[2],beta+Ndir*winit[3]};
         double theta[]={0.,0.,0.,0.};
         gsl_ran_dirichlet(gslran,4,alp,theta);
         wafter[0]=theta[0];
         wafter[1]=theta[1];
         wafter[2]=theta[2];
         wafter[3]=1-theta[0]-theta[1]-theta[2];

         fafter=loglikelyhood(wafter,par);

         // Metropolis-Hastings
         // Move to new position with probability alpha, else reject
         double theta2[]={wprev[0],wprev[1],wprev[2],wprev[3]};
         double alp2[4]={alpha+Ndir*wafter[0],alpha+Ndir*wafter[1],beta+Ndir*wafter[2],beta+Ndir*wafter[3]};
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

               wwindow[counter-1]=wafter;
               if (counter%ws==0){
                  vd wprevmean=wmean;
                  wmean=vdum;
                  for (unsigned int iw=0;iw<ws;iw++){
                     for (unsigned int ib=0;ib<4;ib++){
                        wmean[ib]+=wwindow[iw][ib]/ws;
                     }
                  }

                  dkl=0.;
                  for (int ib=0;ib<4;ib++){
                     dkl+=wmean[ib]*log(wmean[ib]/wprevmean[ib]);
                  }

                  cout << ws << " " << dkl << endl;
                  outf << ws << "\t" << dkl << endl;
                  counter=1;
               }
               else counter++;
            }
            cout.flush();
            iterok++;
         }
         iter++;
      } 

   }

   outf.close();
   exit(2);

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

   vvd
colmean_RW(unsigned int pos,Motif * mot)
{
   vd dum(4,0.);
   double Rmax(.013),rad,theta,phi; // for random walk
   double fprev,fafter;
   unsigned int nsample(20); // vary nsample for auto-correlation test, nmean is the size of the sliding window
   void * par[2]={(void *)mot,&pos};

   unsigned int nmean(150); // vary nsample for auto-correlation test, nmean is the size of the sliding window
   double dcutoff(0.01); // loop ending requirement
   unsigned int nforcvmean=20;

   vd winit(4,0.);
   winit[0]=conca*exp(mot->matprec[pos][0]);
   winit[1]=conct*exp(mot->matprec[pos][1]);
   winit[2]=concc*exp(mot->matprec[pos][2]);
   winit[3]=concg*exp(mot->matprec[pos][3]);
   vd wfinal(4,0.);

   vd wbestfinal(4,0.);

   // DIFFERENT PATHS FOR AVERAGING
   for (unsigned int imean=0;imean<nforcvmean;imean++){

      vd f(nmean,0.); // for f mean computation, taken every nsample steps
      vvd w(nmean,dum); // for w mean computation, taken every nsample steps
      vint nrej(nmean,0); // for rejection rate, taken every step
      double distKL(10.); // loop ending requirement
      vd wprev(4,0.),wafter(4,0.),eps(4,0.),wmean(4,0.),wprevmean(4,0.);
      vd wbest(4,0.);
      double fbest(0.);

      // Starting point
      wprev=winit;

      fprev=likelyhood(wprev,par);
      f[0]=fprev;
      w[0]=wprev;
      wprevmean=wprev;
      nrej[0]=1;

      // =========
      // MAIN LOOP
      // =========

      double alph,prob; // for Metropolis Hastings
      // counter counts accepted trials 
      // tc counts accepted trials every nsample steps modulo nmean
      // iter counts every step
      unsigned int counter(1),tc(0),iter(1);

      while (abs(distKL)>dcutoff) {

         // ===========
         // Random walk
         // ===========
         //  PROPOSAL = UNIFORM IN A BALL
         rad=Rmax*gsl_rng_uniform (gslran);
         theta=M_PI*gsl_rng_uniform (gslran);
         phi=2*M_PI*gsl_rng_uniform (gslran);
         eps[0]=rad*sin(theta)*cos(phi);
         eps[1]=rad*sin(theta)*sin(phi);
         eps[2]=rad*cos(theta);
         eps[3]=-(eps[0]+eps[1]+eps[2]);

         wafter=wprev+eps;

         // check if new probabilities are consistent
         int wtest(0);
         for (int i=0;i<4;i++) if (wafter[i]<0 || wafter[i]>1) wtest=1;
         if (wtest) continue;


         fafter=likelyhood(wafter,par);

         // Metropolis-Hastings algorithm
         alph=min(1.,fafter/fprev);
         prob=gsl_rng_uniform (gslran);

         if (prob<=alph){
            // we define nsample to avoid auto-correlation
            // we want one point every nsample, nmean times
            tc=counter%(nmean*nsample);
            if (tc%nsample==0){
               tc/=nsample;
               unsigned int prevcounter;
               if (tc==0) prevcounter=nmean-1;
               else prevcounter=tc-1;
               f[tc]=fafter;
               w[tc]=wafter;

               if (fafter>fbest){
                  fbest=fafter;
                  wbest=wafter;
               }

               vd fnum(4,0.); 
               for (int itc=0;itc<nmean;itc++){
                  for (int ib=0;ib<3;ib++){
                     fnum[ib]+=w[itc][ib]/nmean;
                  }  
               }

               if (tc==0){
                  distKL=0.;
                  for (int ib=0;ib<3;ib++){
                     wmean[ib]=fnum[ib];
                     distKL+=wmean[ib]*log(wmean[ib]/wprevmean[ib]);
                  }
                  wmean[3]=1.-wmean[0]-wmean[1]-wmean[2];
                  distKL+=wmean[3]*log(wmean[3]/wprevmean[3]);
                  wprevmean=wmean;
               }
            }
            counter++;
            wprev=wafter;
            fprev=fafter;
         }

         iter++;
      }

      for (unsigned int ib=0;ib<4;ib++){
         wfinal[ib]+=wmean[ib]/nforcvmean;
         wbestfinal[ib]+=wbest[ib]/nforcvmean;
      }
   }

   vd resmean;
   resmean.push_back(log(wfinal[0]/conca));
   resmean.push_back(log(wfinal[1]/conct));
   resmean.push_back(log(wfinal[2]/concc));
   resmean.push_back(log(wfinal[3]/concg));
   vd resmax;
   resmax.push_back(log(wbestfinal[0]/conca));
   resmax.push_back(log(wbestfinal[1]/conct));
   resmax.push_back(log(wbestfinal[2]/concc));
   resmax.push_back(log(wbestfinal[3]/concg));

   vvd coltot;
   coltot.push_back(resmean);
   coltot.push_back(resmax);

   return coltot;
}

   vd
colmean_dirichlet(unsigned int pos,Motif * mot)
{
   unsigned int Ndir(0);
   for (ivma ivm=mot->seqs.begin();ivm!=mot->seqs.end();ivm++){
      Ndir+=accumulate(ivm->matches.begin(),ivm->matches.end(),0)/2; 

   }
   cout << "Ndir=" << Ndir << endl;

   unsigned int nsample(20);
   cout << "nsample=" << nsample << endl;

   //unsigned int nmean(100);
   //cout << "nmean=" << nmean << endl;
   double dkl=1.;
   double dcutoff(1e-4);
   cout << "dcutoff=" << dcutoff << endl;

   // INIT
   vd winit(4,0.);
   winit[0]=conca*exp(mot->matprec[pos][0]);
   winit[1]=conct*exp(mot->matprec[pos][1]);
   winit[2]=concc*exp(mot->matprec[pos][2]);
   winit[3]=1.-winit[0]-winit[1]-winit[2];

   double fprev,fafter;
   void * par[2]={(void *)mot,&pos};

   vd wprev(4,0.),wafter(4,0.),wmean(4,0.);

   // Starting point
   double alpinit[4]={alpha+Ndir*winit[0],alpha+Ndir*winit[1],beta+Ndir*winit[2],beta+Ndir*winit[3]};
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
   unsigned int counter(1),iterok(1);
   double logalph,prob; // for Metropolis Hastings
   // BURN-IN
   unsigned int bicount(0),bicutoff(nsample);

   //while (counter<=nmean) {
   while (dkl>=dcutoff) {

      double alp[4]={alpha+Ndir*wprev[0],alpha+Ndir*wprev[1],beta+Ndir*wprev[2],beta+Ndir*wprev[3]};
      double theta[]={0.,0.,0.,0.};
      gsl_ran_dirichlet(gslran,4,alp,theta);
      wafter[0]=theta[0];
      wafter[1]=theta[1];
      wafter[2]=theta[2];
      wafter[3]=1-theta[0]-theta[1]-theta[2];

      fafter=loglikelyhood(wafter,par);

      // Metropolis-Hastings
      // Move to new position with probability alpha, else reject
      double theta2[]={wprev[0],wprev[1],wprev[2],wprev[3]};
      double alp2[4]={alpha+Ndir*wafter[0],alpha+Ndir*wafter[1],beta+Ndir*wafter[2],beta+Ndir*wafter[3]};
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

            vd wprevmean=wmean;
            for (unsigned int ib=0;ib<4;ib++){
               //wmean[ib]+=wafter[ib]/nmean;
               wmean[ib]=(counter-1.)/counter*wprevmean[ib]+wafter[ib]/counter;
            }

            dkl=0.;
            for (int ib=0;ib<4;ib++){
               dkl+=wmean[ib]*log(wmean[ib]/wprevmean[ib]);
            }

            counter++;
         }

         iterok++;
      }
   } 

   cout << "counter=" << counter << endl;
   vd res;
   res.push_back(log(wmean[0]/conca));
   res.push_back(log(wmean[1]/conct));
   res.push_back(log(wmean[2]/concc));
   res.push_back(log(wmean[3]/concg));

   return res;
}

   vd
colmean_dirichlet_metropolis(unsigned int pos,Motif * mot)
{
   unsigned int Ndir(30);
   vd dum(4,0.);
   vd winit(4,0.);
   winit[0]=conca*exp(mot->matprec[pos][0]);
   winit[1]=conct*exp(mot->matprec[pos][1]);
   winit[2]=concc*exp(mot->matprec[pos][2]);
   winit[3]=concg*exp(mot->matprec[pos][3]);
   double fprev,fafter;
   void * par[2]={(void *)mot,&pos};

   unsigned int nmean(80); // vary nsample for auto-correlation test, nmean is the size of the sliding window

   vd wfinal(4,0.);
   unsigned int numtotmean(5000);
   for (unsigned int totmean=0;totmean<numtotmean;totmean++){

      vd wprev(4,0.),wafter(4,0.);
      // Starting point
      wprev=winit;

      fprev=likelyhood(wprev,par);
      // =========
      // MAIN LOOP
      // =========

      double alph,prob; // for Metropolis Hastings
      // counter counts accepted trials 
      // tc counts accepted trials every nsample steps modulo nmean
      // iter counts every step
      unsigned int counter(2);
      while (counter<=nmean) {

         double theta[]={0.,0.,0.,0.};
         //double alp[4]={Ndir*wprev[0],Ndir*wprev[1],Ndir*wprev[2],Ndir*wprev[3]};
         //double alp[4]={Ndir*winit[0],Ndir*winit[1],Ndir*winit[2],Ndir*winit[3]};
         double alp[4]={alpha+Ndir*wprev[0],alpha+Ndir*wprev[1],beta+Ndir*wprev[2],beta+Ndir*wprev[3]};
         gsl_ran_dirichlet(gslran,4,alp,theta);
         wafter[0]=theta[0];
         wafter[1]=theta[1];
         wafter[2]=theta[2];
         wafter[3]=1-theta[0]-theta[1]-theta[2];

         fafter=likelyhood(wafter,par);

         // Metropolis-Hastings algorithm
         // Move to new position with probability alpha, else reject
         double theta2[]={wprev[0],wprev[1],wprev[2],wprev[3]};
         double alp2[4]={alpha+Ndir*wafter[0],alpha+Ndir*wafter[1],beta+Ndir*wafter[2],beta+Ndir*wafter[3]};
         alph=min(1.,fafter*gsl_ran_dirichlet_pdf(4,alp2,theta2)/(fprev*gsl_ran_dirichlet_pdf(4,alp,theta)));

         //alph=min(1.,fafter/fprev);
         prob=gsl_rng_uniform (gslran);

         if (prob<=alph){
            // we define nsample to avoid auto-correlation
            // we want one point every nsample, nmean times
            for (unsigned int ib=0;ib<4;ib++){
               wprev[ib]=(counter-1.)/counter*wprev[ib]+wafter[ib]/counter;
            }

            counter++;
            wprev=wafter;
            fprev=fafter;
         }
      }

      for (unsigned int ib=0;ib<4;ib++){
         wfinal[ib]+=wprev[ib]/numtotmean;
      }
   }

   vd res;
   res.push_back(log(wfinal[0]/conca));
   res.push_back(log(wfinal[1]/conct));
   res.push_back(log(wfinal[2]/concc));
   res.push_back(log(wfinal[3]/concg));

   return res;
}

// t-test between consecutive samples until convergence is reached
   vd
colmean_t_test(unsigned int pos,Motif * mot)
{
   vd dum(4,0.);
   double Rmax(.013),rad,theta,phi; // for random walk
   double fprev,fafter;
   unsigned int nsample(20); // vary nsample for auto-correlation test, nmean is the size of the sliding window
   void * par[2]={(void *)mot,&pos};

   double dcutoff(0.01); // loop ending requirement
   unsigned int nforcvmean=20;

   vd winit(4,0.);
   winit[0]=conca*exp(mot->matprec[pos][0]);
   winit[1]=conct*exp(mot->matprec[pos][1]);
   winit[2]=concc*exp(mot->matprec[pos][2]);
   winit[3]=concg*exp(mot->matprec[pos][3]);
   vd wfinal(4,0.);

   // DIFFERENT PATHS FOR AVERAGING
   for (unsigned int imean=1;imean<=nforcvmean;imean++){


      unsigned int nmean; // number of points for mean
      nmean=2*10*imean;

      cout << "List of 2x" << nmean/2 << " values." << endl;

      vd f(nmean,0.); // for f mean computation, taken every nsample steps
      vvd w(nmean,dum); // for w mean computation, taken every nsample steps
      double distKL(10.); // loop ending requirement
      vd wprev(4,0.),wafter(4,0.),eps(4,0.);
      vd wmean1(4,0.),wmean2(4,0.),wsd1(4,0.),wsd2(4,0.);
      vd samplevar(4,0.);

      // Starting point
      wprev=winit;

      fprev=likelyhood(wprev,par);
      f[0]=fprev;
      w[0]=wprev;

      // =========
      // MAIN LOOP
      // =========

      double alph,prob; // for Metropolis Hastings
      // counter counts accepted trials 
      // tc counts accepted trials every nsample steps modulo nmean
      // iter counts every step
      unsigned int counter(1),tc(0),iter(1);

      int gotonext(1);
      while (gotonext) {

         // ===========
         // Random walk
         // ===========
         //  PROPOSAL = UNIFORM IN A BALL
         rad=Rmax*gsl_rng_uniform (gslran);
         theta=M_PI*gsl_rng_uniform (gslran);
         phi=2*M_PI*gsl_rng_uniform (gslran);
         eps[0]=rad*sin(theta)*cos(phi);
         eps[1]=rad*sin(theta)*sin(phi);
         eps[2]=rad*cos(theta);
         eps[3]=-(eps[0]+eps[1]+eps[2]);

         wafter=wprev+eps;

         // check if new probabilities are consistent
         int wtest(0);
         for (int i=0;i<4;i++) if (wafter[i]<0 || wafter[i]>1) wtest=1;
         if (wtest) continue;


         fafter=likelyhood(wafter,par);

         // Metropolis-Hastings algorithm
         alph=min(1.,fafter/fprev);
         prob=gsl_rng_uniform (gslran);

         if (prob<=alph){
            // we define nsample to avoid auto-correlation
            // we want one point every nsample, nmean times
            tc=counter%(nmean*nsample);
            if (tc%nsample==0){
               tc/=nsample;
               unsigned int prevcounter;
               if (tc==0) prevcounter=nmean-1;
               else prevcounter=tc-1;
               f[tc]=fafter;
               w[tc]=wafter;

               if (tc==0){
                  // two-sided t-test, confidence 99%
                  double tcutoff,tval;
                  tcutoff=gsl_cdf_tdist_Pinv((0.99+1)/2.,nmean-2);

                  int isaccepted(1);
                  for (int ib=0;ib<3;ib++){
                     // !! Here we compute sample standard deviations
                     for (int itc=0;itc<nmean/2;itc++){
                        wmean1[ib]+=w[itc][ib];
                        wsd1[ib]+=pow(w[itc][ib],2);
                     }
                     wsd1[ib]-=pow(wmean1[ib],2)/(nmean/2);
                     wsd1[ib]/=(nmean/2-1);
                     wmean1[ib]/=(nmean/2);

                     for (int itc=nmean/2;itc<nmean;itc++){
                        wmean2[ib]+=w[itc][ib];
                        wsd2[ib]+=pow(w[itc][ib],2);
                     }
                     wsd2[ib]-=pow(wmean2[ib],2)/(nmean/2);
                     wsd2[ib]/=(nmean/2-1);
                     wmean2[ib]/=(nmean/2);

                     samplevar[ib]=sqrt((wsd1[ib]+wsd2[ib])/nmean);

                     tval=(wmean1[ib]-wmean2[ib])/samplevar[ib];

                     cout << "base " << ib << " with tval " << abs(tval) << " for cutoff " << tcutoff << ": ";
                     bool isok=(abs(tval)>tcutoff);
                     cout << isok << endl; 

                     if (!isok) isaccepted=0;

                  }
                  if (isaccepted) return wmean1;
                  gotonext=0;
               }
            }
            counter++;
            wprev=wafter;
            fprev=fafter;
         }

         iter++;
      }
   }

   vd resmean;
   resmean.push_back(log(wfinal[0]/conca));
   resmean.push_back(log(wfinal[1]/conct));
   resmean.push_back(log(wfinal[2]/concc));
   resmean.push_back(log(wfinal[3]/concg));

   return resmean;
}

// Tests convergence for different number of paths for total mean
   vd
colmean_test_numiter(unsigned int pos,Motif * mot)
{
   //ofstream outcv("test_numiter.dat");
   ofstream outcv("dum.dat");

   vd dum(4,0.);
   vd winit(4,0.);
   winit[0]=conca*exp(mot->matprec[pos][0]);
   winit[1]=conct*exp(mot->matprec[pos][1]);
   winit[2]=concc*exp(mot->matprec[pos][2]);
   winit[3]=concg*exp(mot->matprec[pos][3]);
   vd wprevmean(4,0.);
   wprevmean=winit;
   double Rmax(.013),rad,theta,phi; // for random walk
   double fprev,fafter;
   unsigned int nsample(20); // vary nsample for auto-correlation test, nmean is the size of the sliding window
   void * par[2]={(void *)mot,&pos};
   vd res;

   //for (unsigned int nforcvmean=1;nforcvmean<=3001;nforcvmean+=50){
   for (unsigned int nforcvmean=1000;nforcvmean<=1001;nforcvmean+=50){

      unsigned int nmean(nforcvmean); // vary nsample for auto-correlation test, nmean is the size of the sliding window
      vd f(nmean,0.); // for f mean computation, taken every nsample steps
      //vd dw(nmean,0.); // integration volumes at each step
      vvd w(nmean,dum); // for w mean computation, taken every nsample steps
      vint nrej(nmean,0); // for rejection rate, taken every step

      vd wfinal(4,0.);
      unsigned int numtotmean(1000);
      for (unsigned int totmean=0;totmean<numtotmean;totmean++){

         vd wprev(4,0.),wafter(4,0.),eps(4,0.),wmean(4,0.);
         // Starting point
         wprev=winit;

         fprev=likelyhood(wprev,par);
         f[0]=fprev;
         w[0]=wprev;
         nrej[0]=1;
         // =========
         // MAIN LOOP
         // =========

         double alph,prob; // for Metropolis Hastings
         // counter counts accepted trials 
         // tc counts accepted trials every nsample steps modulo nmean
         // iter counts every step
         unsigned int counter(1),tc(0),iter(1);
         while (1) {

            // ===========
            // Random walk
            // ===========
            //  PROPOSAL = UNIFORM IN A BALL
            rad=Rmax*gsl_rng_uniform (gslran);
            theta=M_PI*gsl_rng_uniform (gslran);
            phi=2*M_PI*gsl_rng_uniform (gslran);
            eps[0]=rad*sin(theta)*cos(phi);
            eps[1]=rad*sin(theta)*sin(phi);
            eps[2]=rad*cos(theta);
            eps[3]=-(eps[0]+eps[1]+eps[2]);

            wafter=wprev+eps;

            // check if new probabilities are consistent
            int wtest(0);
            for (int i=0;i<4;i++) if (wafter[i]<0 || wafter[i]>1) wtest=1;
            if (wtest) continue;


            fafter=likelyhood(wafter,par);

            // Metropolis-Hastings algorithm
            alph=min(1.,fafter/fprev);
            prob=gsl_rng_uniform (gslran);

            if (prob<=alph){
               // we define nsample to avoid auto-correlation
               // we want one point every nsample, nmean times
               tc=counter%(nmean*nsample);
               if (tc%nsample==0){
                  tc/=nsample;
                  unsigned int prevcounter;
                  if (tc==0) prevcounter=nmean-1;
                  else prevcounter=tc-1;
                  f[tc]=fafter;
                  w[tc]=wafter;

                  if (tc==0){
                     for (int itc=0;itc<nmean;itc++){
                        for (int ib=0;ib<3;ib++){
                           wmean[ib]+=w[itc][ib]/nmean;
                        }  
                     }
                     wmean[3]=1.-wmean[0]-wmean[1]-wmean[2];
                     break;
                  }
               }
               counter++;
               wprev=wafter;
               fprev=fafter;
            }
         }

         for (unsigned int ib=0;ib<4;ib++){
            wfinal[ib]+=wmean[ib]/numtotmean;
         }
      }
      double prevdistKL=0.;
      for (int ib=0;ib<4;ib++){
         prevdistKL+=wfinal[ib]*log(wfinal[ib]/wprevmean[ib]);
      }

      wprevmean=wfinal;

      double totdistKL=0.;
      for (int ib=0;ib<4;ib++){
         totdistKL+=wfinal[ib]*log(wfinal[ib]/winit[ib]);
      }

      outcv << nforcvmean << "\t";
      outcv << prevdistKL << "\t";
      outcv << totdistKL << "\t";
      for (int ib=0;ib<4;ib++){
         outcv << wfinal[ib] << "\t";
      }
      outcv << endl;
      res=wfinal;
   }


   outcv.close();
   //   res.push_back(log(wmean[0]/conca));
   //   res.push_back(log(wmean[1]/conct));
   //   res.push_back(log(wmean[2]/concc));
   //   res.push_back(log(wmean[3]/concg));

   return res;
}

// Tests convergence for different number of paths for total mean
   vd
colmean_test_numiter_w_dirichlet_metropolis(unsigned int pos,Motif * mot)
{
   ofstream outcv("test_numiter_dir.dat");
   //ofstream outcv("dum.dat");

   //FOR WINIT
   //unsigned int Ndir(20);
   //FOR WPREV
   unsigned int Ndir(30);
   vd dum(4,0.);
   vd winit(4,0.);
   winit[0]=conca*exp(mot->matprec[pos][0]);
   winit[1]=conct*exp(mot->matprec[pos][1]);
   winit[2]=concc*exp(mot->matprec[pos][2]);
   winit[3]=concg*exp(mot->matprec[pos][3]);
   double fprev,fafter;
   void * par[2]={(void *)mot,&pos};
   vd res;

   for (unsigned int nforcvmean=1;nforcvmean<=150;nforcvmean+=1){
      unsigned int nmean(nforcvmean); // vary nsample for auto-correlation test, nmean is the size of the sliding window
      vd wfinal(4,0.);
      unsigned int numtotmean(5000);

      for (unsigned int totmean=0;totmean<numtotmean;totmean++){

         vd wprev(4,0.),wafter(4,0.);
         // Starting point
         wprev=winit;

         fprev=likelyhood(wprev,par);
         // =========
         // MAIN LOOP
         // =========

         double alph,prob; // for Metropolis Hastings
         // counter counts accepted trials 
         // tc counts accepted trials every nsample steps modulo nmean
         // iter counts every step
         unsigned int counter(2);
         while (counter<=nmean) {

            double theta[]={0.,0.,0.,0.};
            // IF USING WPREV, THEN USE PROPOSAL IN REJECTION BELOW
            double alp[4]={alpha+Ndir*wprev[0],alpha+Ndir*wprev[1],beta+Ndir*wprev[2],beta+Ndir*wprev[3]};
            //double alp[4]={Ndir*winit[0],Ndir*winit[1],Ndir*winit[2],Ndir*winit[3]};
            gsl_ran_dirichlet(gslran,4,alp,theta);
            wafter[0]=theta[0];
            wafter[1]=theta[1];
            wafter[2]=theta[2];
            wafter[3]=1-theta[0]-theta[1]-theta[2];

            // check if new probabilities are consistent
            //            int wtest(0);
            //            for (int i=0;i<4;i++) if (wafter[i]<1e-10) wtest=1;
            //            if (wtest) continue;

            fafter=likelyhood(wafter,par);

            // Metropolis-Hastings
            // Move to new position with probability alpha, else reject
            double theta2[]={wprev[0],wprev[1],wprev[2],wprev[3]};
            double alp2[4]={alpha+Ndir*wafter[0],alpha+Ndir*wafter[1],beta+Ndir*wafter[2],beta+Ndir*wafter[3]};
            alph=min(1.,fafter*gsl_ran_dirichlet_pdf(4,alp2,theta2)/(fprev*gsl_ran_dirichlet_pdf(4,alp,theta)));

            //alph=min(1.,fafter/fprev);
            prob=gsl_rng_uniform (gslran);

            if (prob<=alph){
               // we define nsample to avoid auto-correlation
               // we want one point every nsample, nmean times
               for (unsigned int ib=0;ib<4;ib++){
                  wprev[ib]=(counter-1.)/counter*wprev[ib]+wafter[ib]/counter;
               }

               counter++;
               wprev=wafter;
               fprev=fafter;
            }
         }

         for (unsigned int ib=0;ib<4;ib++){
            wfinal[ib]+=wprev[ib]/numtotmean;
         }
      }

      double totdistKL=0.;
      for (int ib=0;ib<4;ib++){
         totdistKL+=wfinal[ib]*log(wfinal[ib]/winit[ib]);
      }

      outcv << nforcvmean << "\t";
      outcv << totdistKL << "\t";
      for (int ib=0;ib<4;ib++){
         outcv << wfinal[ib] << "\t";
      }
      outcv << endl;
      res=wfinal;
   }

   outcv.close();

   return res;

}

// Tests convergence for different number of paths for total mean
   vd
colmean_test_numiter_w_dirichlet(unsigned int pos,Motif * mot)
{
   ofstream outcv("test_numiter_dir.dat");

   unsigned int Ndir(0);
   for (ivma ivm=mot->seqs.begin();ivm!=mot->seqs.end();ivm++){
      Ndir+=accumulate(ivm->matches.begin(),ivm->matches.end(),0)/2; 

   }
   cout << "Ndir=" << Ndir << endl;

   unsigned int nsample(20);
   cout << "nsample=" << nsample << endl;

   //FOR WINIT
   //unsigned int Ndir(20);
   //FOR WPREV
   //unsigned int Ndir(1500);
   vd dum(4,0.);
   vd winit(4,0.);
   winit[0]=conca*exp(mot->matprec[pos][0]);
   winit[1]=conct*exp(mot->matprec[pos][1]);
   winit[2]=concc*exp(mot->matprec[pos][2]);
   //winit[3]=concg*exp(mot->matprec[pos][3]);
   winit[3]=1.-winit[0]-winit[1]-winit[2];
   double fprev,fafter;
   void * par[2]={(void *)mot,&pos};
   vd res;

   outcv << 0 << "\t";
   outcv << 0 << "\t";
   outcv << 0 << "\t";
   for (int ib=0;ib<4;ib++){
      outcv << winit[ib] << "\t";
   }
   outcv << endl;

   vd wfinalprev(4,0.);
   wfinalprev=winit;

   for (unsigned int nforcvmean=2;nforcvmean<=150;nforcvmean+=1){
      unsigned int nmean(nforcvmean); // vary nsample for auto-correlation test, nmean is the size of the sliding window
      vd wfinal(4,0.);
      unsigned int numtotmean(1);

      for (unsigned int totmean=1;totmean<=numtotmean;totmean++){

         double fmean(0.);
         vd wmean(4,0.),wprev(4,0.),wafter(4,0.);
         vd wtmpmean(4,0.);
         // Starting point
         //wprev=winit;
         //
         //fprev=likelyhood(wprev,par);

         double alpinit[4]={alpha+Ndir*winit[0],alpha+Ndir*winit[1],beta+Ndir*winit[2],beta+Ndir*winit[3]};
         double thetainit[]={0.,0.,0.,0.};
         gsl_ran_dirichlet(gslran,4,alpinit,thetainit);
         wprev[0]=thetainit[0];
         wprev[1]=thetainit[1];
         wprev[2]=thetainit[2];
         wprev[3]=1-thetainit[0]-thetainit[1]-thetainit[2];
         fprev=loglikelyhood(wprev,par);

         double totdistKL=0.;
         for (int ib=0;ib<4;ib++){
            totdistKL+=wprev[ib]*log(wprev[ib]/winit[ib]);
         }

         outcv << 1 << "\t";
         outcv << totdistKL << "\t";
         outcv << totdistKL << "\t";
         for (int ib=0;ib<4;ib++){
            outcv << wprev[ib] << "\t";
         }
         outcv << endl;

         // =========
         // MAIN LOOP
         // =========

         double logalph,prob; // for Metropolis Hastings
         // counter counts accepted trials 
         unsigned int counter(1),iterok(1);
         // BURN-IN
         unsigned int bicount(0),bicutoff(nsample);

         while (counter<=nmean) {

            double alp[4]={alpha+Ndir*wprev[0],alpha+Ndir*wprev[1],beta+Ndir*wprev[2],beta+Ndir*wprev[3]};
            double theta[]={0.,0.,0.,0.};
            gsl_ran_dirichlet(gslran,4,alp,theta);
            wafter[0]=theta[0];
            wafter[1]=theta[1];
            wafter[2]=theta[2];
            wafter[3]=1-theta[0]-theta[1]-theta[2];

            fafter=loglikelyhood(wafter,par);

            // Metropolis-Hastings
            // Move to new position with probability alpha, else reject
            double theta2[]={wprev[0],wprev[1],wprev[2],wprev[3]};
            double alp2[4]={alpha+Ndir*wafter[0],alpha+Ndir*wafter[1],beta+Ndir*wafter[2],beta+Ndir*wafter[3]};
            logalph=min(0.,fafter+log(gsl_ran_dirichlet_pdf(4,alp2,theta2))-(fprev+log(gsl_ran_dirichlet_pdf(4,alp,theta))));
            //alph=min(1.,fafter/fprev);
            prob=gsl_rng_uniform (gslran);

            if (log(prob)<=logalph){

               // we define nsample to avoid auto-correlation
               // we want one point every nsample, nmean times
               wprev=wafter;
               fprev=fafter;

               if (bicount<bicutoff){
                  bicount++;
                  continue;
               }

               if (iterok%nsample==0){
                  for (unsigned int ib=0;ib<4;ib++){
                     wmean[ib]+=wprev[ib]/nmean;
                  }
                  counter++;
               }
               iterok++;
            }
         }

         for (unsigned int ib=0;ib<4;ib++){
            wfinal[ib]+=wmean[ib]/numtotmean;
         }
      }

      double totdistKL=0.;
      for (int ib=0;ib<4;ib++){
         totdistKL+=wfinal[ib]*log(wfinal[ib]/winit[ib]);
      }
      double prevdistKL=0.;
      for (int ib=0;ib<4;ib++){
         prevdistKL+=wfinal[ib]*log(wfinal[ib]/wfinalprev[ib]);
      }
      wfinalprev=wfinal;

      outcv << nforcvmean << "\t";
      outcv << totdistKL << "\t";
      outcv << prevdistKL << "\t";
      for (int ib=0;ib<4;ib++){
         outcv << wfinal[ib] << "\t";
      }
      outcv << endl;
      res=wfinal;
   }

   outcv.close();

   return res;
}

   vd
colmean_testCV(unsigned int pos,Motif * mot)
{
   unsigned int nforcvmean=20;

   ofstream outcv("testCV.dat");
   vint vmean;
   vmean.push_back(20);
   vmean.push_back(40);
   for (int iex=50;iex<1050;iex+=50) vmean.push_back(iex);

   vd vdcutoff;
   vdcutoff.push_back(.0008);
   vdcutoff.push_back(.0009);
   for (int iex=-3;iex<=-1;iex++){
      for (int i=1;i<10;i++) vdcutoff.push_back(i*pow(10,iex));
   }

   vd dum(4,0.);
   vd winit(4,0.);
   winit[0]=conca*exp(mot->matprec[pos][0]);
   winit[1]=conct*exp(mot->matprec[pos][1]);
   winit[2]=concc*exp(mot->matprec[pos][2]);
   winit[3]=concg*exp(mot->matprec[pos][3]);
   double Rmax(.013),rad,theta,phi; // for random walk
   double fprev,fafter;
   unsigned int nsample(20); // vary nsample for auto-correlation test, nmean is the size of the sliding window
   void * par[2]={(void *)mot,&pos};


   for (ivint ivi=vmean.begin();ivi!=vmean.end();ivi++){
      for (ivd iv=vdcutoff.begin();iv!=vdcutoff.end();iv++){
         vd wfinal(4,0.);
         for (unsigned int imean=0;imean<nforcvmean;imean++){

            unsigned int nmean(*ivi); // vary nsample for auto-correlation test, nmean is the size of the sliding window
            double dcutoff(*iv); // loop ending requirement
            vd f(nmean,0.); // for f mean computation, taken every nsample steps
            //vd dw(nmean,0.); // integration volumes at each step
            vvd w(nmean,dum); // for w mean computation, taken every nsample steps
            vint nrej(nmean,0); // for rejection rate, taken every step
            double distKL(10.); // loop ending requirement
            vd wprev(4,0.),wafter(4,0.),eps(4,0.),wmean(4,0.),wprevmean(4,0.);

            vd wbest(4,0.);
            double fbest(0.);

            // Starting point
            wprev=winit;

            fprev=likelyhood(wprev,par);
            f[0]=fprev;
            w[0]=wprev;
            wprevmean=wprev;
            nrej[0]=1;

            // =========
            // MAIN LOOP
            // =========

            double alph,prob; // for Metropolis Hastings
            // counter counts accepted trials 
            // tc counts accepted trials every nsample steps modulo nmean
            // iter counts every step
            unsigned int counter(1),tc(0),iter(1);

            //ofstream outl("likelyhood.dat");
            while (abs(distKL)>dcutoff) {

               // ===========
               // Random walk
               // ===========
               //  PROPOSAL = UNIFORM IN A BALL
               rad=Rmax*gsl_rng_uniform (gslran);
               theta=M_PI*gsl_rng_uniform (gslran);
               phi=2*M_PI*gsl_rng_uniform (gslran);
               eps[0]=rad*sin(theta)*cos(phi);
               eps[1]=rad*sin(theta)*sin(phi);
               eps[2]=rad*cos(theta);
               eps[3]=-(eps[0]+eps[1]+eps[2]);

               wafter=wprev+eps;

               // check if new probabilities are consistent
               int wtest(0);
               for (int i=0;i<4;i++) if (wafter[i]<0 || wafter[i]>1) wtest=1;
               if (wtest) continue;


               fafter=likelyhood(wafter,par);

               // Metropolis-Hastings algorithm
               alph=min(1.,fafter/fprev);
               prob=gsl_rng_uniform (gslran);

               if (prob<=alph){
                  // we define nsample to avoid auto-correlation
                  // we want one point every nsample, nmean times
                  tc=counter%(nmean*nsample);
                  if (tc%nsample==0){
                     tc/=nsample;
                     unsigned int prevcounter;
                     if (tc==0) prevcounter=nmean-1;
                     else prevcounter=tc-1;
                     f[tc]=fafter;
                     w[tc]=wafter;
                     //            for (int ib=0;ib<4;ib++){
                     //               outl << w[tc][ib] << " ";
                     //            }  
                     //            outl << fafter <<endl;


                     if (fafter>fbest){
                        fbest=fafter;
                        wbest=wafter;
                     }

                     vd fnum(4,0.); 
                     double fdenom(0.);
                     for (int itc=0;itc<nmean;itc++){
                        for (int ib=0;ib<3;ib++){
                           fnum[ib]+=w[itc][ib]*f[itc];
                        }  
                        fdenom+=f[itc];
                     }

                     if (tc==0){
                        distKL=0.;
                        for (int ib=0;ib<3;ib++){
                           wmean[ib]=fnum[ib]/fdenom;
                           distKL+=wmean[ib]*log(wmean[ib]/wprevmean[ib]);
                        }
                        wmean[3]=1.-wmean[0]-wmean[1]-wmean[2];
                        distKL+=wmean[3]*log(wmean[3]/wprevmean[3]);
                        wprevmean=wmean;
                     }
                  }
                  counter++;
                  wprev=wafter;
                  fprev=fafter;
               }

               iter++;
            }
            //outl.close();
            //

            for (unsigned int ib=0;ib<4;ib++){
               wfinal[ib]+=wmean[ib]/nforcvmean;
            }
         }

         double totdistKL=0.;
         for (int ib=0;ib<4;ib++){
            totdistKL+=wfinal[ib]*log(wfinal[ib]/winit[ib]);
         }
         outcv << *ivi << "\t";
         outcv << *iv << "\t";
         outcv << totdistKL << endl;

      }
   }

   vd res;
   //   res.push_back(log(wmean[0]/conca));
   //   res.push_back(log(wmean[1]/conct));
   //   res.push_back(log(wmean[2]/concc));
   //   res.push_back(log(wmean[3]/concg));

   return res;
}
   
   void
colmean_testCV_w_dirichlet(unsigned int pos,Motif * mot)
{
   unsigned int Ndir(0);
   for (ivma ivm=mot->seqs.begin();ivm!=mot->seqs.end();ivm++){
      Ndir+=accumulate(ivm->matches.begin(),ivm->matches.end(),0)/2; 

   }
   cout << "Ndir=" << Ndir << endl;

   vd wprev(4,0.),wafter(4,0.);
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
   unsigned int nsample(30);
   void * par[2]={(void *)mot,&pos};

   //    // This is for autocorrelation computation
   ofstream outf("testCV.dat");
   //      ofstream outf("autocorrelation_0.4_sphere.dat"); 
   for (unsigned int nmean=10;nmean<1500;nmean+=50){
      
      vvd wmean,wvar; 

      // Starting point
      double alpinit[4]={alpha+Ndir*winit[0],alpha+Ndir*winit[1],beta+Ndir*winit[2],beta+Ndir*winit[3]};
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

      double logalph,prob; // for Metropolis Hastings
      // BURN-IN
      unsigned int counter(1),iter(1),iterok(1);

      // here we don't need more than nmean
      while (iterok<nmean*nsample) {

         double alp[4]={alpha+Ndir*wprev[0],alpha+Ndir*wprev[1],beta+Ndir*wprev[2],beta+Ndir*wprev[3]};
         double theta[]={0.,0.,0.,0.};
         gsl_ran_dirichlet(gslran,4,alp,theta);
         wafter[0]=theta[0];
         wafter[1]=theta[1];
         wafter[2]=theta[2];
         wafter[3]=1-theta[0]-theta[1]-theta[2];

         // check if new probabilities are consistent, else reject
         int wtest(0);
         for (int i=0;i<4;i++) if (wafter[i]<=0 || wafter[i]>1) wtest=1;
         if (wtest)  continue;

         fafter=loglikelyhood(wafter,par);

         // Metropolis-Hastings
         // Move to new position with probability alpha, else reject
         double theta2[]={wprev[0],wprev[1],wprev[2],wprev[3]};
         double alp2[4]={alpha+Ndir*wafter[0],alpha+Ndir*wafter[1],beta+Ndir*wafter[2],beta+Ndir*wafter[3]};
         logalph=min(0.,fafter+log(gsl_ran_dirichlet_pdf(4,alp2,theta2))-(fprev+log(gsl_ran_dirichlet_pdf(4,alp,theta))));
         prob=gsl_rng_uniform (gslran);

         if (log(prob)<=logalph){
            // auto-correlation with nsample size
            if (iterok%nsample==0){
               wmean.push_back(wafter);
               vd wvartmp(4,0.);
               for (unsigned int ib=0;ib<4;ib++){
                  wvartmp[ib]=pow(wafter[ib],2);
               }
               wvar.push_back(wvartmp);
               counter++;
            }
            iterok++;

            wprev=wafter;
            fprev=fafter;
         }
      }
      outf << nmean << "\t";
      vd wmeanf(4,0.),wvarf(4,0.);
      for (ivvd iv=wmean.begin();iv!=wmean.end();iv++){
         for (unsigned int ib=0;ib!=4;ib++){
            wmeanf[ib]+=(*iv)[ib]/counter;
         }
      }
      for (ivvd iv=wvar.begin();iv!=wvar.end();iv++){
         for (unsigned int ib=0;ib!=4;ib++){
            wvarf[ib]+=(*iv)[ib]/counter;
         }
      }
      for (unsigned int ib=0;ib!=4;ib++){
         wvarf[ib]-=pow(wmeanf[ib],2);
         wvarf[ib]*=counter/(counter-1);
      }
      
      for (unsigned int ib=0;ib!=4;ib++){
         outf << wmeanf[ib] << "\t" << wvarf[ib] << "\t";
      }
      outf << "\n";
      cout << nmean << "\t" << wmeanf;
   }
   outf.close();

   return;
}

   vd
colmean_test_likelihood_w_dirichlet(unsigned int pos,Motif * mot)
{
   unsigned int Ndir(0);
   //               for (ivma ivm=mot->seqs.begin();ivm!=mot->seqs.end();ivm++){
   //                  Ndir+=accumulate(ivm->matches.begin(),ivm->matches.end(),0)/2; 
   //
   //               }
   Ndir=0;
   cout << "Ndir=" << Ndir << endl;

   vd dum(4,0.);
   vd winit(4,0.);
   winit[0]=conca*exp(mot->matprec[pos][0]);
   winit[1]=conct*exp(mot->matprec[pos][1]);
   winit[2]=concc*exp(mot->matprec[pos][2]);
   winit[3]=1.-winit[0]-winit[1]-winit[2];
   double fprev,fafter;
   void * par[2]={(void *)mot,&pos};

   unsigned int nmean(20000); // vary nsample for auto-correlation test, nmean is the size of the sliding window

   vd w(4,0.);
   double f(0.);
   // Starting point
   w=winit;
   f=likelyhood(w,par);

   // =========
   // MAIN LOOP
   // =========

   // counter counts accepted trials 
   unsigned int counter(1);
   ofstream outl("likelyhood.dat");
   while (counter<=nmean) {

      double theta[]={0.,0.,0.,0.};
      //double alp[4]={alpha,alpha,beta,beta};
      double alp[4]={1.,1.,1.,1.};
      //double alp[4]={alpha+Ndir*winit[0],alpha+Ndir*winit[1],beta+Ndir*winit[2],beta+Ndir*winit[3]};
      gsl_ran_dirichlet(gslran,4,alp,theta);
      w[0]=theta[0];
      w[1]=theta[1];
      w[2]=theta[2];
      w[3]=1.-theta[0]-theta[1]-theta[2];

      f=loglikelyhood(w,par);

      outl << w[0] << " ";
      outl << w[1] << " ";
      outl << w[2] << " ";
      outl << w[3] << " ";
      outl << f << "\n";

      counter++;
   } 
   outl.close();
   vd res;

   return res;
}

   vd
colmean_test_rejection_w_dirichlet_metropolis(unsigned int pos,Motif * mot)
{
   double fprev,fafter;
   unsigned int nsample(1),nmean(300); // vary nsample for auto-correlation test, nmean is the size of the sliding window
   void * par[2]={(void *)mot,&pos};

   vd winit(4,0.);
   // Starting point
   winit[0]=conca*exp(mot->matprec[pos][0]);
   winit[1]=conct*exp(mot->matprec[pos][1]);
   winit[2]=concc*exp(mot->matprec[pos][2]);
   winit[3]=concg*exp(mot->matprec[pos][3]);

   // This is for rejection rate computation
   ofstream outf("rejection_rate_ball_dirichlet.dat"); 
   unsigned int Ndir;
   for (unsigned int Ndir=1;Ndir<3000;Ndir+=100){
      cout << "Ndir=" << Ndir << endl;
      vd wprev(4,0.),wafter(4,0.),eps(4,0.),wmean(4,0.);
      vd dum(4,0.);
      int nrej(0);

      // Starting point
      //
      // PB : FAR TOO BIG LIKELIHOOD TO BEGIN WITH!
      //      wprev=winit;
      //      fprev=likelyhood(wprev,par);
      //      cout << fprev << endl;

      double alpinit[4]={alpha+Ndir*winit[0],alpha+Ndir*winit[1],beta+Ndir*winit[2],beta+Ndir*winit[3]};
      //double alpinit[4]={alpha,alpha,beta,beta};
      //double alp[4]={Ndir*alpha,Ndir*alpha,Ndir*beta,Ndir*beta};
      //double alp[4]={alpha+Ndir*winit[0],alpha+Ndir*winit[1],beta+Ndir*winit[2],beta+Ndir*winit[3]};
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

      double logalph,prob; // for Metropolis Hastings
      // counter counts accepted trials every nsample steps modulo nmean
      // iter counts every step
      // iterok counts accepted
      unsigned int iter(1),iterok(1);

      // here we don't need more than nmean
      while (iter<=nmean) {

         // IF USING WPREV, THEN USE PROPOSAL IN REJECTION BELOW
         //double alp[4]={Ndir*wprev[0],Ndir*wprev[1],Ndir*wprev[2],Ndir*wprev[3]};
         double alp[4]={alpha+Ndir*wprev[0],alpha+Ndir*wprev[1],beta+Ndir*wprev[2],beta+Ndir*wprev[3]};
         //double alp[4]={Ndir*alpha,Ndir*alpha,Ndir*beta,Ndir*beta};
         //double alp[4]={alpha+Ndir*winit[0],alpha+Ndir*winit[1],beta+Ndir*winit[2],beta+Ndir*winit[3]};
         double theta[]={0.,0.,0.,0.};
         gsl_ran_dirichlet(gslran,4,alp,theta);
         wafter[0]=theta[0];
         wafter[1]=theta[1];
         wafter[2]=theta[2];
         wafter[3]=1-theta[0]-theta[1]-theta[2];

         fafter=loglikelyhood(wafter,par);

         // Metropolis-Hastings
         // Move to new position with probability alpha, else reject
         double theta2[]={wprev[0],wprev[1],wprev[2],wprev[3]};
         double alp2[4]={alpha+Ndir*wafter[0],alpha+Ndir*wafter[1],beta+Ndir*wafter[2],beta+Ndir*wafter[3]};
         logalph=min(0.,fafter+log(gsl_ran_dirichlet_pdf(4,alp2,theta2))-(fprev+log(gsl_ran_dirichlet_pdf(4,alp,theta))));
         //alph=min(1.,fafter/fprev);
         prob=gsl_rng_uniform (gslran);


         if (log(prob)>logalph) nrej++;
         else {
            // every nsample to avoid auto-correlation
            wprev=wafter;
            for (unsigned int ib=0;ib<4;ib++){
               wmean[ib]+=wprev[ib]/nmean;
            }
            fprev=fafter;
            iterok++;
         }

         iter++;
      }
      outf << Ndir << "\t" << (double)nrej/iter << endl;
      cout << "iterOK=" << iterok << endl;
   }
   outf.close();

   exit(2);
   vd res;

   return res;
}

   void
colmean_test_autocorr_w_dirichlet(unsigned int pos,Motif * mot)
{
   unsigned int Ndir(0);
   for (ivma ivm=mot->seqs.begin();ivm!=mot->seqs.end();ivm++){
      Ndir+=accumulate(ivm->matches.begin(),ivm->matches.end(),0)/2; 

   }
   cout << "Ndir=" << Ndir << endl;

   vd wprev(4,0.),wafter(4,0.);
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
   unsigned int nsample(50),nmean(5000);
   void * par[2]={(void *)mot,&pos};

   //    // This is for autocorrelation computation
   ofstream outf("autocorrelation_dir.dat"); 
   //      ofstream outf("autocorrelation_0.4_sphere.dat"); 
   for (unsigned int isample=1;isample<nsample;isample++){
      
      double fnfk(0.),fn(0.),fk(0.),fn2(0.),fk2(0.),fprevmean(0.);
      vd wnwk(4,0.),wn(4,0),wk(4,0),wprevmean(4,0.); 
      unsigned int ncorr(0);

      // Starting point
      double alpinit[4]={alpha+Ndir*winit[0],alpha+Ndir*winit[1],beta+Ndir*winit[2],beta+Ndir*winit[3]};
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

      double logalph,prob; // for Metropolis Hastings
      // BURN-IN
      unsigned int counter(1),iter(1),iterok(1);

      // here we don't need more than nmean
      while (iterok<nmean*isample) {

         double alp[4]={alpha+Ndir*wprev[0],alpha+Ndir*wprev[1],beta+Ndir*wprev[2],beta+Ndir*wprev[3]};
         double theta[]={0.,0.,0.,0.};
         gsl_ran_dirichlet(gslran,4,alp,theta);
         wafter[0]=theta[0];
         wafter[1]=theta[1];
         wafter[2]=theta[2];
         wafter[3]=1-theta[0]-theta[1]-theta[2];

         // check if new probabilities are consistent, else reject
         int wtest(0);
         for (int i=0;i<4;i++) if (wafter[i]<=0 || wafter[i]>1) wtest=1;
         if (wtest)  continue;

         fafter=loglikelyhood(wafter,par);

         // Metropolis-Hastings
         // Move to new position with probability alpha, else reject
         double theta2[]={wprev[0],wprev[1],wprev[2],wprev[3]};
         double alp2[4]={alpha+Ndir*wafter[0],alpha+Ndir*wafter[1],beta+Ndir*wafter[2],beta+Ndir*wafter[3]};
         logalph=min(0.,fafter+log(gsl_ran_dirichlet_pdf(4,alp2,theta2))-(fprev+log(gsl_ran_dirichlet_pdf(4,alp,theta))));
         prob=gsl_rng_uniform (gslran);

         if (log(prob)<=logalph){
            // auto-correlation with nsample size
            if (iterok%isample==0){
               fnfk+=fafter*fprevmean;
               fn+=fafter;
               fn2+=pow(fafter,2);
               fk+=fprevmean;
               fk2+=pow(fprevmean,2);
               fprevmean=fafter;

               for (unsigned int i=0;i<4;i++){
                  wnwk[i]+=wafter[i]*wprevmean[i];
                  wn[i]+=wafter[i];
                  wk[i]+=wprevmean[i];
               }
               wprevmean=wafter;

               ncorr++;
               counter++;
            }
            iterok++;

            wprev=wafter;
            fprev=fafter;
         }

      }
      outf << isample << "\t" << 1./ncorr*(fnfk-fn*fk/ncorr) << "\t";
      cout << isample << "\t" << 1./ncorr*(fnfk-fn*fk/ncorr) << "\t";
      for (unsigned int i=0;i<4;i++){
         outf << 1./ncorr*(wnwk[i]-wn[i]*wk[i]/ncorr) << "\t";
         cout << 1./ncorr*(wnwk[i]-wn[i]*wk[i]/ncorr) << "\t";
      }
      outf << "\n";
      cout << "\n";
   }
   outf.close();

   return;
}

   vd
colmean_test_rejection_n_autocorr(unsigned int pos,Motif * mot)
{
   double Rmax(0.35),rad,theta,phi; // for random walk
   vd wprev(4,0.),wafter(4,0.),eps(4,0.);
   double fprev,fafter;
   unsigned int nsample(10),nmean(200); // vary nsample for auto-correlation test, nmean is the size of the sliding window
   void * par[2]={(void *)mot,&pos};

   //    // This is for autocorrelation computation
   ofstream outf("autocorrelation_0.35_ball.dat"); 
   //      ofstream outf("autocorrelation_0.4_sphere.dat"); 
   for (unsigned int isample=5;isample<100;isample++){
      double fnfk(0.),fn(0.),fk(0.),fn2(0.),fk2(0.),fprevmean(0.);
      unsigned int ncorr(0);
      nsample=isample;

      double fmean(0.);

      // vd fmean(nmean,0.); // for f mean computation, taken every nsample steps
      vd dum(4,0.);
      vvd wmean(nmean,dum); // for w mean computation, taken every nsample steps

      // This is for rejection rate computation
      //   ofstream outf("rejection_rate_ball.dat");    
      //   for (int ipow=1;ipow<400;ipow++){}
      //      Rmax=pow(50.,-(double)ipow/100);
      int nrej(0);

      // Starting point
      wprev[0]=conca*exp(mot->matprec[pos][0]);
      wprev[1]=conct*exp(mot->matprec[pos][1]);
      wprev[2]=concc*exp(mot->matprec[pos][2]);
      wprev[3]=concg*exp(mot->matprec[pos][3]);

      fprev=likelyhood(wprev,par);
      //fmean[0]=fprev;
      fmean=fprev;
      fprevmean=fprev;
      wmean[0]=wprev;

      // =========
      // MAIN LOOP
      // =========

      double alph,prob; // for Metropolis Hastings
      // counter counts accepted trials every nsample steps modulo nmean
      // rejcounter counts all trials modulo nmean
      // iter counts every step
      // iterok counts accepted
      unsigned int counter(1),iter(1),iterok(1);

      // here we don't need more than nmean
      while (iterok<nmean*nsample) {

         // Random walk
         //  UNCOMMENT FOR BALL
         rad=Rmax*gsl_rng_uniform (gslran);
         //      rad=Rmax;
         theta=M_PI*gsl_rng_uniform (gslran);
         phi=2*M_PI*gsl_rng_uniform (gslran);
         eps[0]=rad*sin(theta)*cos(phi);
         eps[1]=rad*sin(theta)*sin(phi);
         eps[2]=rad*cos(theta);
         eps[3]=-(eps[0]+eps[1]+eps[2]);

         wafter=wprev+eps;

         // check if new probabilities are consistent, else reject
         int wtest(0);
         for (int i=0;i<4;i++) if (wafter[i]<0 || wafter[i]>1) wtest=1;
         if (wtest){
            nrej++;
            iter++;
            continue;
         }

         fafter=likelyhood(wafter,par);

         // Metropolis-Hastings
         // Move to new position with probability alpha, else reject
         alph=min(1.,fafter/fprev);
         prob=gsl_rng_uniform (gslran);

         if (prob>alph) nrej++;
         else {
            // every nsample to avoid auto-correlation
            cout << iterok << endl;
            if (iterok%nsample==0){
               //            // uncomment for autocorrelation computation
               //                              fnfk+=fmean[counter]*fmean[prevcounter];
               //                              fn+=fmean[counter];
               //                              fn2+=pow(fmean[counter],2);
               //                              fk+=fmean[prevcounter];
               //                              fk2+=pow(fmean[prevcounter],2);
               fnfk+=fafter*fprevmean;
               fn+=fafter;
               fn2+=pow(fafter,2);
               fk+=fprevmean;
               fk2+=pow(fprevmean,2);
               fprevmean=fafter;
               ncorr++;
               counter++;
               cout << fafter << endl;
               cout << fprevmean << endl;
            }
            iterok++;

            wprev=wafter;
            fprev=fafter;
         }

         iter++;

      }
      //outf << Rmax << "\t" << (double)nrej/iter << endl;
      outf << nsample << "\t" << (fnfk-fn*fk/ncorr)/sqrt((fk2-pow(fk,2)/ncorr)*(fn2-pow(fn,2)/ncorr)) << endl;
   }
   outf.close();

   exit(2);
   vd res;
   //   res.push_back(log(w0/conca));
   //   res.push_back(log(w1/conct));
   //   res.push_back(log(w2/concc));
   //   res.push_back(log(w3/concg));

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

   Motalign
mahamming(unsigned int pos, Sequence & seq, Motif & mot,int sens,unsigned int numhamm)
{
   Motalign ma;
   unsigned int motwidth=mot.motwidth;
   ma.seq_start=seq.iseqs[0].begin();
   ma.seq_stop=seq.iseqs[0].end();
   civvint imap=seq.imaps.begin()+1;
   int seqnum=1;
   ma.matches=vint(nbspecies,0);
   vint seqdum=vint(motwidth,4);
   ivint posdum;
   ma.matchespos.push_back(seq.iseqs[0].begin()+pos);
   vint seqmel=vint(seq.iseqs[0].begin()+pos,seq.iseqs[0].begin()+pos+motwidth);
   if (sens==1){
      ma.alignseq.push_back(seqmel);
      ma.strand=1;
   } else {
      ma.alignseq.push_back(reversecomp(seqmel));
      ma.strand=-1;
   }
   ma.matches[0]=1;
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
            vint seqstart=vint(startsite,startsite+motwidth);
            if (scorefhamming(seqmel,seqstart)<=numhamm){
               vint::const_iterator endci=(*is).end()-motwidth+1;
               unsigned int sh=shifthamming((*is).begin()+start,seqmel,endci,2*neighbext);
               vint seqt=vint((*is).begin()+start+sh,(*is).begin()+start+sh+motwidth);
               if (sens==1){
                  ma.alignseq.push_back(seqt);
               } else if (sens==-1){
                  ma.alignseq.push_back(reversecomp(seqt));
               }
               ma.matchespos.push_back((*is).begin()+start+sh);
               hasmatch=1;
               break;
            }
         }
         if (hasmatch){
            ma.matches[seqnum]=1;
         } else {
            ma.alignseq.push_back(seqdum);
            ma.matchespos.push_back(posdum);
            ma.matches[seqnum]=0;
         }
      }
      else {
         ma.alignseq.push_back(seqdum);
         ma.matchespos.push_back(posdum);
         ma.matches[seqnum]=0;
      }
      seqnum++;
      imap++;
   }
   //   cout << "\n";
   return ma;
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
      mot1.matenergy=mattoenergy(mot1.matprec);
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
loadmotswnames ( const char * filename, vmot & mots )
{
   if (!filename){
      cout << "Please give a motifs file. Exiting..." << endl;
      exit(1);
   }

   ifstream ffmotifs;
   ffmotifs.open(filename);

   string dum;
   getline(ffmotifs,dum);
   //fmotifs >> dum;

   unsigned int i(0);//!!to comply with cpp convention
   while (! ffmotifs.eof()){
      stringstream fmotifs(dum);
      Motif mot1;
      mot1.index=i;
      fmotifs >> mot1.name;
      fmotifs >> mot1.bsinit;
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
      mot1.matenergy=mattoenergy(mot1.matprec);
      mot1.matprecrevcomp=reversecomp(mot1.matprec);
      fmotifs >> mot1.distmot;
      fmotifs >> dum;
      mot1.motscorethr2=mot1.motwidth*scorethr2/10;
      mot1.motscorethr=mot1.motwidth*(scorethr2-1)/10;
      mot1.motscorethrcons=mot1.motwidth*(scorethr2-1)/10;
      mots.push_back(mot1);
      getline(ffmotifs,dum);
      i++;
   }
   ffmotifs.close();

   if (i==0){//No motifs
      cout << "No motifs" << endl;
      exit(1);
   } else if (i<nbmots_for_score-1) {
      nbmots_for_score=i;
      //cout << "Changed nbmots_for_score to value " << i << endl;
   }
}		/* -----  end of function loadmots  ----- */

   void
loadjaspardb ( vmot & mots )
{
   ifstream fmotifs;
   //fmotifs.open("/home/santolin/these/files/jaspar/jaspar_pwm_wname.dat");
   //fmotifs.open("/home/santolin/these/files/jaspar/jaspar_pwm_wname+personal.dat");
   fmotifs.open("/home/santolin/these/files/jaspar/jaspar+motgen+known.dat");
   string dum;
   fmotifs >> dum;
   unsigned int i(0);//!!to comply withh cpp convention
   while (! fmotifs.eof()){
      Motif mot1;
      mot1.index=i;
      mot1.name=dum;
      mot1.id=mot1.name.substr(0,mot1.name.find("_"));
      mot1.name=mot1.name.substr(mot1.name.find("_")+1,mot1.name.size());
      fmotifs >> mot1.bsinit;
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
      mot1.matenergy=mattoenergy(mot1.matprec);
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

   fmotifs.open("/home/santolin/these/files/transfac/matrices/all/pwms.dat");
   fmotifs >> dum;
   while (! fmotifs.eof()){
      Motif mot1;
      mot1.index=i;
      mot1.name=dum;
      mot1.id=mot1.name;
      fmotifs >> mot1.bsinit;
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
      mot1.matenergy=mattoenergy(mot1.matprec);
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
      score+=nbmots[imot];//*log((*ivm).lambdatrain/(*ivm).lambda);
      imot++;
   }
}

   void
GroupInstance::compscoreweight(vmot & lmots,unsigned int nbmots_score)
{
   score=0;
   unsigned int imot=0;
   for (ivmot ivm=lmots.begin();ivm!=lmots.begin()+nbmots_score;ivm++){
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
         //         for (ivTSS ivt=(*ivg).TSSs.begin();ivt!=(*ivg).TSSs.end();ivt++){}
         //            if ((*ivt).besttss.gene==*ivs){
         if ((*ivg).besttss.gene==*ivs){
            if ((*ivg).score>scoregene){
               scoregene=(*ivg).score;
            }
            //               break;
            //            }
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
   double x0, x = 0.1, r_expected = sqrt (5.0);
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


   fconf << "computed priors : alpha=" << alpha << ", beta=" << beta endl;
   gsl_root_fdfsolver_free (s);
   return status;
}

   void 
displaymat(vvd & mat)
{
   cout.precision(4);
   for (int i=0;i<4;i++){
      for (int j=0;j<mat.size();j++){
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



