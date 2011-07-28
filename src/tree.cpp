#include <iostream>
#include<iomanip>
#include <fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include	<cmath>
#include	<gsl/gsl_vector.h>
#include	<gsl/gsl_matrix.h>
#include	<gsl/gsl_blas.h>
#include	<gsl/gsl_odeiv.h>
#include<gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include <functional>
#include <numeric>

#include "tree.hpp"
#include "const.hpp"
#include "random.hpp" 
#include "scangen.hpp" 

using namespace std;

clock_t start,finish;
double dif;

vnoe treedist;

double pa;
double pc;
const double pback[4]={conca,conct,concc,concg};
const double punif[4]={0.25,0.25,0.25,0.25};
gsl_ran_discrete_t * gslback=gsl_ran_discrete_preproc (4,pback);
gsl_ran_discrete_t * gslunif=gsl_ran_discrete_preproc (4,punif);

const double kappa=2.0;

const double integr_step=0.01;//0.01;//=0.001; //!! 0.001 instead of 0.01
double fat,fac,fag,fta,ftc,ftg,fca,fct,fcg,fga,fgt,fgc;

noeud::noeud(int e1,int e2, int n, double p1, double p2)
{
   esp1=e1;
   esp2=e2;
   noe=n;
   if (evolutionary_model==2){
      prox1=p1;
      prox2=p2;
   } else if (evolutionary_model==1){
      double correction=0.5+4.0*conca*concc;
      //double correction=1.0;
      prox1=exp(-p1/correction);
      prox2=exp(-p2/correction);
   }
}

vpgslmat vtransi;
gsl_matrix * instrates;
gsl_vector * pnoe1out;
gsl_vector * pnoe2out;
gsl_matrix * id;
gsl_matrix * m1;
gsl_matrix * m2;
gsl_matrix * m3;
gsl_matrix * m4;
gsl_matrix * pij;
gsl_matrix * pijp;

gsl_vector * proba2;
int noemax;
vd dlca;//species distances to common last common ancestor

/*
   gsl_matrix *
   instant_rates_halpern (vd w, double time_step)
   {
   gsl_matrix * rates = gsl_matrix_calloc(4,4);
   pa=conca;
   pc=concc;

   double w0=w[0];
   double w1=w[1];
   double w2=w[2];
   double w3=w[3];
   if (w0<0 || w1<0 || w2<0 || w3<0){
   cout << "Problem in instant_rates_halpern" << endl;
   exit(1);
   }

   double fat=proba_fixation_rel(w1/w0);
   double fac=proba_fixation_rel(pa*w2/pc/w0);
   double fag=proba_fixation_rel(pa*w3/pc/w0);

   double fta=proba_fixation_rel(w0/w1);
   double ftc=proba_fixation_rel(pa*w2/pc/w1);
   double ftg=proba_fixation_rel(pa*w3/pc/w1);

   double fca=proba_fixation_rel(pc*w0/pa/w2);
   double fct=proba_fixation_rel(pc*w1/pa/w2);
   double fcg=proba_fixation_rel(w3/w2);

   double fga=proba_fixation_rel(pc*w0/pa/w3);
   double fgt=proba_fixation_rel(pc*w1/pa/w3);
   double fgc=proba_fixation_rel(w2/w3);

   double prefact=1.0/(4*kappa*pa*pc+0.5);

   gsl_matrix_set(m1, 0, 0, -(pa*fat+pc*fac+pc*kappa*fag));
   gsl_matrix_set(m1, 0, 1, pa*fta);
   gsl_matrix_set(m1, 0, 2, pa*fca);
   gsl_matrix_set(m1, 0, 3, pa*kappa*fga);

   gsl_matrix_set(m1, 1, 0, pa*fat);
   gsl_matrix_set(m1, 1, 1, -(pa*fta+pc*kappa*ftc+pc*ftg));
   gsl_matrix_set(m1, 1, 2, pa*kappa*fct);
   gsl_matrix_set(m1, 1, 3, pa*fgt);

   gsl_matrix_set(m1, 2, 0, pc*fac);
   gsl_matrix_set(m1, 2, 1, pc*kappa*ftc);
   gsl_matrix_set(m1, 2, 2, -(pa*fca+pa*kappa*fct+pc*fcg));
   gsl_matrix_set(m1, 2, 3, pc*fgc);

   gsl_matrix_set(m1, 3, 0, pc*kappa*fag);
   gsl_matrix_set(m1, 3, 1, pc*ftg);
   gsl_matrix_set(m1, 3, 2, pc*fcg);
   gsl_matrix_set(m1, 3, 3, -(pa*kappa*fga+pa*fgt+pc*fgc));

   gsl_matrix_scale(m1,prefact*time_step);
   gsl_matrix * mattemp;

   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,0.5,m1,m1,0.0,m2);
   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1/3.0,m1,m2,0.0,m3);
   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1/4.0,m1,m3,0.0,m4);

//Matrice d'évolution Runge Kutta 4
// Mat(RG4) = Id + h*M + h^2/2*M^2 + h^3/6*M^3 + h^4/24*M^4

//gsl_linalg_exponential_ss(m1, rates, 0.01);

gsl_matrix_memcpy(rates,id);
gsl_matrix_add(rates,m1);
gsl_matrix_add(rates,m2);
gsl_matrix_add(rates,m3);
gsl_matrix_add(rates,m4);


return rates;
}
*/

// This function returns M in
// P(t)=M(t)*P(0)
// where M is exp(R*t), with R the rate matrix
   vvd 
instant_rates_halpern (vd & w, double dist)
{
   gsl_matrix * rates = gsl_matrix_calloc(4,4);
   pa=conca;
   pc=concc;

   double w0=w[0];
   double w1=w[1];
   double w2=w[2];
   double w3=w[3];
   if (w0<0 || w1<0 || w2<0 || w3<0){
      cout << "Problem in instant_rates_halpern" << endl;
      exit(1);
   }

   double fat=proba_fixation_rel(w1/w0);
   double fac=proba_fixation_rel(pa*w2/pc/w0);
   double fag=proba_fixation_rel(pa*w3/pc/w0);

   double fta=proba_fixation_rel(w0/w1);
   double ftc=proba_fixation_rel(pa*w2/pc/w1);
   double ftg=proba_fixation_rel(pa*w3/pc/w1);

   double fca=proba_fixation_rel(pc*w0/pa/w2);
   double fct=proba_fixation_rel(pc*w1/pa/w2);
   double fcg=proba_fixation_rel(w3/w2);

   double fga=proba_fixation_rel(pc*w0/pa/w3);
   double fgt=proba_fixation_rel(pc*w1/pa/w3);
   double fgc=proba_fixation_rel(w2/w3);

   double prefact=1.0/(4*kappa*pa*pc+0.5);

   gsl_matrix_set(m1, 0, 0, -(pa*fat+pc*fac+pc*kappa*fag));
   gsl_matrix_set(m1, 0, 1, pa*fta);
   gsl_matrix_set(m1, 0, 2, pa*fca);
   gsl_matrix_set(m1, 0, 3, pa*kappa*fga);

   gsl_matrix_set(m1, 1, 0, pa*fat);
   gsl_matrix_set(m1, 1, 1, -(pa*fta+pc*kappa*ftc+pc*ftg));
   gsl_matrix_set(m1, 1, 2, pa*kappa*fct);
   gsl_matrix_set(m1, 1, 3, pa*fgt);

   gsl_matrix_set(m1, 2, 0, pc*fac);
   gsl_matrix_set(m1, 2, 1, pc*kappa*ftc);
   gsl_matrix_set(m1, 2, 2, -(pa*fca+pa*kappa*fct+pc*fcg));
   gsl_matrix_set(m1, 2, 3, pc*fgc);

   gsl_matrix_set(m1, 3, 0, pc*kappa*fag);
   gsl_matrix_set(m1, 3, 1, pc*ftg);
   gsl_matrix_set(m1, 3, 2, pc*fcg);
   gsl_matrix_set(m1, 3, 3, -(pa*kappa*fga+pa*fgt+pc*fgc));

   gsl_matrix_scale(m1,prefact*dist);
   gsl_linalg_exponential_ss(m1, rates, 1e-3);

   //Matrice d'évolution Runge Kutta 4
   // Mat(RG4) = Id + h*M + h^2/2*M^2 + h^3/6*M^3 + h^4/24*M^4

   //   gsl_matrix * mattemp;
   //
   //   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,0.5,m1,m1,0.0,m2);
   //   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1/3.0,m1,m2,0.0,m3);
   //   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1/4.0,m1,m3,0.0,m4);



   // rates=exp(m1), with precision of 1e-3
   //gsl_matrix_memcpy(rates,id);
   //gsl_matrix_add(rates,m1);
   //gsl_matrix_add(rates,m2);
   //gsl_matrix_add(rates,m3);
   //gsl_matrix_add(rates,m4);

   vd dum(4,0.);
   vvd r(4,dum);    
   for (unsigned int col=0;col<4;col++){
      for (unsigned int row=0;row<4;row++){
         r[col][row]=gsl_matrix_get(rates,row,col);
      }
   }

   return r;
}

   vvd 
instant_rates_felsen(vd & w, double dist)
{

   vd dum(4,0.);
   vvd M(4,dum);

   double correction=0.5+4.0*conca*concc;
   double prox=exp(-dist/correction);

   for (unsigned int col=0;col<4;col++){
      for (unsigned int row=0;row<4;row++){
         if (row==col) M[col][row]=prox+(1-prox)*w[row];
         else M[col][row]=(1.-prox)*w[row];
      }
   }

   return M;

}
//probs is the initial base prob in the tree, freq its selection frequency, dist

   vd 
evolvedist_felsen(vd & probs,vd & freqs, double dist)
{

   vd pf(4,0.);
   vvd M=instant_rates_felsen(freqs,dist);
   for (unsigned int row=0;row<4;row++){
      for (unsigned int k=0;k<4;k++){
         pf[row]+=M[k][row]*probs[k];
      }
   }
   return pf;
}
   
   vd 
evolvedist_felsen_backwards(vd & probs,vd & freqs, double dist)
{

   vd pf(4,0.);
   vvd M=instant_rates_felsen(freqs,dist);
   double sum=0;
   for (unsigned int row=0;row<4;row++){
      for (unsigned int k=0;k<4;k++){
         pf[row]+=M[row][k]*probs[k];
      }
      sum+=pf[row];
   }
   for (unsigned int row=0;row<4;row++){
      pf[row]/=sum;
   }
   return pf;
}

   vd 
evolvedist_halpern(vd probs,vd freqs, double dist)
{

   vd pf(4,0.);
   vvd M=instant_rates_halpern(freqs,dist);
   for (unsigned int row=0;row<4;row++){
      for (unsigned int k=0;k<4;k++){
         pf[row]+=M[k][row]*probs[k];
      }
   }
   return pf;
}

   vd 
evolvedist(vd probs,vd freqs, double dist)
{
   if (evolutionary_model==1) return evolvedist_felsen(probs,freqs,dist);
   else if (evolutionary_model==2) return evolvedist_halpern(probs,freqs,dist);
}

   void
evolvedisttest()
{
   vd freqs;
   freqs.push_back(conca);
   freqs.push_back(conct);
   freqs.push_back(concc);
   freqs.push_back(concg);
   for (int base=0;base<4;base++){
      ofstream outf("evol_cg.dat");
      vd initprob(4,0);
      initprob[2]=1.;
      for (double dist=0;dist<5.;dist+=0.01){
         vd pdisp;
         pdisp=evolvedist(initprob,freqs,dist);
         outf << dist << " ";
         for (ivd iv=pdisp.begin();iv!=pdisp.end();iv++){
            outf << *iv << " ";
         }
         outf << endl;
      }
      outf.close();
   }
}

// recursive function
// 
   void
loopsites(int n,vint& site,vvint& sites)
{
   if (n==0)
   {
      sites.push_back(site);
   }
   else
   {
      for(int j=0;j<4;j++)
      {
         vint sitetemp(site);
         if (n==10)
         {
            sitetemp.clear();
            sitetemp.push_back(j);
         }
         else
         {
            sitetemp.push_back(j);
         }
         loopsites(n-1,sitetemp,sites);
      }
   }
   return;
}

void evolvemean(Motif & mot)
{
   //generates initital site
   vint initsite;
   vvd initprob;
   vint sitegiven;
   sitegiven.push_back(2);
   sitegiven.push_back(1);
   sitegiven.push_back(3);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(1);
   sitegiven.push_back(1);//!!!!delete it for mus
   vd vpback;
   vpback.push_back(conca);
   vpback.push_back(conct);
   vpback.push_back(concc);
   vpback.push_back(concg);
   //   for (ivvd im=mot.matfreq.begin();im!=mot.matfreq.end();im++){
   //      vd probcol(4,0.);
   //      unsigned int base=gsl_ran_discrete (gslran,gslunif);
   //      probcol[base]=1.;
   //      initprob.push_back(probcol);
   //      initsite.push_back(base);
   //   }
   for (ivint is=sitegiven.begin();is!=sitegiven.end();is++){
      vd probcol(4,0.);
      unsigned int base=(*is);
      probcol[base]=1.;
      initprob.push_back(probcol);
      initsite.push_back(base);
   }

   vint::const_iterator ivc=initsite.begin();
   double ic1;
   ic1=scoref(ivc,mot.matprec);
   ivc=initsite.begin();
   double icr1(scoref(ivc,mot.matprecrevcomp));
   if (icr1>ic1) ic1=icr1;
   cout << vinttostring(initsite) << ", score= " << ic1 << endl;

   //generates the set of sites for ic distribution
   vvint sites;
   vint site;
   loopsites(width,site,sites);
   int numbersites(sites.size());

   vd means; 
   vd derivs;
   vd vars;
   double meantemp(ic1);
   double disttemp(0.);
   double distances[]={0.001,0.01,0.1,0.2,0.3,0.5,0.6,0.7,0.8,0.9,1,2,3,5,10};
   int numdist=sizeof(distances)/sizeof(double);

   for (int id=0;id<numdist;id++){

      double dist;
      dist=distances[id];
      int col(0);
      vvd probdist;

      //generates site probability at distance dist
      for (ivvd im=mot.matfreq.begin();im!=mot.matfreq.end();im++){
         vd probs=evolvedist(initprob[col],*im,dist);
         //vd probs=evolvedist(initprob[col],vpback,dist);
         probdist.push_back(probs);
         col++;
      }

      //computes mean of ic distribution
      double mean=0.;
      double var=0.;
      double deriv=0.;
      for (ivvint ivv=sites.begin();ivv!=sites.end();ivv++){
         ivc=(*ivv).begin();
         double ic;
         ic=scoref(ivc,mot.matprec);
         ivc=(*ivv).begin();
         double icr(scoref(ivc,mot.matprecrevcomp));
         if (icr>ic) ic=icr;

         double product(1);
         col=0;
         for (ivint ivi=(*ivv).begin();ivi!=(*ivv).end();ivi++){
            product*=probdist[col][(*ivi)];
            col++;
         }
         mean+=ic*product;
         var+=product*pow(ic,2);
      }
      var-=pow(mean,2);

      deriv=(mean-meantemp)/(dist-disttemp);

      cout << dist << " " << mean << " " << var <<  " " << deriv << endl;

      derivs.push_back(deriv);
      means.push_back(mean);
      vars.push_back(var);
      meantemp=mean;
      disttemp=dist;
   }


   system("if ! test -d icdist;then mkdir icdist;fi;");      
   ofstream outf("icdist/icmean.dat");
   for (int j=0;j<numdist;j++){
      outf << distances[j] << "\t" << means[j] << "\t" << vars[j] << "\t" << derivs[j] << "\n";
   }
   outf.close();

   return;
}

   void 
evolvebarrier(Motif & mot)
{

   double sth=-10.;
   //generates initital site
   vint initsite;
   vvd initprob;
   vint sitegiven;
   sitegiven.push_back(2);
   sitegiven.push_back(1);
   sitegiven.push_back(3);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(1);
   vd vpback;
   vpback.push_back(conca);
   vpback.push_back(conct);
   vpback.push_back(concc);
   vpback.push_back(concg);

   for (ivint is=sitegiven.begin();is!=sitegiven.end();is++){
      vd probcol(4,0.);
      unsigned int base=(*is);
      probcol[base]=1.;
      initprob.push_back(probcol);
      initsite.push_back(base);
   }

   int numbersites(300000);
   double binsize;
   binsize=0.05;
   double minbin(-50.);
   double maxbin(10.);
   int numberbins(floor((maxbin-minbin)/binsize));
   vvd probsinfo;

   vvd distinfo;
   double distances[]={0.001,0.1,1,10};
   int numdist=sizeof(distances)/sizeof(double);
   for (int id=0;id<numdist;id++){
      double dist;
      dist=distances[id];
      int col(0);
      vvd probdist;
      //generates site probability at distance dist
      for (ivvd im=mot.matfreq.begin();im!=mot.matfreq.end();im++){
         vd probs=evolvedist(initprob[col],*im,dist);
         //vd probs=evolvedist(initprob[col],vpback,dist);
         probdist.push_back(probs);
         col++;
      }
      //generates a set of sites for ic distribution
      vvint sites;
      for (int i=0;i<numbersites;i++){
         vint site;
         for (int j=0;j<width;j++){
            double pdist[4]={probdist[j][0],probdist[j][1],probdist[j][2],probdist[j][3]};
            gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
            unsigned int base=gsl_ran_discrete (gslran,gsldist);
            gsl_ran_discrete_free (gsldist);
            site.push_back(base);
         }
         sites.push_back(site);
      }
      //generates ic distribution at distance dist
      vd icdist;
      for (ivvint ivv=sites.begin();ivv!=sites.end();ivv++){
         vint::const_iterator ivc=(*ivv).begin();
         double ic;
         ic=scoref(ivc,mot.matprec);
         double icr(scoref(ivc,mot.matprecrevcomp));
         if (icr>ic) ic=icr;
         icdist.push_back(ic);
      }
      sort(icdist.begin(),icdist.end());
      distinfo.push_back(icdist);

      vd histinfo(numberbins,0.);
      for (ivd icd=icdist.begin();icd!=icdist.end();icd++){
         double ictemp;
         if ((*icd)<minbin) ictemp=minbin;
         else if ((*icd)>maxbin) ictemp=maxbin;
         else ictemp=*icd;
         int index;
         index=floor((ictemp-minbin)/binsize);
         histinfo[index]+=(double)1/numbersites;
         //histinfo[index]++;
      }
      probsinfo.push_back(histinfo);
   }

   vint::const_iterator ivc=initsite.begin();
   double ic1;
   ic1=scoref(ivc,mot.matprec);
   double icr1(scoref(ivc,mot.matprecrevcomp));
   if (icr1>ic1) ic1=icr1;
   cout << vinttostring(initsite) << ", score= " << ic1 << endl;


   system("if ! test -d icdist;then mkdir icdist;fi;");      
   ofstream outf("icdist/icprobs.dat");
   for (int i=0;i<numberbins;i++){
      double scorebin;
      scorebin=(binsize*i+minbin)+binsize/2;
      outf << scorebin << "\t";
      for (int j=0;j<numdist;j++){
         outf << probsinfo[j][i] << "\t";
      }
      outf << endl;
   }
   outf.close();
   outf.open("icdist/icdist.dat");
   for (int i=0;i<numbersites;i++){
      for (int j=0;j<numdist;j++){
         outf << distinfo[j][i] << "\t";
      }
      outf << endl;
   }
   outf.close();

   return;


}

   void
evolvebystep(Motif & mot)
{
   double sth=4.;
   //generates initital site
   vint initsite;
   vvd initprob;
   vint sitegiven;
   sitegiven.push_back(0);
   sitegiven.push_back(2);
   sitegiven.push_back(1);
   sitegiven.push_back(3);
   sitegiven.push_back(1);
   sitegiven.push_back(1);
   sitegiven.push_back(0);
   sitegiven.push_back(3);
   sitegiven.push_back(3);
   sitegiven.push_back(0);//!!!!delete it for mus
   vd vpback;
   vpback.push_back(conca);
   vpback.push_back(conct);
   vpback.push_back(concc);
   vpback.push_back(concg);
   //   for (ivvd im=mot.matfreq.begin();im!=mot.matfreq.end();im++){
   //      vd probcol(4,0.);
   //      double pmat[4]={(*im)[0],(*im)[1],(*im)[2],(*im)[3]};
   //      gsl_ran_discrete_t * gslmat=gsl_ran_discrete_preproc (4,pmat);
   //      unsigned int base=gsl_ran_discrete (gslran,gslmat);
   //      gsl_ran_discrete_free (gslmat);
   //      probcol[base]=1.;
   //      initprob.push_back(probcol);
   //      initsite.push_back(base);
   //   }
   for (ivint is=sitegiven.begin();is!=sitegiven.end();is++){
      vd probcol(4,0.);
      unsigned int base=(*is);
      probcol[base]=1.;
      initprob.push_back(probcol);
      initsite.push_back(base);
   }

   vint::const_iterator ivc=initsite.begin();
   double ic1(scoref(ivc,mot.matprec));
   ivc=initsite.begin();
   double icr1(scoref(ivc,mot.matprecrevcomp));
   if (icr1>ic1) ic1=icr1;
   cout << "Init: " << vinttostring(initsite) << ", score= " << ic1 << endl;

   int numbermut(100);
   double binsize;
   binsize=0.05;
   double minbin(-40.);
   double maxbin(15.);
   int numberbins(floor((maxbin-minbin)/binsize));
   vvd probsinfo;

   ofstream outf("icdist/icstep.dat");

   for (int j=0;j<10000;j++){
      vint testsite=initsite;
      double ic,ictemp;
      for (int i=0;i<numbermut;i++)
      {
         double psite[width];
         for (int j=0;j<10;j++){
            psite[j]=(double)1/width;
         }
         gsl_ran_discrete_t * gslsite=gsl_ran_discrete_preproc (10,psite);
         unsigned int numbase=gsl_ran_discrete (gslran,gslsite);
         gsl_ran_discrete_free (gslsite);

         unsigned int oldbase;
         oldbase=testsite[numbase];

         double pbase[4];
         for (int j=0;j<4;j++){
            pbase[j]=mot.matfreq[numbase][j];
         }
         //      vd thisbase(4,0.);
         //      thisbase[oldbase]=1.;
         //      vd vpbase=evolvedist(thisbase,mot.matfreq[numbase],0.01);
         //      double pbase[4]={vpbase[0],vpbase[1],vpbase[2],vpbase[3]};

         gsl_ran_discrete_t * gslbase=gsl_ran_discrete_preproc (4,pbase);
         unsigned int newbase=gsl_ran_discrete (gslran,gslbase);
         gsl_ran_discrete_free (gslbase);

         //cout << "base number " << numbase << " was " << oldbase << " and is now " << newbase  << endl;   

         initsite[numbase]=newbase;

         ivc=testsite.begin();
         ic=scoref(ivc,mot.matprec);
         ivc=testsite.begin();
         double icr(scoref(ivc,mot.matprecrevcomp));
         if (icr>ic) ic=icr;

         if (ic<sth){
            testsite[numbase]=oldbase;
            ic=ictemp;
            //cout << "IC: " << ic << endl;
         }
         ictemp=ic;
         //cout << "Evol " << i <<  ": " << vinttostring(initsite) << ", score= " << ic << endl;
      }
      outf << ic << endl;
   }
   outf.close();
   exit(9);

   //      //generates ic distribution at distance dist
   //      vd icdist;
   //      for (ivvint ivv=sites.begin();ivv!=sites.end();ivv++){
   //         vint::const_iterator ivc=(*ivv).begin();
   //         double ic;
   //         ic=scoref(ivc,mot.matprec);
   //          ivc=(*ivv).begin();
   //         double icr(scoref(ivc,mot.matprecrevcomp));
   //         if (icr>ic) ic=icr;
   //         icdist.push_back(ic);
   //      }
   //      sort(icdist.begin(),icdist.end());
   //      distinfo.push_back(icdist);
   //
   //      vd histinfo(numberbins,0.);
   //      for (ivd icd=icdist.begin();icd!=icdist.end();icd++){
   //         double ictemp;
   //         if ((*icd)<minbin) ictemp=minbin;
   //         else if ((*icd)>maxbin) ictemp=maxbin;
   //         else ictemp=*icd;
   //         int index;
   //         index=floor((ictemp-minbin)/binsize);
   //         histinfo[index]+=(double)1/numbersites;
   //         //histinfo[index]++;
   //      }
   //      probsinfo.push_back(histinfo);
   //   }
   //
   //
   //
   //   system("if ! test -d icdist;then mkdir icdist;fi;");      
   //   ofstream outf("icdist/icprobs.dat");
   //   for (int i=0;i<numberbins;i++){
   //      double scorebin;
   //      scorebin=(binsize*i+minbin)+binsize/2;
   //      outf << scorebin << "\t";
   //      for (int j=0;j<numdist;j++){
   //         outf << probsinfo[j][i] << "\t";
   //      }
   //      outf << endl;
   //   }
   //   outf.close();
   //   outf.open("icdist/icdist.dat");
   //   for (int i=0;i<numbersites;i++){
   //      for (int j=0;j<numdist;j++){
   //         outf << distinfo[j][i] << "\t";
   //      }
   //      outf << endl;
   //   }
   //   outf.close();

   return;

   }

bool inversesort (double d1,double d2) { return (d1>d2); }
double factorial(int n)
{
   double t=1;
   for (int i=n;i>1;i--){
      t*=i;
   }
   return t;
}
bool operator<(const icprob & info1,const icprob & info2)
{
   return info1.ic < info2.ic; 
}

void evolveic(Motif & mot)
{
   //generates initital site
   vint initsite;
   vvd initprob;
   vint sitegiven;
   sitegiven.push_back(2);
   sitegiven.push_back(1);
   sitegiven.push_back(3);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(1);
   sitegiven.push_back(1);//!!!!delete it for mus
   vd vpback;
   vpback.push_back(conca);
   vpback.push_back(conct);
   vpback.push_back(concc);
   vpback.push_back(concg);
   //   for (ivvd im=mot.matfreq.begin();im!=mot.matfreq.end();im++){
   //      vd probcol(4,0.);
   //      unsigned int base=gsl_ran_discrete (gslran,gslunif);
   //      probcol[base]=1.;
   //      initprob.push_back(probcol);
   //      initsite.push_back(base);
   //   }


   //Whole set of sites for ic distribution
   vvint sites;
   vint site;
   loopsites(width,site,sites);
   int numbersites(sites.size());

   //initial site probabilities
   for (ivint is=sitegiven.begin();is!=sitegiven.end();is++){
      vd probcol(4,0.);
      unsigned int base=(*is);
      probcol[base]=1.;
      initprob.push_back(probcol);
      initsite.push_back(base);
   }

   //Initial site score
   vint::const_iterator ivc=initsite.begin();
   double ic1;
   ic1=scoref(ivc,mot.matenergy);
   //   ic1=scoref(ivc,mot.matprec);
   //   ivc=initsite.begin();
   //   double icr1(scoref(ivc,mot.matprecrevcomp));
   //   if (icr1>ic1) ic1=icr1;
   cout << vinttostring(initsite) << ", score= " << ic1 << endl;

   double binsize;
   binsize=0.1;
   double minbin(0.);
   double maxbin(60.);
   //   binsize=0.5;
   //   double minbin(-40.);
   //   double maxbin(15.);
   int numberbins(floor((maxbin-minbin)/binsize));
   vvd probsinfo;
   vvic distinfo;

   double distances[]={0.01,0.1,1,10};
   int numdist=sizeof(distances)/sizeof(double);

   for (int id=0;id<numdist;id++){

      double dist;
      dist=distances[id];
      vvd probdist;
      int col(0);
      //generates site probability at distance dist
      for (ivvd im=mot.matfreq.begin();im!=mot.matfreq.end();im++){
         //vd probs=evolvedist(initprob[col],*im,dist);
         vd probs=evolvedist(initprob[col],vpback,dist);
         probdist.push_back(probs);
         col++;
      }

      //generates ic distribution at distance dist
      vic vinfo;
      for (ivvint ivv=sites.begin();ivv!=sites.end();ivv++){
         icprob info;

         vint::const_iterator ivc=(*ivv).begin();
         double ic;
         ic=scoref(ivc,mot.matenergy);
         //         ic=scoref(ivc,mot.matprec);
         //         ivc=(*ivv).begin();
         //         double icr(scoref(ivc,mot.matprecrevcomp));
         //         if (icr>ic) ic=icr;

         info.ic=ic;

         double weight=1;
         col=0;
         for (ivint iv=(*ivv).begin();iv!=(*ivv).end();iv++){
            weight*=probdist[col][*iv];
            col++;
         }
         info.prob=weight;
         vinfo.push_back(info);
      }
      sort(vinfo.begin(),vinfo.end());
      distinfo.push_back(vinfo);

      vd histinfo(numberbins,0.);
      vint histnumperbin(numberbins,0);
      for (ivic icd=vinfo.begin();icd!=vinfo.end();icd++){
         double icdi=(*icd).ic;
         double ictemp;
         if (icdi<minbin) ictemp=minbin;
         else if (icdi>maxbin) ictemp=maxbin;
         else ictemp=icdi;
         int index;
         index=floor((ictemp-minbin)/binsize);
         histinfo[index]+=(*icd).prob;///binsize;
         //histnumperbin[index]++;
      }
      //      double histarea(0);
      //      for (int ihist=0;ihist<numberbins;ihist++){
      //         histarea+=histinfo[ihist];
      //      }
      //      for (int ihist=0;ihist<numberbins;ihist++){
      //         histinfo[ihist]/=histarea;
      //      }
      probsinfo.push_back(histinfo);
   }


   system("if ! test -d icdist;then mkdir icdist;fi;");      
   ofstream outf("icdist/icprobs.dat");
   for (int i=0;i<numberbins;i++){
      double scorebin;
      scorebin=(binsize*i+minbin)+binsize/2;
      outf << scorebin << "\t";
      for (int j=0;j<numdist;j++){
         outf << probsinfo[j][i] << "\t";
      }
      outf << endl;
   }
   //   for (int i=0;i<numbersites;i++){
   //      outf << distinfo[0][i].ic << "\t";
   //      for (int j=0;j<numdist;j++){
   //         outf << distinfo[j][i].prob << "\t";
   //      }
   //      outf << endl;
   //   }
   outf.close();
   //   ofstream outf("icdist/icprobs.dat");
   //   for (int i=0;i<numberbins;i++){
   //      double scorebin;
   //      scorebin=(binsize*i+minbin)+binsize/2;
   //      outf << scorebin << "\t";
   //      for (int j=0;j<numdist;j++){
   //         outf << probsinfo[j][i] << "\t";
   //      }
   //      outf << endl;
   //   }
   //   outf.close();
   //   outf.open("icdist/icdist.dat");
   //   for (int i=0;i<numbersites;i++){
   //      for (int j=0;j<numdist;j++){
   //         outf << distinfo[j][i] << "\t";
   //      }
   //      outf << endl;
   //   }
   //   outf.close();

   return;

}

/*
   void evolveicapprox(Motif & mot)
   {
//generates initital site
vint initsite;
vvd initprob;
vint sitegiven;
sitegiven.push_back(2);
sitegiven.push_back(1);
sitegiven.push_back(3);
sitegiven.push_back(2);
sitegiven.push_back(2);
sitegiven.push_back(2);
sitegiven.push_back(2);
sitegiven.push_back(2);
sitegiven.push_back(1);
sitegiven.push_back(1);//!!!!delete it for mus
vd vpback;
vpback.push_back(conca);
vpback.push_back(conct);
vpback.push_back(concc);
vpback.push_back(concg);
//   for (ivvd im=mot.matfreq.begin();im!=mot.matfreq.end();im++){
//      vd probcol(4,0.);
//      unsigned int base=gsl_ran_discrete (gslran,gslunif);
//      probcol[base]=1.;
//      initprob.push_back(probcol);
//      initsite.push_back(base);
//   }

//draw a random base in the site
double psite[width];
for (int j=0;j<width;j++){
psite[j]=(double)1/width;
}
gsl_ran_discrete_t * gslsite=gsl_ran_discrete_preproc (width,psite);


//Whole set of sites for ic distribution
vvint sites;
vint site;
loopsites(width,site,sites);
int numbersites(sites.size());

//initial site probabilities
for (ivint is=sitegiven.begin();is!=sitegiven.end();is++){
vd probcol(4,0.);
unsigned int base=(*is);
probcol[base]=1.;
initprob.push_back(probcol);
initsite.push_back(base);
}

//Initial site score
vint::const_iterator ivc=initsite.begin();
double ic1;
ic1=scoref(ivc,mot.matprec);
ivc=initsite.begin();
double icr1(scoref(ivc,mot.matprecrevcomp));
if (icr1>ic1) ic1=icr1;
cout << vinttostring(initsite) << ", score= " << ic1 << endl;

double binsize;
binsize=0.05;
double minbin(-40.);
double maxbin(15.);
int numberbins(floor((maxbin-minbin)/binsize));
vvd probsinfo;

vvd distinfo;
double distances[]={0.001,0.1,1,10};
int numdist=sizeof(distances)/sizeof(double);

system("if ! test -d icdist;then mkdir icdist;fi;");      
system("if ! test -d icdist/rankprobs;then mkdir icdist/rankprobs;fi;");      
int fcount(1);
for (int id=0;id<numdist;id++){

   ostringstream os;
   os << "icdist/rankprobs/rankprobs_" << fcount << ".dat";
   ofstream outf(os.str().c_str());
   double dist;
   dist=distances[id];

   vvd probdist;
   vvd rankmat;//each element is a base number (row=rank-1,column=site position)
   int col(0);
   //generates site probability at distance dist
   for (ivvd im=mot.matfreq.begin();im!=mot.matfreq.end();im++){
      vd probs=evolvedist(initprob[col],*im,dist);
      //vd probs=evolvedist(initprob[col],vpback,dist);
      probdist.push_back(probs);
      vd probs2=probs;
      sort(probs2.begin(),probs2.end(),inversesort);
      vd rankcol;
      for (ivd ivp2=probs2.begin();ivp2!=probs2.end();ivp2++){
         int ipos(0);
         for (ivd ivp=probs.begin();ivp!=probs.end();ivp++){
            if ((*ivp2)==(*ivp)){
               *ivp=-1;
               break;
            }
            ipos++;
         }
         rankcol.push_back(ipos);
      }
      rankmat.push_back(rankcol);//rankcol);
      col++;
   }

   //we draw the most interesting sites 
   vvint ranksites;

   vint bestsite;
   vint worstsite;
   for (ivvd ivv=rankmat.begin();ivv!=rankmat.end();ivv++){
      bestsite.push_back((*ivv)[0]);
      worstsite.push_back((*ivv)[3]);
   }
   ranksites.push_back(bestsite);
   ranksites.push_back(worstsite);
   vint::const_iterator ivc=bestsite.begin();
   double ic;
   ic=scoref(ivc,mot.matprec);
   ivc=bestsite.begin();
   double icr(scoref(ivc,mot.matprecrevcomp));
   if (icr>ic) ic=icr;
   // icdist.push_back(ic);

   double weight;
   weight=1;   
   for (int j=0;j<width;j++){
      weight*=probdist[j][bestsite[j]];
   }
   if (weight<1.) weight=0.;
   outf << ic << " " << weight << endl;

   ivc=worstsite.begin();
   ic=scoref(ivc,mot.matprec);
   ivc=bestsite.begin();
   icr=scoref(ivc,mot.matprecrevcomp);
   if (icr>ic) ic=icr;
   // icdist.push_back(ic);
   weight=1;   
   for (int j=0;j<width;j++){
      weight*=probdist[j][worstsite[j]];
   }
   if (weight<1.) weight=0.;
   outf << ic << " " << weight << endl;

   for (int i=0;i<width;i++){
      bestsite[i]=rankmat[i][1];
      worstsite[i]=rankmat[i][2];
      ranksites.push_back(bestsite);
      ranksites.push_back(worstsite);

      ivc=bestsite.begin();
      ic=scoref(ivc,mot.matprec);
      ivc=bestsite.begin();
      icr=scoref(ivc,mot.matprecrevcomp);
      if (icr>ic) ic=icr;
      // icdist.push_back(ic);
      double weight;
      weight=pow(3,i+1)*factorial(width)/factorial(width-i-1);   
      for (int j=0;j<width;j++){
         weight*=probdist[j][bestsite[j]];
      }
      if (weight<1.) weight=0.;
      outf << ic << " " << weight << endl;

      ivc=worstsite.begin();
      ic=scoref(ivc,mot.matprec);
      ivc=bestsite.begin();
      icr=scoref(ivc,mot.matprecrevcomp);
      if (icr>ic) ic=icr;
      // icdist.push_back(ic);
      weight=pow(3,i+1)*factorial(width)/factorial(width-i-1);   
      for (int j=0;j<width;j++){
         weight*=probdist[j][worstsite[j]];
      }
      if (weight<1.) weight=0.;
      outf << ic << " " << weight << endl;
   }
   //
   //      for (int i=0;i<width;i++){
   //         vint currbest=bestsite;
   //         vint currworst=worstsite;
   //         currbest[i]=rankmat[i][1];
   //         currworst[i]=rankmat[i][2];
   //         ranksites.push_back(currbest);
   //         ranksites.push_back(currworst);
   //         for (int j=0;j<width;j++){
   //            vint currbest2=currbest;
   //            vint currworst2=currworst;
   //            if (j!=i){
   //               currbest2[j]=rankmat[j][1];
   //               currworst2[j]=rankmat[j][2];
   //               ranksites.push_back(currbest2);
   //               ranksites.push_back(currworst2);
   //               for (int k=0;k<width;k++){
   //                  vint currbest3=currbest2;
   //                  vint currworst3=currworst2;
   //                  if (k!=j && k!=i){
   //                     currbest3[k]=rankmat[k][1];
   //                     currworst3[k]=rankmat[k][2];
   //                     ranksites.push_back(currbest3);
   //                     ranksites.push_back(currworst3);
   //                  }
   //               }
   //            }
   //         }
   //      }

   //generates a set of sites for ic distribution
   //      vvint sites;
   //      for (int i=0;i<numbersites;i++){
   //         vint site;
   //         for (int j=0;j<width;j++){
   //            double pdist[4]={probdist[j][0],probdist[j][1],probdist[j][2],probdist[j][3]};
   //            gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
   //            unsigned int base=gsl_ran_discrete (gslran,gsldist);
   //            gsl_ran_discrete_free (gsldist);
   //            site.push_back(base);
   //         }
   //         sites.push_back(site);
   //      }
   //generates ic distribution at distance dist
   // vd icdist;
   //      for (ivvint ivv=ranksites.begin();ivv!=ranksites.end();ivv++){
   //         vint::const_iterator ivc=(*ivv).begin();
   //         double ic;
   //         ic=scoref(ivc,mot.matprec);
   //         ivc=(*ivv).begin();
   //         double icr(scoref(ivc,mot.matprecrevcomp));
   //         if (icr>ic) ic=icr;
   //        // icdist.push_back(ic);
   //
   //         double weight;
   //         weight=pow(4,width);   
   //         for (int j=0;j<width;j++){
   //            weight*=probdist[j][(*ivv)[j]];
   //         }
   //         outf << ic << " " << weight << endl;
   //         //cout << currbest << " " << currworst << endl;exit(9);
   //      }
   outf.close();
   fcount++;


   //sort(icdist.begin(),icdist.end());
   //distinfo.push_back(icdist);

   //      vd histinfo(numberbins,0.);
   //      for (ivd icd=icdist.begin();icd!=icdist.end();icd++){
   //         double ictemp;
   //         if ((*icd)<minbin) ictemp=minbin;
   //         else if ((*icd)>maxbin) ictemp=maxbin;
   //         else ictemp=*icd;
   //         int index;
   //         index=floor((ictemp-minbin)/binsize);
   //         histinfo[index]+=(double)1/numbersites;
   //         //histinfo[index]++;
   //      }
   //      probsinfo.push_back(histinfo);

}


//   ofstream outf("icdist/icprobs.dat");
//   for (int i=0;i<numberbins;i++){
//      double scorebin;
//      scorebin=(binsize*i+minbin)+binsize/2;
//      outf << scorebin << "\t";
//      for (int j=0;j<numdist;j++){
//         outf << probsinfo[j][i] << "\t";
//      }
//      outf << endl;
//   }
//   outf.close();
//   outf.open("icdist/icdist.dat");
//   for (int i=0;i<numbersites;i++){
//      for (int j=0;j<numdist;j++){
//         outf << distinfo[j][i] << "\t";
//      }
//      outf << endl;
//   }
//   outf.close();

return;

}
*/

void evolveic2(Motif & mot)
{
   //generates initital site
   vint initsite;
   vvd initprob;
   vint sitegiven;
   sitegiven.push_back(2);
   sitegiven.push_back(1);
   sitegiven.push_back(3);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(2);
   sitegiven.push_back(1);
   sitegiven.push_back(1);//!!!!delete it for mus
   vd vpback;
   vpback.push_back(conca);
   vpback.push_back(conct);
   vpback.push_back(concc);
   vpback.push_back(concg);

   for (ivint is=sitegiven.begin();is!=sitegiven.end();is++){
      vd probcol(4,0.);
      unsigned int base=(*is);
      probcol[base]=1.;
      initprob.push_back(probcol);
      initsite.push_back(base);
   }

   vint::const_iterator ivc=initsite.begin();
   double ic0;
   ic0=scoref(ivc,mot.matprec);
   ivc=initsite.begin();
   double icr0(scoref(ivc,mot.matprecrevcomp));
   if (icr0>ic0) ic0=icr0;
   cout << vinttostring(initsite) << ", score= " << ic0 << endl;

   //generates the set of sites for ic distribution
   vvint sites;
   vint site;
   loopsites(width,site,sites);
   int numbersites(sites.size());

   double binsize;
   binsize=0.05;
   double minbin(-25.);
   double maxbin(10.);
   int numberbins(floor((maxbin-minbin)/binsize));

   vvd probsinfo;
   vvd distinfo;
   double dist0=1;
   double dist1=1;
   double dist2=1;



   int col(0);
   vvd halfprobs;
   vd halfweights;
   //generates site probability at dist0
   for (ivvd im=mot.matfreq.begin();im!=mot.matfreq.end();im++){
      vd halfprob=evolvedist(initprob[col],*im,dist0);
      col++;
      halfprobs.push_back(halfprob);
   }


   ofstream outf("heatmap.dat"); 
   int icount(0);
   //generates ic distribution at distance dist
   for (ivvint ivv1=sites.begin();ivv1!=sites.end();ivv1++){//loop on s1
      for (ivvint ivv2=sites.begin();ivv2!=sites.end();ivv2++){//loop on s2
         double totweight=0;
         double ic1,ic2;
         for (ivvint ivv=sites.begin();ivv!=sites.end();ivv++){//loop on s
            double weight=1;
            vvd halfinitprobs;
            col=0;
            //this is p(s)
            for (ivint iv=(*ivv).begin();iv!=(*ivv).end();iv++){
               weight*=halfprobs[col][*iv];
               vd probcol(4,0.);
               probcol[*iv]=1.;
               halfinitprobs.push_back(probcol);
            }
            //generates final probabilities
            vvd pf1;
            vvd pf2;
            col=0;
            for (ivvd im=mot.matfreq.begin();im!=mot.matfreq.end();im++){
               vd p1=evolvedist(halfinitprobs[col],*im,dist1);
               vd p2=evolvedist(halfinitprobs[col],*im,dist2);
               col++;
               pf1.push_back(p1);
               pf2.push_back(p2);
            }

            double weight1=1;
            col=0;
            for (ivint iv1=(*ivv1).begin();iv1!=(*ivv1).end();iv1++){
               weight1*=pf1[col][*iv1];
            }
            double weight2=1;
            col=0;
            for (ivint iv2=(*ivv2).begin();iv2!=(*ivv2).end();iv2++){
               weight2*=pf2[col][*iv2];
            }
            totweight+=weight*weight1*weight2;

            vint::const_iterator ivc1=(*ivv1).begin();
            ic1=scoref(ivc1,mot.matprec);
            ivc1=(*ivv1).begin();
            double icr1(scoref(ivc1,mot.matprecrevcomp));
            if (icr1>ic1) ic1=icr1;

            vint::const_iterator ivc2=(*ivv2).begin();
            ic2=scoref(ivc2,mot.matprec);
            ivc2=(*ivv2).begin();
            double icr2(scoref(ivc2,mot.matprecrevcomp));
            if (icr2>ic2) ic2=icr2;
         }
         outf << ic1 << "\t" << ic2 << "\t" << totweight << endl; 
         cout << icount << " ";
         cout.flush();
         icount++;
      }
   }
   outf.close();

   //   for (ivvint ivv=sites.begin();ivv!=sites.end();ivv++){
   //      vint::const_iterator ivc=(*ivv).begin();
   //      double ic;
   //      ic=scoref(ivc,mot.matprec);
   //      ivc=(*ivv).begin();
   //      double icr(scoref(ivc,mot.matprecrevcomp));
   //      if (icr>ic) ic=icr;
   //      icdist.push_back(ic);
   //   }
   //   sort(icdist.begin(),icdist.end());
   //
   //   vd histinfo(numberbins,0.);
   //   for (ivd icd=icdist.begin();icd!=icdist.end();icd++){
   //      double ictemp;
   //      if ((*icd)<minbin) ictemp=minbin;
   //      else if ((*icd)>maxbin) ictemp=maxbin;
   //      else ictemp=*icd;
   //      int index;
   //      index=floor((ictemp-minbin)/binsize);
   //      histinfo[index]+=(double)1/numbersites;
   //      //histinfo[index]++;
   //   }
   //   probsinfo.push_back(histinfo);
   //
   //
   //vint::const_iterator ivc=initsite.begin();
   //double ic1;
   //ic1=scoref(ivc,mot.matprec);
   //      ivc=(*ivv).begin();
   //double icr1(scoref(ivc,mot.matprecrevcomp));
   //if (icr1>ic1) ic1=icr1;
   //cout << vinttostring(initsite) << ", score= " << ic1 << endl;
   //
   //
   //   system("if ! test -d icdist;then mkdir icdist;fi;");      
   //   ofstream outf("icdist/icprobs.dat");
   //   for (int i=0;i<numberbins;i++){
   //      double scorebin;
   //      scorebin=(binsize*i+minbin)+binsize/2;
   //      outf << scorebin << "\t";
   //      for (int j=0;j<numdist;j++){
   //         outf << probsinfo[j][i] << "\t";
   //      }
   //      outf << endl;
   //   }
   //   outf.close();
   //   outf.open("icdist/icdist.dat");
   //   for (int i=0;i<numbersites;i++){
   //      for (int j=0;j<numdist;j++){
   //         outf << distinfo[j][i] << "\t";
   //      }
   //      outf << endl;
   //   }
   //   outf.close();

   return;


}


   vint
evolve(unsigned int base)
{
   vint bases;
   inittreedist();
   gsl_matrix * probatree=gsl_matrix_alloc(4,noemax+1); 

   for (unsigned int i=0;i<nbspecies;i++){ 
      gsl_matrix_set(probatree,0,i,0);
      gsl_matrix_set(probatree,1,i,0);
      gsl_matrix_set(probatree,2,i,0);
      gsl_matrix_set(probatree,3,i,0);
   }
   gsl_matrix_set(probatree,base,noemax,1);
   cout << endl;

   for (ivnoe iv=treedist.end()-1;iv>=treedist.begin();iv--){
      int noe=iv->noe;
      int n1=iv->esp1;
      int n2=iv->esp2;
      double prox1=iv->prox1;
      double prox2=iv->prox2;
      //cout << n1 << " " << n2 << " " << noe << " " << prox1 << " " << prox2 << " " << endl;
      gsl_vector_view pnoe=gsl_matrix_column(probatree,noe);
      gsl_vector_view pnoe1=gsl_matrix_column(probatree,n1);
      gsl_vector_view pnoe2=gsl_matrix_column(probatree,n2);
      gsl_ran_discrete_t * g;
      double p[4];
      double p1[4];
      double p2[4];
      double sum(0);
      for (int i=0;i<4;i++){
         p[i]=gsl_vector_get(&pnoe.vector,i);
         sum+=p[i];
      }
      double conc[4];
      conc[0]=conca;
      conc[1]=conct;
      conc[2]=concc;
      conc[3]=concg;
      for (int i=0;i<4;i++){
         p1[i]=prox1*p[i]+(1-prox1)*sum*conc[i];
         p2[i]=prox2*p[i]+(1-prox2)*sum*conc[i];
         gsl_matrix_set(probatree,i,n1,p1[i]);
         gsl_matrix_set(probatree,i,n2,p2[i]);
         //         cout << p[i] << " " << p1[i] << " " << p2[i] << endl;
      }
      //      cout << endl;
      //      unsigned int n[4];
      //         unsigned int basetemp=gsl_ran_discrete (gslran,g);
      ////      double sum;
      //      gsl_blas_ddot(wfull,&pnoe1.vector,&sum);
      //      sum*=(1-prox1);
      //      gsl_vector_set_all(&pnoe.vector,sum);
      //      gsl_blas_daxpy(prox1,&pnoe1.vector,&pnoe.vector);
   }
   for (int i=0;i<4;i++){
      for (int j=0;j<nbspecies;j++){
         cout <<  gsl_matrix_get(probatree,i,j)  << "\t";
      }
      cout << gsl_matrix_get(probatree,i,noemax) << endl;
   }

   return bases;
}

// we have dP/Dt=M*P, with M the instant rates matrix
// this solves to P(t)=exp(M*t)*P0
// approximated ro P(t+dt)=(1+Mdt+...)*P(t) with RK4
//                             VVV
//                          instrates 
// so that P(t+n*dt)=(instrates^n)*P(t)
// first step: pij=id <=> instrates^0
// then pij=instrates^(n-1)
// vtransi gives the transition matrix at different time scales
// each time scale can be common to several branches
//

// Given a PWM and a starting site on Common Ancestor, returns a set of aligned sites
   vvint
evolve_forward(Motif & mot,vint & site)
{

   vint idum(nbspecies,0);
   vvint results(mot.motwidth,idum);

   for (unsigned int pos=0;pos<mot.motwidth;pos++)
   {
      gsl_vector *w;
      w = gsl_vector_alloc (4);
      gsl_vector_set (w, 0, mot.matfreq[pos][0]);
      gsl_vector_set (w, 1, mot.matfreq[pos][1]);
      gsl_vector_set (w, 2, mot.matfreq[pos][2]);
      gsl_vector_set (w, 3, mot.matfreq[pos][3]);
      gsl_matrix * pmattemp;

      // DRAW EVOLVED BASES FOLLOWING THE TREE

      // INIT
      double pinit[4];
      gsl_ran_discrete_t * gslinit;
      //      double pinit[4]= { gsl_vector_get(w,0),
      //         gsl_vector_get(w,1),
      //         gsl_vector_get(w,2),
      //         gsl_vector_get(w,3) };
      //      gsl_ran_discrete_t * gslinit=gsl_ran_discrete_preproc (4,pinit);
      //      unsigned int base=gsl_ran_discrete (gslran,gslinit);
      //      gsl_ran_discrete_free (gslinit);
      unsigned int base=site[pos];
      gsl_matrix * probatree=gsl_matrix_alloc(4,noemax+1); 
      gsl_matrix_set_zero(probatree);
      gsl_vector_view pnoeinit=gsl_matrix_column(probatree,noemax);
      gsl_vector_set_zero(&pnoeinit.vector);
      gsl_vector_set(&pnoeinit.vector,base,1);

      vint result(nbspecies,0);
      for (vnoe::reverse_iterator iv=treedist.rbegin();iv!=treedist.rend();++iv){

         int n1=iv->esp1;
         int n2=iv->esp2;
         if (evolutionary_model==2){

            if (instant_rates(w,instrates)) exit(1);
            //      for (unsigned int i=0;i<4;i++){
            //         for (unsigned int j=0;j<4;j++){
            //            cout << gsl_matrix_get(instrates,i,j) << " ";
            //         }
            //         cout << endl;
            //      }
            //      cout << endl;

            gsl_matrix_memcpy(pij,id);

            if (species=="droso"){
               for (unsigned int i=1;i<117;i++){
                  //integration
                  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
                  pmattemp=pij;
                  pij=pijp;
                  pijp=pmattemp;
                  //      gsl_matrix_memcpy(pij,pijp);
                  if (i==1){
                     gsl_matrix_memcpy(vtransi[10],pij);
                  } else if (i==2){
                     gsl_matrix_memcpy(vtransi[0],pij);
                     gsl_matrix_memcpy(vtransi[11],pij);
                  } else if (i==3){
                     gsl_matrix_memcpy(vtransi[1],pij);
                     gsl_matrix_memcpy(vtransi[3],pij);
                  } else if (i==4){
                     gsl_matrix_memcpy(vtransi[7],pij);
                  } else if (i==7){
                     gsl_matrix_memcpy(vtransi[2],pij);
                     gsl_matrix_memcpy(vtransi[6],pij);
                  } else if (i==10){
                     gsl_matrix_memcpy(vtransi[5],pij);
                     gsl_matrix_memcpy(vtransi[17],pij);
                  } else if (i==11){
                     gsl_matrix_memcpy(vtransi[4],pij);
                  } else if (i==14){
                     gsl_matrix_memcpy(vtransi[20],pij);
                  } else if (i==15){
                     gsl_matrix_memcpy(vtransi[21],pij);
                  } else if (i==35){
                     gsl_matrix_memcpy(vtransi[12],pij);
                     gsl_matrix_memcpy(vtransi[14],pij);
                  } else if (i==42){
                     gsl_matrix_memcpy(vtransi[19],pij);
                  } else if (i==46){
                     gsl_matrix_memcpy(vtransi[16],pij);
                  } else if (i==49){
                     gsl_matrix_memcpy(vtransi[15],pij);
                  } else if (i==58){
                     gsl_matrix_memcpy(vtransi[9],pij);
                  } else if (i==68){
                     gsl_matrix_memcpy(vtransi[13],pij);
                  } else if (i==82){
                     gsl_matrix_memcpy(vtransi[8],pij);
                  } else if (i==116){
                     gsl_matrix_memcpy(vtransi[18],pij);
                  }
               }
            }
            else if (species=="eutherian"){
               // approx, integr_step=0.01
               for (unsigned int i=1;i<26;i++){
                  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
                  //            cout << "t=" << i*integr_step << endl;
                  //            for (unsigned int row=0;row<4;row++){
                  //               for (unsigned int col=0;col<4;col++){
                  //                  cout << gsl_matrix_get(pij,row,col) << " ";
                  //               }
                  //               cout << endl;
                  //            }
                  //            cout << endl;
                  pmattemp=pij;
                  pij=pijp;
                  pijp=pmattemp;
                  //      gsl_matrix_memcpy(pij,pijp);
                  if (i==1){
                     gsl_matrix_memcpy(vtransi[15],pij);
                     gsl_matrix_memcpy(vtransi[2],pij);
                     gsl_matrix_memcpy(vtransi[3],pij);
                     gsl_matrix_memcpy(vtransi[5],pij);
                     gsl_matrix_memcpy(vtransi[7],pij);
                  } else if (i==2){
                     gsl_matrix_memcpy(vtransi[4],pij);
                     gsl_matrix_memcpy(vtransi[16],pij);
                  } else if (i==3){
                     gsl_matrix_memcpy(vtransi[17],pij);
                  } else if (i==6){
                     gsl_matrix_memcpy(vtransi[6],pij);
                  } else if (i==8){
                     gsl_matrix_memcpy(vtransi[0],pij);
                     gsl_matrix_memcpy(vtransi[10],pij);
                     gsl_matrix_memcpy(vtransi[11],pij);
                     gsl_matrix_memcpy(vtransi[13],pij);
                     gsl_matrix_memcpy(vtransi[1],pij);
                  } else if (i==11){
                     gsl_matrix_memcpy(vtransi[9],pij);
                     gsl_matrix_memcpy(vtransi[14],pij);
                  } else if (i==15){
                     gsl_matrix_memcpy(vtransi[12],pij);
                  } else if (i==25){
                     gsl_matrix_memcpy(vtransi[8],pij);
                  }
               }
            }

            gsl_vector_view pnoe=gsl_matrix_column(probatree,iv->noe);

            // NODE 1
            gsl_vector_view pnoe1=gsl_matrix_column(probatree,n1);
            gsl_blas_dgemv(CblasNoTrans,1.0,vtransi[2*(iv->noe-nbspecies)+1],&pnoe.vector,0.0,&pnoe1.vector);

            pinit[0]= gsl_matrix_get(probatree,0,n1);
            pinit[1]= gsl_matrix_get(probatree,1,n1);
            pinit[2]= gsl_matrix_get(probatree,2,n1);
            pinit[3]= gsl_matrix_get(probatree,3,n1);
            gslinit=gsl_ran_discrete_preproc (4,pinit);
            base=gsl_ran_discrete (gslran,gslinit);
            if (n1<nbspecies) result[n1]=base;

            gsl_vector_set_zero(&pnoe1.vector);
            gsl_vector_set(&pnoe1.vector,base,1);

            // NODE 2
            gsl_vector_view pnoe2=gsl_matrix_column(probatree,n2);
            gsl_blas_dgemv(CblasNoTrans,1.0,vtransi[2*(iv->noe-nbspecies)],&pnoe.vector,0.0,&pnoe2.vector);

            pinit[0]= gsl_matrix_get(probatree,0,n2);
            pinit[1]= gsl_matrix_get(probatree,1,n2);
            pinit[2]= gsl_matrix_get(probatree,2,n2);
            pinit[3]= gsl_matrix_get(probatree,3,n2);
            gslinit=gsl_ran_discrete_preproc (4,pinit);
            base=gsl_ran_discrete (gslran,gslinit);
            if (n2<nbspecies) result[n2]=base;

            gsl_vector_set_zero(&pnoe2.vector);
            gsl_vector_set(&pnoe2.vector,base,1);
         }
         else if (evolutionary_model==1){

            double prox1=iv->prox1;
            double prox2=iv->prox2;
            gsl_vector_view pnoe=gsl_matrix_column(probatree,iv->noe);
            gsl_vector_view pnoe1=gsl_matrix_column(probatree,n1);
            gsl_vector_view pnoe2=gsl_matrix_column(probatree,n2);

            // SPECIES 1
            //                     for (unsigned int j=0;j<4;j++){
            //                        cout << gsl_vector_get(&pnoe.vector,j) << " ";
            //                     }
            //                     cout << endl;
            //            cout << prox1 << endl;
            //                     for (unsigned int j=0;j<4;j++){
            //                        cout << gsl_vector_get(w,j) << " ";
            //                     }
            //                     cout << endl;
            //                     for (unsigned int j=0;j<4;j++){
            //                        cout << gsl_vector_get(&pnoe1.vector,j) << " ";
            //                     }
            //                     cout << endl;
            gsl_vector_memcpy(&pnoe1.vector,w);
            //                     for (unsigned int j=0;j<4;j++){
            //                        cout << gsl_vector_get(&pnoe1.vector,j) << " ";
            //                     }
            //                     cout << endl;
            gsl_blas_dscal(1-prox1,&pnoe1.vector);
            //                     for (unsigned int j=0;j<4;j++){
            //                        cout << gsl_vector_get(&pnoe1.vector,j) << " ";
            //                     }
            //                     cout << endl;
            gsl_blas_daxpy(prox1,&pnoe.vector,&pnoe1.vector);
            //                     for (unsigned int j=0;j<4;j++){
            //                        cout << gsl_vector_get(&pnoe1.vector,j) << " ";
            //                     }
            //                     cout << endl;

            pinit[0]= gsl_matrix_get(probatree,0,n1);
            pinit[1]= gsl_matrix_get(probatree,1,n1);
            pinit[2]= gsl_matrix_get(probatree,2,n1);
            pinit[3]= gsl_matrix_get(probatree,3,n1);
            //                        for (unsigned int col=0;col<4;col++){
            //                           cout << pinit[col] << " ";
            //                        }
            //                        cout << endl;
            gslinit=gsl_ran_discrete_preproc (4,pinit);
            base=gsl_ran_discrete (gslran,gslinit);
            if (n1<nbspecies) result[n1]=base;

            gsl_vector_set_zero(&pnoe1.vector);
            gsl_vector_set(&pnoe1.vector,base,1);


            // SPECIES 2   
            gsl_vector_memcpy(&pnoe2.vector,w);
            gsl_blas_dscal(1-prox2,&pnoe2.vector);
            gsl_blas_daxpy(prox2,&pnoe.vector,&pnoe2.vector);

            pinit[0]= gsl_matrix_get(probatree,0,n2);
            pinit[1]= gsl_matrix_get(probatree,1,n2);
            pinit[2]= gsl_matrix_get(probatree,2,n2);
            pinit[3]= gsl_matrix_get(probatree,3,n2);
            gslinit=gsl_ran_discrete_preproc (4,pinit);
            base=gsl_ran_discrete (gslran,gslinit);
            if (n2<nbspecies) result[n2]=base;

            gsl_vector_set_zero(&pnoe2.vector);
            gsl_vector_set(&pnoe2.vector,base,1);

            //         for (unsigned int row=0;row<4;row++){
            //            for (unsigned int col=0;col<noemax+1;col++){
            //               cout << gsl_matrix_get(probatree,row,col) << " ";
            //            }
            //            cout << endl;
            //         }
            //         cout << endl;
         }

      }

      results[pos]=result;
      gsl_matrix_free(probatree);
      gsl_ran_discrete_free (gslinit);
   }

   vint idum2(mot.motwidth,0);
   vvint alignseq(nbspecies,idum2);
   for (unsigned int i=0;i<mot.motwidth;i++){
      for (unsigned int j=0;j<nbspecies;j++){
         alignseq[j][i]=results[i][j];
      }
   }

   return alignseq;
}

   void
test_evol_models_on_tree(Motif & refmot,unsigned int evol1,unsigned int evol2,ofstream & outf)
{

   Motif mot;
   mot.matprec=refmot.matprec;
   mot.seqs.clear();
   int iter;
   int logiter=0;

   for (unsigned int i=0;i<300;i++){

      vint initsite=refmot.refinstances.iseqs[i];

      Motalign ma;

      evolutionary_model=evol1;
      inittreedist();

      ma.alignseq=evolve_forward(refmot,initsite);

      vint matches(nbspecies,1);
      ma.matches=matches;
      //ma.print();
      mot.seqs.push_back(ma);

      int numinst=mot.seqs.size(); 
      iter=(int)pow(10,0.2*logiter);
      if (numinst<iter) continue;
      logiter++;

      //getchar();
      evolutionary_model=evol2;
      evolutionary_model=2;

      start=clock();
      mot.compprec_MCMC();
      //mot.compprec();
      finish=clock();
      dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
      printf ("---> time(mean): %f milliseconds.\n", dif );
      //printf ("---> time(opti): %f milliseconds.\n", dif );
      start=clock();
      //mot.compprec();
      finish=clock();
      dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
      //printf ("---> time(opti): %f milliseconds.\n", dif );
      mot.matfreq=mattofreq(mot.matprec);

      double distkl=0.;
      for (int pos=0;pos<mot.motwidth;pos++){
         for (int ib=0;ib<4;ib++){
            distkl+=mot.matfreq[pos][ib]*log(mot.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            //            cout << mot.matfreq[pos][ib] << " " <<
            //               refmot.matfreq[pos][ib] << " " <<
            //               distkl << endl;

         }
      }
      cout << evol1 << "_" << evol2 << " " << i+1 << " " << distkl << endl;
      //      displaymat(mot.matfreq);
      outf << evol1 << "_" << evol2 << " " << i+1 << " " << distkl << endl;
   }

}

   void
test_evol_models_on_tree_ALL(Motif & refmot,ofstream & outf)
{

   Motif mot1,mot2;
   mot1.matprec=refmot.matprec;
   mot2.matprec=refmot.matprec;
   int iter;
   int logiter=0;

   for (unsigned int i=0;i<refmot.refinstances.iseqs.size();i++){

      vint initsite=refmot.refinstances.iseqs[i];
      //cout << vinttostring(initsite) << endl;

      Motalign ma1,ma2;

      evolutionary_model=1;
      inittreedist();
      ma1.alignseq=evolve_forward(refmot,initsite);
      evolutionary_model=2;
      inittreedist();
      ma2.alignseq=evolve_forward(refmot,initsite);

      vint matches(nbspecies,1);
      ma1.matches=matches;
      ma2.matches=matches;
      mot1.seqs.push_back(ma1);
      mot2.seqs.push_back(ma2);

      int numinst=mot1.seqs.size(); 
      iter=(int)pow(10,0.2*logiter);
      if (numinst<iter) continue;
      logiter++;

      Motif mot11(mot1),mot12(mot1),mot21(mot2),mot22(mot2);
      start=clock();

      evolutionary_model=1;
      inittreedist();

      mot11.compprec_MCMC();
      mot11.matfreq=mattofreq(mot11.matprec);
      double distkl=0.;
      for (int pos=0;pos<refmot.motwidth;pos++){
         for (int ib=0;ib<4;ib++){
            distkl+=mot11.matfreq[pos][ib]*log(mot11.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            //cout << "1_1 " << mot11.matfreq[pos][ib] << " " << refmot.matfreq[pos][ib] << " " << distkl << endl;
         }
      }
      outf << "1_1" << " " << i+1 << " " << distkl << endl;

      mot21.compprec_MCMC();
      mot21.matfreq=mattofreq(mot21.matprec);
      distkl=0.;
      for (int pos=0;pos<refmot.motwidth;pos++){
         for (int ib=0;ib<4;ib++){
            distkl+=mot21.matfreq[pos][ib]*log(mot21.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            //cout << "2_1 " << mot21.matfreq[pos][ib] << " " << refmot.matfreq[pos][ib] << " " << distkl << endl;
         }
      }
      outf << "2_1" << " " << i+1 << " " << distkl << endl;

      evolutionary_model=2;
      inittreedist();

      mot12.compprec_MCMC();
      mot12.matfreq=mattofreq(mot12.matprec);
      distkl=0.;
      for (int pos=0;pos<refmot.motwidth;pos++){
         for (int ib=0;ib<4;ib++){
            distkl+=mot12.matfreq[pos][ib]*log(mot12.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            //cout << "1_2 " << mot12.matfreq[pos][ib] << " " << refmot.matfreq[pos][ib] << " " << distkl << endl;
         }
      }
      outf << "1_2" << " " << i+1 << " " << distkl << endl;

      mot22.compprec_MCMC();
      mot22.matfreq=mattofreq(mot22.matprec);
      distkl=0.;
      for (int pos=0;pos<refmot.motwidth;pos++){
         for (int ib=0;ib<4;ib++){
            distkl+=mot22.matfreq[pos][ib]*log(mot22.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
            //cout << "2_2 " << mot22.matfreq[pos][ib] << " " << refmot.matfreq[pos][ib] << " " << distkl << endl;
         }
      }
      outf << "2_2" << " " << i+1 << " " << distkl << endl;

      finish=clock();
      dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
      //   printf ("---> time(%d): %f milliseconds.\n", i+1, dif );
   }

}

   void
test_felsen_on_all_species(Motif & refmot,ofstream & outf)
{

   evolutionary_model=1;
   inittreedist();
   Motif mot1;
   Motalign ma1;
   vint matches(nbspecies,1);
   ma1.matches=matches;
   mot1.matprec=refmot.matprec;

   for (unsigned int i=0;i<refmot.refinstances.iseqs.size();i++){
      vint initsite=refmot.refinstances.iseqs[i];
      //cout << vinttostring(initsite) << endl;
      ma1.alignseq=evolve_forward(refmot,initsite);
      mot1.seqs.push_back(ma1);
   }

   start=clock();

   int iter(0);
   int logiter(0);
   for (unsigned int n=0;n<nbspecies;n++)
   {
      for (unsigned int k=0;k<5;k++){
         // RANDOM SHUFFLING
         random_shuffle(mot1.seqs.begin(),mot1.seqs.end());
         mot1.refinstances.seqs.clear();
         mot1.refinstances.iseqs.clear();
         logiter=0;
         for (int i=0;i<mot1.seqs.size();i++){
            if (mot1.seqs[i].matches[n]){
               mot1.refinstances.iseqs.push_back(mot1.seqs[i].alignseq[n]);
               mot1.refinstances.seqs.push_back(vinttostring(mot1.seqs[i].alignseq[n]));

               int numinst=mot1.refinstances.seqs.size(); 
               iter=(int)pow(10,0.1*logiter);
               if (numinst<iter) continue;
               logiter++;

               //               cout << "# of sites for ConsRef: " <<   mot1.refinstances.iseqs.size() << 
               //                  " (species=" << n+1 << "/" << nbspecies << ", iter=" << k+1 << "/" << 5 << ")" << endl;
               //               cout.flush();
               mot1.comprefmot();
               double distkl=0.;
               for (int pos=0;pos<mot1.motwidth;pos++){
                  for (int ib=0;ib<4;ib++){
                     distkl+=mot1.matfreq[pos][ib]*log(mot1.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
                  }
               }
               outf << "ConsRef_" << numtospecies(n) << " " << numinst << " " << distkl << endl;
               outf.flush();
            }
         }
      }
   }
   finish=clock();
   dif = 1000*(double)(finish - start) / CLOCKS_PER_SEC;
   //printf ("---> time(%d): %f milliseconds.\n", i+1, dif );

}

// test different models of evolution with synthetic TFBS
// (reconstruction of a PWM going down and up of a tree)
   void
phylotest(Motif & mot)
{

   //evolvebarrier(mot);   
   //evolvebystep(mot);   
   //evolvemean(mot);
   //evolveic(mot);
   Motif refmot=mot;
   ofstream outf("testcv.dat");

   cout << "Computing reference (common ancestor) PWM..." << endl;
   for (unsigned int i=0;i<1300;i++){
      // INITIAL SITE
      vint initsite(refmot.motwidth,4);
      while (scoref(initsite,refmot.matprec)<refmot.motscorethr2){
         for (unsigned int pos=0;pos<mot.motwidth;pos++){
            double pinit[4]= { mot.matfreq[pos][0],
               mot.matfreq[pos][1],
               mot.matfreq[pos][2],
               mot.matfreq[pos][3] };
            gsl_ran_discrete_t * gslinit=gsl_ran_discrete_preproc (4,pinit);
            unsigned int base=gsl_ran_discrete (gslran,gslinit);
            gsl_ran_discrete_free (gslinit);
            initsite[pos]=base;
         }
      }
      refmot.refinstances.iseqs.push_back(initsite);
      refmot.refinstances.seqs.push_back(vinttostring(initsite));
      //cout << vinttostring(initsite) << endl;
   }

   // UNCOMMENT TO USE REFMOT AS CONVERGENCE GOAL
   refmot.comprefmot();

   //   displaymat(refmot.matfreq);
   //   displaymat(mot.matfreq);

   Sequence initRefInstances=refmot.refinstances;
   Motif tmpmot;
   tmpmot=refmot;
   // /*
   cout << "Testing convergence on common ancestor..." << endl;
   for (unsigned int k=0;k<10;k++)
   {
      // RANDOM SHUFFLING
      random_shuffle(initRefInstances.iseqs.begin(),initRefInstances.iseqs.end());
      tmpmot.refinstances.iseqs.clear();
      tmpmot.refinstances.seqs.clear();
      int logiter=0;
      int iter;
      //for (ivvint iv=initRefInstances.iseqs.begin();iv!=initRefInstances.iseqs.begin()+300;iv++){}
      for (ivvint iv=initRefInstances.iseqs.begin();iv!=initRefInstances.iseqs.end();iv++){
         tmpmot.refinstances.iseqs.push_back(*iv);
         tmpmot.refinstances.seqs.push_back(vinttostring(*iv));

         // 0.1 step in log scale
         int numinst=tmpmot.refinstances.seqs.size(); 
         iter=(int)pow(10,0.1*logiter);
         if (numinst<iter) continue;
         logiter++;

         tmpmot.comprefmot();
         double distkl=0.;
         for (int pos=0;pos<tmpmot.motwidth;pos++){
            for (int ib=0;ib<4;ib++){
               distkl+=tmpmot.matfreq[pos][ib]*log(tmpmot.matfreq[pos][ib]/refmot.matfreq[pos][ib]);
               //cout << tmpmot.matfreq[pos][ib] << " " << refmot.matfreq[pos][ib] << " " << distkl << endl;
            }
         }
         outf << "AllRef " << numinst << " " << distkl << endl;
      }
   }

   cout << "Testing convergence on aligned species..." << endl;
   for (unsigned int k=0;k<10;k++){
      random_shuffle(refmot.refinstances.iseqs.begin(),refmot.refinstances.iseqs.end());
      test_felsen_on_all_species(refmot,outf);
   }

   cout << "Testing evolution models..." << endl;
   for (unsigned int k=0;k<10;k++)
   {
      random_shuffle(refmot.refinstances.iseqs.begin(),refmot.refinstances.iseqs.end());
      test_evol_models_on_tree_ALL(refmot,outf);
   }
   //   test_evol_models_on_tree(refmot,1,1,outf);
   //   test_evol_models_on_tree(refmot,1,2,outf);
   //   test_evol_models_on_tree(refmot,2,1,outf);
   //   test_evol_models_on_tree(refmot,2,2,outf);
   outf.close();

}

   void
evolvesite(Motif & mot)
{

//   vd dtest;
//   dtest.push_back(1.);
//   dtest.push_back(0.);
//   dtest.push_back(0.);
//   dtest.push_back(0.);
//
//   for (double dist=0.1;dist<10;dist+=0.1){
//      vd vtemp=evolvedist_felsen(dtest,mot.matfreq[1],dist);
//      cout << evolvedist_felsen_backwards(vtemp,mot.matfreq[1],dist) << endl;
//   }
//   exit(9);

   //number of sites for IC stat:
   unsigned int numforstat(50);

   ofstream outf("evol_site.dat");
   double time_step=0.01; // for site evolution with barrier

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
      distances[0]=0;
      distances[1]=0.1587;
      distances[2]=0.4653;
      distances[3]=0.4662;
      distances[4]=0.4708;
      distances[5]=0.496;
      distances[6]=0.5512;
      distances[7]=0.5512;
      distances[8]=0.5397;
      distances[9]=0.4971;
   }


   evolutionary_model=2;
   inittreedist();
   vd dum(4,0.);
   vvd vdum(4,dum);
   vvd wdum(mot.motwidth,dum);
   vvvd M_h(mot.motwidth,vdum);
   vvvd M_f(mot.motwidth,vdum);

   for (unsigned int i=0;i<mot.motwidth;i++){
      M_h[i]=instant_rates_halpern (mot.matfreq[i], time_step);
      M_f[i]=instant_rates_felsen (mot.matfreq[i], time_step);
      //cout << M[i] << endl;
   }

   double initscore;
//   for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){
//      ivm->print();
//   }

   for (unsigned int n=1;n<nbspecies;n++){
      cout << numtospecies(n) << endl;
      for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){

         if (!ivm->matches[n]) continue;

         initscore=scoref(ivm->alignseq[0],mot.matprec);

         vint site(ivm->alignseq[0]);

         // previous m initialized to Id
         vvd pm(vdum);
         pm[0][0]=1.;
         pm[1][1]=1.;
         pm[2][2]=1.;
         pm[3][3]=1.;

         // initial probability initialized to TFBS
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

         // WITH BARRIER
         ///*
         vvint sites_w_barrier_h(numforstat,site);
         vvint sites_w_barrier_f(numforstat,site);
         //
         for (double dist=time_step;dist<=distances[n];dist+=time_step){

            // HALPERN:
            for (ivvint ivv=sites_w_barrier_h.begin();ivv!=sites_w_barrier_h.end();ivv++){
               site=*ivv;
               vint tempsite(mot.motwidth,4);
               while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
                  for (unsigned int i=0;i<mot.motwidth;i++){
                     vd psite_init(4,0.);
                     psite_init[site[i]]=1.;
                     vd psite_evol(4,0.);
                     for (unsigned int row=0;row<4;row++){
                        for (unsigned int k=0;k<4;k++){
                           psite_evol[row]+=M_h[i][k][row]*psite_init[k];
                        }
                     }
                     double pdist[4]={psite_evol[0],psite_evol[1],psite_evol[2],psite_evol[3]};
                     gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
                     unsigned int base=gsl_ran_discrete (gslran,gsldist);
                     gsl_ran_discrete_free (gsldist);
                     (*ivv)[i]=base;
                  }
                  tempsite=*ivv;
               }
               *ivv=tempsite;
            }
            // FELSEN:
            for (ivvint ivv=sites_w_barrier_f.begin();ivv!=sites_w_barrier_f.end();ivv++){
               site=*ivv;
               vint tempsite(mot.motwidth,4);
               while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
                  for (unsigned int i=0;i<mot.motwidth;i++){
                     vd psite_init(4,0.);
                     psite_init[site[i]]=1.;
                     vd psite_evol(4,0.);
                     for (unsigned int row=0;row<4;row++){
                        for (unsigned int k=0;k<4;k++){
                           psite_evol[row]+=M_f[i][k][row]*psite_init[k];
                        }
                     }
                     double pdist[4]={psite_evol[0],psite_evol[1],psite_evol[2],psite_evol[3]};
                     gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
                     unsigned int base=gsl_ran_discrete (gslran,gsldist);
                     gsl_ran_discrete_free (gsldist);
                     (*ivv)[i]=base;
                  }
                  tempsite=*ivv;
               }
               *ivv=tempsite;
            }
         }
          //*/


         // SCORES COMPUTATION
         //
         // First, the "true" TFBS
         //

         outf << numtospecies(n) << " ";
         outf << "Data" << " ";
         outf << scoref(ivm->alignseq[n],mot.matprec)-initscore << endl;


         // draw numforstat sites for ICs computation
         //cout << "data " << vinttostring(ivm->alignseq[n]) << endl;
         for (unsigned int ks=0;ks<numforstat;ks++){

            vint dumsite(mot.motwidth,4);

            // HALPERN
            vint tempsite(dumsite);
            while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
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
            outf << scoref(site,mot.matprec)-initscore << endl;

//            if (ks==1)
//            cout << vinttostring(ivm->alignseq[n]) << " " << vinttostring(site) << " (Halpern)" << endl;

            //cout << "halpern " << vinttostring(site) << endl;

            // FELSEN
            tempsite=dumsite;
            while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
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
            outf << scoref(site,mot.matprec)-initscore << endl;
            
//            if (ks==1)
//            cout << vinttostring(ivm->alignseq[n]) << " " << vinttostring(site) << " (Felsen)" << endl;

            // HALPERN WITH BARRIERS
         //   /*
               outf << numtospecies(n) << " ";
               outf << "Halpern_w_barrier" << " ";
               outf << scoref(sites_w_barrier_h[ks],mot.matprec)-initscore << endl;

            // FELSEN WITH BARRIERS
            outf << numtospecies(n) << " ";
            outf << "Felsen_w_barrier" << " ";
            outf << scoref(sites_w_barrier_f[ks],mot.matprec)-initscore << endl;
            //cout << "felsen " << vinttostring(site) << endl;
//            */


         }
         //         getchar();
      }
   }

   outf.close();
}

   void
evolvesite(vmot & mots)
{

   unsigned int isref=0;
   if (mots.size()>1) isref=1;

   Motif refmot=mots[0];


   //number of sites for IC stat:
   unsigned int numforstat(100);

   ofstream outf("evol_site.dat");
   double time_step=0.01; // for site evolution with barrier

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
      distances[0]=0;
      distances[1]=0.1587;
      distances[2]=0.4653;
      distances[3]=0.4662;
      distances[4]=0.4708;
      distances[5]=0.496;
      distances[6]=0.5512;
      distances[7]=0.5512;
      distances[8]=0.5397;
      distances[9]=0.4971;
   }

   for (ivmot ivmo=mots.begin()+isref;ivmo!=mots.end();ivmo++){

      Motif mot=(*ivmo);

      cout << mot.name << endl;

      evolutionary_model=2;
      inittreedist();
      vd dum(4,0.);
      vvd vdum(4,dum);
      vvd wdum(mot.motwidth,dum);
      vvvd M_h(mot.motwidth,vdum);
      vvvd M_f(mot.motwidth,vdum);

      for (unsigned int i=0;i<mot.motwidth;i++){
         M_h[i]=instant_rates_halpern (mot.matfreq[i], time_step);
         M_f[i]=instant_rates_felsen (mot.matfreq[i], time_step);
         //cout << M[i] << endl;
      }

      double initscore;

      for (unsigned int n=1;n<nbspecies;n++){
         cout << numtospecies(n) << endl;

         //for (ivma ivm=refmot.seqs.begin();ivm!=min(refmot.seqs.begin()+numforstat,refmot.seqs.end());ivm++){}
         for (ivma ivm=refmot.seqs.begin();ivm!=refmot.seqs.end();ivm++){

            unsigned int goon=1;
            //            for (ivma ivm1=refmot.seqs.begin();ivm1!=refmot.seqs.end();ivm1++){
            //               if (ivm1->seq_start==ivm->seq_start && ivm1->seq_stop==ivm->seq_stop){
            //                  goon=1;
            //               }
            //            }
            //            if (goon==0) cout << "NO OVERLAP" << endl;
            //
            for (ivmot ivm1=mots.begin()+isref;ivm1!=mots.end();ivm1++){
               if (mot.name!=ivm1->name && scoref(ivm->alignseq[0],ivm1->matprec)>scoref(ivm->alignseq[0],mot.matprec)){
                  goon=0;
               }
            }

            if (!goon) continue;



            if (!ivm->matches[n]) continue;

            initscore=scoref(ivm->alignseq[0],refmot.matprec);

            vint site(ivm->alignseq[0]);

            // previous m initialized to Id
            vvd pm(vdum);
            pm[0][0]=1.;
            pm[1][1]=1.;
            pm[2][2]=1.;
            pm[3][3]=1.;

            // initial probability initialized to TFBS
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

            // WITH BARRIER
            /*
               vvint sites_w_barrier_h(numforstat,site);
               vvint sites_w_barrier_f(numforstat,site);
            //
            for (double dist=time_step;dist<=distances[n];dist+=time_step){

            // HALPERN:
            for (ivvint ivv=sites_w_barrier_h.begin();ivv!=sites_w_barrier_h.end();ivv++){
            site=*ivv;
            vint tempsite(mot.motwidth,4);
            while (scoref(tempsite,refmot.matprec)<refmot.motscorethr2){
            for (unsigned int i=0;i<mot.motwidth;i++){
            vd psite_init(4,0.);
            psite_init[site[i]]=1.;
            vd psite_evol(4,0.);
            for (unsigned int row=0;row<4;row++){
            for (unsigned int k=0;k<4;k++){
            psite_evol[row]+=M_h[i][k][row]*psite_init[k];
            }
            }
            double pdist[4]={psite_evol[0],psite_evol[1],psite_evol[2],psite_evol[3]};
            gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
            unsigned int base=gsl_ran_discrete (gslran,gsldist);
            gsl_ran_discrete_free (gsldist);
            (*ivv)[i]=base;
            }
            tempsite=*ivv;
            }
             *ivv=tempsite;
             }
            // FELSEN:
            for (ivvint ivv=sites_w_barrier_f.begin();ivv!=sites_w_barrier_f.end();ivv++){
            site=*ivv;
            vint tempsite(mot.motwidth,4);
            while (scoref(tempsite,refmot.matprec)<refmot.motscorethr2){
            for (unsigned int i=0;i<mot.motwidth;i++){
            vd psite_init(4,0.);
            psite_init[site[i]]=1.;
            vd psite_evol(4,0.);
            for (unsigned int row=0;row<4;row++){
            for (unsigned int k=0;k<4;k++){
            psite_evol[row]+=M_f[i][k][row]*psite_init[k];
            }
            }
            double pdist[4]={psite_evol[0],psite_evol[1],psite_evol[2],psite_evol[3]};
            gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
            unsigned int base=gsl_ran_discrete (gslran,gsldist);
            gsl_ran_discrete_free (gsldist);
            (*ivv)[i]=base;
            }
            tempsite=*ivv;
            }
             *ivv=tempsite;
             }
             }
             */


            // SCORES COMPUTATION
            //
            // First, the "true" TFBS
            //

            outf << numtospecies(n) << " ";
            outf << "Data" << " ";
            outf << scoref(ivm->alignseq[n],refmot.matprec)-initscore << endl;


            // draw numforstat sites for ICs computation
            //cout << "data " << vinttostring(ivm->alignseq[n]) << endl;
            for (unsigned int ks=0;ks<numforstat;ks++){

               vint dumsite(mot.motwidth,4);

               // HALPERN
               vint tempsite(dumsite);
               while (scoref(tempsite,refmot.matprec)<refmot.motscorethr2){
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
               outf << scoref(site,refmot.matprec)-initscore << endl;

               //cout << "halpern " << vinttostring(site) << endl;
               // FELSEN
               tempsite=dumsite;
               while (scoref(tempsite,refmot.matprec)<refmot.motscorethr2){
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
               outf << scoref(site,refmot.matprec)-initscore << endl;

               /*
               // HALPERN WITH BARRIERS
               outf << numtospecies(n) << " ";
               outf << "Halpern_w_barrier" << " ";
               outf << scoref(sites_w_barrier_h[k],refmot.matprec)-initscore << endl;


               // FELSEN WITH BARRIERS
               outf << numtospecies(n) << " ";
               outf << "Felsen_w_barrier" << " ";
               outf << scoref(sites_w_barrier_f[k],refmot.matprec)-initscore << endl;
               //cout << "felsen " << vinttostring(site) << endl;
               */

            }
            //         getchar();
         }
      }

      cout << endl;
   }

   outf.close();
}

// evolve base with different phylogenetic models
   void
evolvebase(Motif & mot)
{
   ofstream outf("evol_base.dat");
   double time_step=0.01;

   //unsigned int pos=3;
   //unsigned int base=0;
   for (unsigned int pos=0;pos<mot.motwidth;pos++){
      for (unsigned int base=0;base<4;base++){
         vd w=mot.matfreq[pos];
         vd pdisp(4,0);
         pdisp[base]=1;

         vd initprob(4,0);
         initprob[base]=1.;
         gsl_vector * pnoe=gsl_vector_alloc(4);
         gsl_vector_set_zero(pnoe);
         gsl_vector_set(pnoe,base,1.0);

         evolutionary_model=2;
         inittreedist();
         //         gsl_matrix * M=gsl_matrix_alloc(4,4);
         //         gsl_matrix_memcpy(pij,id);
         //         gsl_matrix * pmattemp;

         //         M=instant_rates_halpern (w, time_step);

         vd dum(4,0.);
         vvd m(4,dum);
         vvd pm(4,dum);
         m=instant_rates_halpern (w, time_step);
         //         for (unsigned int row=0;row<4;row++){
         //            for (unsigned int col=0;col<4;col++){
         //               m[col][row]=gsl_matrix_get(M,row,col);
         //            }
         //         }
         pm[0][0]=1.;
         pm[1][1]=1.;
         pm[2][2]=1.;
         pm[3][3]=1.;
         double distkl;

         for (double dist=time_step;dist<2.;dist+=time_step){

            // compute M^n      
            //      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,M,pij,0.0,pijp);
            //      gsl_blas_dgemv(CblasNoTrans,1.0,pijp,pnoe,0.0,pnoe);
            //      gsl_matrix_memcpy(pij,pijp);
            //      pdisp[0]= gsl_vector_get(pnoe,0);
            //      pdisp[1]= gsl_vector_get(pnoe,1);
            //      pdisp[2]= gsl_vector_get(pnoe,2);
            //      pdisp[3]= gsl_vector_get(pnoe,3);
            //            vvd pm2(4,dum);
            //            for (unsigned int row=0;row<4;row++){
            //               for (unsigned int col=0;col<4;col++){
            //                  for (unsigned int k=0;k<4;k++){
            //                     pm2[col][row]+=m[k][row]*pm[col][k];
            //                  }
            //               }
            //            }
            //            pm=pm2;
            //            //      cout << pm << endl;
            //            vd pdisp2(4,0.);
            //            for (unsigned int row=0;row<4;row++){
            //               for (unsigned int k=0;k<4;k++){
            //                  pdisp2[row]+=pm[k][row]*pdisp[k];
            //               }
            //            }
            //            pdisp=pdisp2;

            pdisp=evolvedist_halpern(initprob,w,dist);

            distkl=0;
            outf << "Halpern" << " ";
            outf << pos << " ";
            outf << base << " ";
            outf << dist << " ";
            for (int ib=0;ib<4;ib++){
               distkl+=pdisp[ib]*log(pdisp[ib]/w[ib]);
            }
            outf << distkl << endl;

            //      for (unsigned int row=0;row<4;row++){
            //         for (unsigned int col=0;col<4;col++){
            //            //cout << gsl_matrix_get(pijp,row,col) << " ";
            //            cout << pm[col][row] << " ";
            //         }
            //         cout << endl;
            //      }
            //      cout << endl;
            //      for (unsigned int col=0;col<4;col++){
            //         //cout << gsl_vector_get(pnoe,col) << " ";
            //         cout << pdisp[col] << " ";
            //      }
            //      cout << endl;
            //      for (unsigned int col=0;col<4;col++){
            //         cout << w[col] << " ";
            //      }
            //      cout << endl;
            //      cout << endl;
            //getchar();


            pdisp=evolvedist_felsen(initprob,w,dist);
            outf << "Felsen" << " ";
            outf << pos << " ";
            outf << base << " ";
            outf << dist << " ";
            distkl=0;
            for (int ib=0;ib<4;ib++){
               distkl+=pdisp[ib]*log(pdisp[ib]/w[ib]);
            }

            outf << distkl << endl;
         }

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
            distances[0]=0;
            distances[1]=0.1587;
            distances[2]=0.4653;
            distances[3]=0.4662;
            distances[4]=0.4708;
            distances[5]=0.496;
            distances[6]=0.5512;
            distances[7]=0.5512;
            distances[8]=0.5397;
            distances[9]=0.4971;
         }

         for (unsigned int n=1;n<nbspecies;n++){
            vint counts(4,0);
            for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){
               if (ivm->matches[n] && ivm->alignseq[0][pos]==base){
                  counts[ivm->alignseq[n][pos]]++;
               }
            }
            unsigned int totcounts=0;
            totcounts=accumulate(counts.begin(),counts.end(),0); 
            if (totcounts<10) continue;
            cout << numtospecies(n) << " " << totcounts << endl;

            pdisp[0]=(counts[0]+alpha)/(totcounts+2*alpha+2*beta);
            pdisp[1]=(counts[1]+alpha)/(totcounts+2*alpha+2*beta);
            pdisp[2]=(counts[2]+beta)/(totcounts+2*alpha+2*beta);
            pdisp[3]=(counts[3]+beta)/(totcounts+2*alpha+2*beta);

            outf << numtospecies(n) << " ";
            outf << pos << " ";
            outf << base << " ";
            outf << distances[n] << " ";
            distkl=0;
            for (int ib=0;ib<4;ib++){
               distkl+=pdisp[ib]*log(pdisp[ib]/w[ib]);
            }
            outf << distkl << endl;

         }
      }
   }
   outf.close();
}


// fit evol dist per base on phylogenetic models
   void
fitdistperbase(Motif & mot)
{
   ofstream outf("fitdist_base.dat");
   double time_step=0.01;
   double distkl,distkl_f;
   double taubest,btotbest(1e6);
   vd dum(4,0.);
   vd w(dum),pdisp(dum);

   //distances to ref
   double distance;
   vd distances(nbspecies,0);
   vd newdist(nbspecies,0);
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
      distances[0]=0;
      distances[1]=0.1587;
      distances[2]=0.4653;
      distances[3]=0.4662;
      distances[4]=0.4708;
      distances[5]=0.496;
      distances[6]=0.5512;
      distances[7]=0.5512;
      distances[8]=0.5397;
      distances[9]=0.4971;
   }

   evolutionary_model=2;
   inittreedist();


   ofstream chi2("chi2.dat");
   for (double tau=0.1;tau<10;tau+=0.05){
      double btot(0);
      for (unsigned int pos=0;pos<mot.motwidth;pos++){
         for (unsigned int base=0;base<4;base++){

            // modified distances

            for (unsigned int id=0;id<nbspecies;id++){
               newdist[id]=distances[id]/tau;
            }

            // compute the quantity b2= sum ( distkl_obs - distkl_felsen ) ^ 2

            double b2(0);
            w=mot.matfreq[pos];
            vd initprob(4,0);
            initprob[base]=1.;

            for (unsigned int n=1;n<nbspecies;n++){

               // 1. distkl_obs

               vint counts(4,0);
               for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){
                  if (ivm->matches[n] && ivm->alignseq[0][pos]==base){
                     counts[ivm->alignseq[n][pos]]++;
                  }
               }
               unsigned int totcounts=0;
               totcounts=accumulate(counts.begin(),counts.end(),0); 
               //               if (totcounts<10) continue;

               pdisp[0]=(counts[0]+alpha)/(totcounts+2*alpha+2*beta);
               pdisp[1]=(counts[1]+alpha)/(totcounts+2*alpha+2*beta);
               pdisp[2]=(counts[2]+beta)/(totcounts+2*alpha+2*beta);
               pdisp[3]=(counts[3]+beta)/(totcounts+2*alpha+2*beta);

               distkl=0;
               for (int ib=0;ib<4;ib++){
                  distkl+=pdisp[ib]*log(pdisp[ib]/w[ib]);
               }

               // 2. distkl_f
               pdisp=evolvedist_felsen(initprob,w,newdist[n]);
               //pdisp=evolvedist_halpern(initprob,w,newdist[n]);
               distkl_f=0;
               for (int ib=0;ib<4;ib++){
                  distkl_f+=pdisp[ib]*log(pdisp[ib]/w[ib]);
               }

               b2+=pow(distkl-distkl_f,2);
            }

            btot+=b2;
         }
      }
      cout << tau << "\t" << btot << endl;
      chi2 << tau << "\t" << btot << endl;

      if (btot<btotbest){
         taubest=tau;
         btotbest=btot;
      }
   }
   chi2.close();

   for (unsigned int pos=0;pos<mot.motwidth;pos++){
      for (unsigned int base=0;base<4;base++){
         w=mot.matfreq[pos];

         vd initprob(4,0);
         initprob[base]=1.;

         for (unsigned int id=0;id<nbspecies;id++){
            newdist[id]=distances[id]/taubest;
         }


         for (double dist=time_step;dist<2.;dist+=time_step){

            distance=dist/taubest;
            pdisp=evolvedist_halpern(initprob,w,distance);
            outf << "Halpern" << " ";
            outf << pos << " ";
            outf << base << " ";
            outf << distance << " ";
            distkl=0;
            for (int ib=0;ib<4;ib++){
               distkl+=pdisp[ib]*log(pdisp[ib]/w[ib]);
            }
            outf << distkl << endl;

            pdisp=evolvedist_felsen(initprob,w,distance);
            outf << "Felsen" << " ";
            outf << pos << " ";
            outf << base << " ";
            outf << distance << " ";
            distkl=0;
            for (int ib=0;ib<4;ib++){
               distkl+=pdisp[ib]*log(pdisp[ib]/w[ib]);
            }

            outf << distkl << endl;
         }


         for (unsigned int n=1;n<nbspecies;n++){
            vint counts(4,0);
            for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){
               if (ivm->matches[n] && ivm->alignseq[0][pos]==base){
                  counts[ivm->alignseq[n][pos]]++;
               }
            }
            unsigned int totcounts=0;
            totcounts=accumulate(counts.begin(),counts.end(),0); 
            if (totcounts<10) continue;
            cout << numtospecies(n) << " " << totcounts << endl;

            pdisp[0]=(counts[0]+alpha)/(totcounts+2*alpha+2*beta);
            pdisp[1]=(counts[1]+alpha)/(totcounts+2*alpha+2*beta);
            pdisp[2]=(counts[2]+beta)/(totcounts+2*alpha+2*beta);
            pdisp[3]=(counts[3]+beta)/(totcounts+2*alpha+2*beta);

            outf << numtospecies(n) << " ";
            outf << pos << " ";
            outf << base << " ";
            outf << newdist[n] << " ";
            distkl=0;
            for (int ib=0;ib<4;ib++){
               distkl+=pdisp[ib]*log(pdisp[ib]/w[ib]);
            }
            outf << distkl << endl;

         }
      }
   }
   outf.close();
}

   void
fitdistpersite(Motif & mot)
{

   //number of sites for IC stat:
   unsigned int numforstat(50);

   double time_step=0.01; // for barrier
   vd dum(4,0.);
   vd w(dum),pdisp(dum);

   //distances to ref
   double distance;
   vd distances(nbspecies,0);
   vd newdist(nbspecies,0);
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
      distances[0]=0;
      distances[1]=0.1587;
      distances[2]=0.4653;
      distances[3]=0.4662;
      distances[4]=0.4708;
      distances[5]=0.496;
      distances[6]=0.5512;
      distances[7]=0.5512;
      distances[8]=0.5397;
      distances[9]=0.4971;
   }


   evolutionary_model=2;
   inittreedist();
   vvd vdum(4,dum);
   vvd wdum(mot.motwidth,dum);
   vvvd M_h(mot.motwidth,vdum);
   vvvd M_f(mot.motwidth,vdum);

   for (unsigned int i=0;i<mot.motwidth;i++){
      M_h[i]=instant_rates_halpern (mot.matfreq[i], time_step);
      M_f[i]=instant_rates_felsen (mot.matfreq[i], time_step);
   }

   double initscore;

   // find best tau for Felsen model
   cout << "Finding optimal tau..." << endl;
   double taubest,btotbest(1e6);
   ofstream chi2("chi2.dat");
   for (double tau=0.3;tau<20;tau+=0.05){
      for (unsigned int id=0;id<nbspecies;id++){
         newdist[id]=distances[id]/tau;
      }
      int index=0;
      double btot(0);
      for (unsigned int n=1;n<nbspecies;n++){
         vd distscore(80,0.); // between -4 and +4 with step 0.1
         vd distscore_f(80,0.); // between -4 and +4 with step 0.1

         ivma ivmstop=min(mot.seqs.begin()+numforstat,mot.seqs.end());

         for (ivma ivm=mot.seqs.begin();ivm!=ivmstop;ivm++){

            if (!ivm->matches[n]) continue;

            initscore=scoref(ivm->alignseq[0],mot.matprec);

            // previous m initialized to Id
            vvd pm(vdum);
            pm[0][0]=1.;
            pm[1][1]=1.;
            pm[2][2]=1.;
            pm[3][3]=1.;

            // initial probability initialized to TFBS
            vvd Pdisp_f(wdum);
            vvd Initprob(wdum);
            for (unsigned int i=0;i<mot.motwidth;i++) 
               Initprob[i][ivm->alignseq[0][i]]=1;
            Pdisp_f=Initprob;
            vint site(ivm->alignseq[0]);

            // FELSEN
            for (unsigned int i=0;i<mot.motwidth;i++){
               Pdisp_f[i]=evolvedist_felsen(Initprob[i],mot.matfreq[i],newdist[n]);
            }

            // SCORES COMPUTATION
            // First, the "true" TFBS
            index=(int)(10*(4+(scoref(ivm->alignseq[n],mot.matprec)-initscore)));
            if (index>=0 && index<80) distscore[index]++;
            // draw numforstat sites for ICs computation
            for (unsigned int ks=0;ks<numforstat;ks++){

               vint dumsite(mot.motwidth,4);
               vint tempsite(dumsite);

               // FELSEN
               tempsite=dumsite;
               while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
                  for (unsigned int j=0;j<mot.motwidth;j++){
                     double pdist[4]={Pdisp_f[j][0],Pdisp_f[j][1],Pdisp_f[j][2],Pdisp_f[j][3]};
                     gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
                     unsigned int base=gsl_ran_discrete (gslran,gsldist);
                     gsl_ran_discrete_free (gsldist);
                     site[j]=base;
                  }
                  tempsite=site;
               }
               index=(int)(10*(4+(scoref(site,mot.matprec)-initscore)));
               if (index>=0 && index<80) distscore_f[index]+=1./numforstat;

            }
         }

         for (unsigned int id=0;id!=distscore.size();id++){
            btot+=pow(distscore[id]-distscore_f[id],2);
         }
      }
      chi2 << tau << "\t" << btot << endl << flush;
      cout << tau << "\t" << btot << endl << flush;
      if (btot<btotbest){
         taubest=tau;
         btotbest=btot;
      }
   }
   cout << "Best tau is " << taubest << endl;
   chi2.close();


   cout << "Computing plots at best tau..." << endl;
   ofstream outf("fitdist_site.dat");
   for (unsigned int id=0;id<nbspecies;id++){
      newdist[id]=distances[id]/taubest;
   }
   for (unsigned int n=1;n<nbspecies;n++){
      for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){
         // for (ivma ivm=mot.seqs.begin();ivm!=min(mot.seqs.begin()+numforstat,mot.seqs.end());ivm++){}

         if (!ivm->matches[n]) continue;

         initscore=scoref(ivm->alignseq[0],mot.matprec);

         // previous m initialized to Id
         vvd pm(vdum);
         pm[0][0]=1.;
         pm[1][1]=1.;
         pm[2][2]=1.;
         pm[3][3]=1.;

         // initial probability initialized to TFBS
         vvd Pdisp_f(wdum);
         vvd Pdisp_h(wdum);
         vvd Initprob(wdum);
         for (unsigned int i=0;i<mot.motwidth;i++) 
            Initprob[i][ivm->alignseq[0][i]]=1;
         Pdisp_f=Initprob;
         Pdisp_h=Initprob;


         // HALPERN
         for (unsigned int i=0;i<mot.motwidth;i++){
            Pdisp_h[i]=evolvedist_halpern(Initprob[i],mot.matfreq[i],newdist[n]);
         }

         // FELSEN
         for (unsigned int i=0;i<mot.motwidth;i++){
            Pdisp_f[i]=evolvedist_felsen(Initprob[i],mot.matfreq[i],newdist[n]);
         }

         // WITH BARRIER
         vint site(ivm->alignseq[0]);
         //         vvint sites_w_barrier_h(numforstat,site);
         //         vvint sites_w_barrier_f(numforstat,site);
         //
         //         for (double dist=time_step;dist<=newdist[n];dist+=time_step){
         //
         //            // HALPERN:
         //            for (ivvint ivv=sites_w_barrier_h.begin();ivv!=sites_w_barrier_h.end();ivv++){
         //               site=*ivv;
         //               vint tempsite(mot.motwidth,4);
         //               while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
         //                  for (unsigned int i=0;i<mot.motwidth;i++){
         //                     vd psite_init(4,0.);
         //                     psite_init[site[i]]=1.;
         //                     vd psite_evol(4,0.);
         //                     for (unsigned int row=0;row<4;row++){
         //                        for (unsigned int k=0;k<4;k++){
         //                           psite_evol[row]+=M_h[i][k][row]*psite_init[k];
         //                        }
         //                     }
         //                     double pdist[4]={psite_evol[0],psite_evol[1],psite_evol[2],psite_evol[3]};
         //                     gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
         //                     unsigned int base=gsl_ran_discrete (gslran,gsldist);
         //                     gsl_ran_discrete_free (gsldist);
         //                     (*ivv)[i]=base;
         //                  }
         //                  tempsite=*ivv;
         //               }
         //               *ivv=tempsite;
         //            }
         //            // FELSEN:
         //            for (ivvint ivv=sites_w_barrier_f.begin();ivv!=sites_w_barrier_f.end();ivv++){
         //               site=*ivv;
         //               vint tempsite(mot.motwidth,4);
         //               while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
         //                  for (unsigned int i=0;i<mot.motwidth;i++){
         //                     vd psite_init(4,0.);
         //                     psite_init[site[i]]=1.;
         //                     vd psite_evol(4,0.);
         //                     for (unsigned int row=0;row<4;row++){
         //                        for (unsigned int k=0;k<4;k++){
         //                           psite_evol[row]+=M_f[i][k][row]*psite_init[k];
         //                        }
         //                     }
         //                     double pdist[4]={psite_evol[0],psite_evol[1],psite_evol[2],psite_evol[3]};
         //                     gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
         //                     unsigned int base=gsl_ran_discrete (gslran,gsldist);
         //                     gsl_ran_discrete_free (gsldist);
         //                     (*ivv)[i]=base;
         //                  }
         //                  tempsite=*ivv;
         //               }
         //               *ivv=tempsite;
         //            }
         //         }


         // SCORES COMPUTATION
         //
         // First, the "true" TFBS
         //

         outf << numtospecies(n) << " ";
         outf << "Data" << " ";
         outf << scoref(ivm->alignseq[n],mot.matprec)-initscore << endl;


         // draw numforstat sites for ICs computation
         //cout << "data " << vinttostring(ivm->alignseq[n]) << endl;
         for (unsigned int ks=0;ks<numforstat;ks++){

            vint dumsite(mot.motwidth,4);

            // HALPERN
            vint tempsite(dumsite);
            while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
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
            outf << scoref(site,mot.matprec)-initscore << endl;

            //cout << "halpern " << vinttostring(site) << endl;

            // HALPERN WITH BARRIERS
            //            outf << numtospecies(n) << " ";
            //            outf << "Halpern_w_barrier" << " ";
            //            outf << scoref(sites_w_barrier_h[k],mot.matprec)-initscore << endl;

            // FELSEN
            tempsite=dumsite;
            while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
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
            outf << scoref(site,mot.matprec)-initscore << endl;

            // FELSEN WITH BARRIERS
            //            outf << numtospecies(n) << " ";
            //            outf << "Felsen_w_barrier" << " ";
            //            outf << scoref(sites_w_barrier_f[k],mot.matprec)-initscore << endl;
            //cout << "felsen " << vinttostring(site) << endl;


         }
         //         getchar();
      }
   }
   outf.close();
}

   void
fitkmeanspersite(Motif & mot,const char * filename)
{
      
   ifstream inf;
   inf.open(filename);
   vstring motsfile;
   back_insert_iterator<vstring> dest(motsfile);
   copy(iisstring(inf),iisstring(),dest);
   inf.close();

   unsigned int numk;
   numk=motsfile.size();

   //number of sites for IC stat:
   unsigned int numforstat(80);

   double time_step=0.01; // for barrier
   vd dum(4,0.);
   vd w(dum),pdisp(dum);

   //distances to ref
   double distance;
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
      distances[0]=0;
      distances[1]=0.1587;
      distances[2]=0.4653;
      distances[3]=0.4662;
      distances[4]=0.4708;
      distances[5]=0.496;
      distances[6]=0.5512;
      distances[7]=0.5512;
      distances[8]=0.5397;
      distances[9]=0.4971;
   }


   evolutionary_model=2;
   inittreedist();
   vvd vdum(4,dum);
   vvd wdum(mot.motwidth,dum);
   vvvd M_h(mot.motwidth,vdum);
   vvvd M_f(mot.motwidth,vdum);

   for (unsigned int i=0;i<mot.motwidth;i++){
      M_h[i]=instant_rates_halpern (mot.matfreq[i], time_step);
      M_f[i]=instant_rates_felsen (mot.matfreq[i], time_step);
   }

   double initscore;

   // find best K for Felsen model
   cout << "Finding optimal K..." << endl;
   double kbest,btotbest(1e6);
   ofstream chi2("chi2.dat");
   unsigned int fcount=0,fcountbest(0);
   for (unsigned int k=2;k<numk;k++){

      vmot motmix;
      loadmotswnames(motsfile[fcount].c_str(),motmix);
      for (ivmot ivm=motmix.begin();ivm!=motmix.end();ivm++){
         ivm->setscorethr2meaninfo();
      }

      int index=0;
      double btot(0);
      for (unsigned int n=1;n<nbspecies;n++){
         vd distscore(80,0.); // between -4 and +4 with step 0.1
         vd distscore_f(80,0.); // between -4 and +4 with step 0.1

         //ivma ivmstop=min(mot.seqs.begin()+numforstat,mot.seqs.end());
         ivma ivmstop=mot.seqs.end();

         for (ivma ivm=mot.seqs.begin();ivm!=ivmstop;ivm++){

            if (!ivm->matches[n]) continue;

            initscore=scoref(ivm->alignseq[0],mot.matprec);
            
            double scoretmp=-1e6;
            Motif motevo;
            for (ivmot ivm1=motmix.begin()+1;ivm1!=motmix.end();ivm1++){
               if (scoref(ivm->alignseq[0],ivm1->matprec)>scoretmp){
                  scoretmp=scoref(ivm->alignseq[0],ivm1->matprec);
                  motevo=*ivm1;
               }
            }

            // previous m initialized to Id
            vvd pm(vdum);
            pm[0][0]=1.;
            pm[1][1]=1.;
            pm[2][2]=1.;
            pm[3][3]=1.;

            // initial probability initialized to TFBS
            vvd Pdisp_f(wdum);
            vvd Initprob(wdum);
            for (unsigned int i=0;i<motevo.motwidth;i++) 
               Initprob[i][ivm->alignseq[0][i]]=1;
            Pdisp_f=Initprob;
            vint site(ivm->alignseq[0]);

            // FELSEN
            for (unsigned int i=0;i<mot.motwidth;i++){
               Pdisp_f[i]=evolvedist_felsen(Initprob[i],motevo.matfreq[i],distances[n]);
            }

            // SCORES COMPUTATION
            // First, the "true" TFBS
            index=(int)(10*(4+(scoref(ivm->alignseq[n],mot.matprec)-initscore)));
            if (index>=0 && index<80) distscore[index]++;
            // draw numforstat sites for ICs computation
            for (unsigned int ks=0;ks<numforstat;ks++){

               vint dumsite(motevo.motwidth,4);
               vint tempsite(dumsite);

               // FELSEN
               tempsite=dumsite;
               while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
                  for (unsigned int j=0;j<motevo.motwidth;j++){
                     double pdist[4]={Pdisp_f[j][0],Pdisp_f[j][1],Pdisp_f[j][2],Pdisp_f[j][3]};
                     gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
                     unsigned int base=gsl_ran_discrete (gslran,gsldist);
                     gsl_ran_discrete_free (gsldist);
                     site[j]=base;
                  }
                  tempsite=site;
               }
               index=(int)(10*(4+(scoref(site,mot.matprec)-initscore)));
               if (index>=0 && index<80) distscore_f[index]+=1./numforstat;

            }
         }

         for (unsigned int id=0;id!=distscore.size();id++){
            btot+=pow(distscore[id]-distscore_f[id],2);
         }
      }
      chi2 << k << "\t" << btot << endl << flush;
      cout << k << "\t" << btot << endl << flush;
      if (btot<btotbest){
         kbest=k;
         fcountbest=fcount;
         btotbest=btot;
      }
      fcount++;
   }
   cout << "Best k is " << kbest << endl;
   chi2.close();


   cout << "Computing plots at best k..." << endl;
   ofstream outf("fitkmeans_site.dat");
   vmot motmix;
   loadmotswnames(motsfile[fcountbest].c_str(),motmix);
   for (ivmot ivm=motmix.begin();ivm!=motmix.end();ivm++){
      ivm->setscorethr2meaninfo();
   }
   for (unsigned int n=1;n<nbspecies;n++){
      for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){

         if (!ivm->matches[n]) continue;

         initscore=scoref(ivm->alignseq[0],mot.matprec);

         double scoretmp=-1e6;
         Motif motevo;
         for (ivmot ivm1=motmix.begin()+1;ivm1!=motmix.end();ivm1++){
            if (scoref(ivm->alignseq[0],ivm1->matprec)>scoretmp){
               scoretmp=scoref(ivm->alignseq[0],ivm1->matprec);
               motevo=*ivm1;
            }
         }

         // previous m initialized to Id
         vvd pm(vdum);
         pm[0][0]=1.;
         pm[1][1]=1.;
         pm[2][2]=1.;
         pm[3][3]=1.;

         // initial probability initialized to TFBS
         vvd Pdisp_f(wdum);
         vvd Pdisp_h(wdum);
         vvd Initprob(wdum);
         for (unsigned int i=0;i<mot.motwidth;i++) 
            Initprob[i][ivm->alignseq[0][i]]=1;
         Pdisp_f=Initprob;
         Pdisp_h=Initprob;


         // HALPERN
         for (unsigned int i=0;i<mot.motwidth;i++){
            Pdisp_h[i]=evolvedist_halpern(Initprob[i],motevo.matfreq[i],distances[n]);
         }

         // FELSEN
         for (unsigned int i=0;i<mot.motwidth;i++){
            Pdisp_f[i]=evolvedist_felsen(Initprob[i],motevo.matfreq[i],distances[n]);
         }



         // SCORES COMPUTATION
         //
         // First, the "true" TFBS
         //

         outf << numtospecies(n) << " ";
         outf << "Data" << " ";
         outf << scoref(ivm->alignseq[n],mot.matprec)-initscore << endl;


         // draw numforstat sites for ICs computation
         vint site(ivm->alignseq[0]);
         for (unsigned int ks=0;ks<numforstat;ks++){

            vint dumsite(mot.motwidth,4);

            // HALPERN
            vint tempsite(dumsite);
            while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
               for (unsigned int j=0;j<motevo.motwidth;j++){
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
            outf << scoref(site,mot.matprec)-initscore << endl;

            // FELSEN
            tempsite=dumsite;
            while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
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
            outf << scoref(site,mot.matprec)-initscore << endl;

         }
         //         getchar();
      }
   }
   outf.close();
}
   
   void
fitdistkmeanspersite(Motif & mot,const char * filename)
{
      
   ifstream inf;
   inf.open(filename);
   vstring motsfile;
   back_insert_iterator<vstring> dest(motsfile);
   copy(iisstring(inf),iisstring(),dest);
   inf.close();

   unsigned int numk;
   numk=motsfile.size();

   //number of sites for IC stat:
   unsigned int numforstat(20);

   vd dum(4,0.);
   vd w(dum),pdisp(dum);

   //distances to ref
   double distance;
   vd distances(nbspecies,0);
   vd newdist(nbspecies,0);
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
      distances[0]=0;
      distances[1]=0.1587;
      distances[2]=0.4653;
      distances[3]=0.4662;
      distances[4]=0.4708;
      distances[5]=0.496;
      distances[6]=0.5512;
      distances[7]=0.5512;
      distances[8]=0.5397;
      distances[9]=0.4971;
   }


   evolutionary_model=2;
   inittreedist();
   vvd vdum(4,dum);
   vvd wdum(mot.motwidth,dum);

   double initscore;

   // find best (K,tau) for Felsen model
   cout << "Finding optimal (K,tau)..." << endl;
   double taubest,btotbest(1e6);
   unsigned int kbest;
   ofstream chi2("chi2.dat");
   unsigned int fcountbest(0);
   
   for (double tau=0.2;tau<15;tau+=0.1){
   
      for (unsigned int id=0;id<nbspecies;id++){
         newdist[id]=distances[id]/tau;
      }
   
      unsigned int fcount=0;
   
      for (unsigned int k=2;k<numk;k++){

         vmot motmix;
         loadmotswnames(motsfile[fcount].c_str(),motmix);
         for (ivmot ivm=motmix.begin();ivm!=motmix.end();ivm++){
            ivm->setscorethr2meaninfo();
            ivm->motscorethrcons=ivm->motscorethr2;
         }

         int index=0;
         double btot(0);
         for (unsigned int n=1;n<nbspecies;n++){
            
            vd distscore(80,0.); // between -4 and +4 with step 0.1
            vd distscore_f(80,0.); // between -4 and +4 with step 0.1

            //ivma ivmstop=min(mot.seqs.begin()+numforstat,mot.seqs.end());
            ivma ivmstop=mot.seqs.end();

            for (ivma ivm=mot.seqs.begin();ivm!=ivmstop;ivm++){

               if (!ivm->matches[n]) continue;

               initscore=scoref(ivm->alignseq[0],mot.matprec);

               double scoretmp=-1e6;
               Motif motevo;
               for (ivmot ivm1=motmix.begin()+1;ivm1!=motmix.end();ivm1++){
                  if (scoref(ivm->alignseq[0],ivm1->matprec)>scoretmp){
                     scoretmp=scoref(ivm->alignseq[0],ivm1->matprec);
                     motevo=*ivm1;
                  }
               }

               // previous m initialized to Id
               vvd pm(vdum);
               pm[0][0]=1.;
               pm[1][1]=1.;
               pm[2][2]=1.;
               pm[3][3]=1.;

               // initial probability initialized to TFBS
               vvd Pdisp_f(wdum);
               vvd Initprob(wdum);
               for (unsigned int i=0;i<motevo.motwidth;i++){
                  Initprob[i][ivm->alignseq[0][i]]=1;
               }
               Pdisp_f=Initprob;
             
               vint site(ivm->alignseq[0]);

               // FELSEN
               for (unsigned int i=0;i<mot.motwidth;i++){
                  Pdisp_f[i]=evolvedist_felsen(Initprob[i],motevo.matfreq[i],newdist[n]);
               }

               // SCORES COMPUTATION
               // First, the "true" TFBS
               index=(int)(10*(4+(scoref(ivm->alignseq[n],mot.matprec)-initscore)));
               if (index>=0 && index<80) distscore[index]++;
               // draw numforstat sites for ICs computation
               for (unsigned int ks=0;ks<numforstat;ks++){

                  vint dumsite(motevo.motwidth,4);
                  vint tempsite(dumsite);

                  // FELSEN
                  while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
                     for (unsigned int j=0;j<motevo.motwidth;j++){
                        double pdist[4]={Pdisp_f[j][0],Pdisp_f[j][1],Pdisp_f[j][2],Pdisp_f[j][3]};
                        gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
                        unsigned int base=gsl_ran_discrete (gslran,gsldist);
                        gsl_ran_discrete_free (gsldist);
                        site[j]=base;
                     }
                     tempsite=site;
                  }
                  index=(int)(10*(4+(scoref(site,mot.matprec)-initscore)));
                  if (index>=0 && index<80) distscore_f[index]+=1./numforstat;
               }
            }

            for (unsigned int id=0;id<distscore.size();id++){
               btot+=pow(distscore[id]-distscore_f[id],2);
            }
         }
         chi2 << k << "\t" << tau << "\t" << btot << endl << flush;
         cout << k << "\t" << tau << "\t" << btot << endl << flush;
         if (btot<btotbest){
            kbest=k;
            taubest=tau;
            fcountbest=fcount;
            btotbest=btot;
         }
         fcount++;
      }
   }
   cout << "Best (k,tau) is (" << kbest << "," << taubest << ")" << endl;
   chi2.close();

   numforstat=100;

   cout << "Computing plots at best (k,tau)..." << endl;
   ofstream outf("fitdistkmeans_site.dat");
   for (unsigned int id=0;id<nbspecies;id++){
      newdist[id]=distances[id]/taubest;
   }
   vmot motmix;
   loadmotswnames(motsfile[fcountbest].c_str(),motmix);
   for (ivmot ivm=motmix.begin();ivm!=motmix.end();ivm++){
      ivm->setscorethr2meaninfo();
   }
   for (unsigned int n=1;n<nbspecies;n++){
      for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){

         if (!ivm->matches[n]) continue;

         initscore=scoref(ivm->alignseq[0],mot.matprec);

         double scoretmp=-1e6;
         Motif motevo;
         for (ivmot ivm1=motmix.begin()+1;ivm1!=motmix.end();ivm1++){
            if (scoref(ivm->alignseq[0],ivm1->matprec)>scoretmp){
               scoretmp=scoref(ivm->alignseq[0],ivm1->matprec);
               motevo=*ivm1;
            }
         }

         // previous m initialized to Id
         vvd pm(vdum);
         pm[0][0]=1.;
         pm[1][1]=1.;
         pm[2][2]=1.;
         pm[3][3]=1.;

         // initial probability initialized to TFBS
         vvd Pdisp_f(wdum);
         vvd Pdisp_h(wdum);
         vvd Initprob(wdum);
         for (unsigned int i=0;i<mot.motwidth;i++) 
            Initprob[i][ivm->alignseq[0][i]]=1;
         Pdisp_f=Initprob;
         Pdisp_h=Initprob;


         // HALPERN
         for (unsigned int i=0;i<mot.motwidth;i++){
            Pdisp_h[i]=evolvedist_halpern(Initprob[i],motevo.matfreq[i],newdist[n]);
         }

         // FELSEN
         for (unsigned int i=0;i<mot.motwidth;i++){
            Pdisp_f[i]=evolvedist_felsen(Initprob[i],motevo.matfreq[i],newdist[n]);
         }



         // SCORES COMPUTATION
         //
         // First, the "true" TFBS
         //

         outf << numtospecies(n) << " ";
         outf << "Data" << " ";
         outf << scoref(ivm->alignseq[n],mot.matprec)-initscore << endl;


         // draw numforstat sites for ICs computation
         vint site(ivm->alignseq[0]);
         for (unsigned int ks=0;ks<numforstat;ks++){

            vint dumsite(mot.motwidth,4);

            // HALPERN
            vint tempsite(dumsite);
            while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
               for (unsigned int j=0;j<motevo.motwidth;j++){
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
            outf << scoref(site,mot.matprec)-initscore << endl;

            // FELSEN
            tempsite=dumsite;
            while (scoref(tempsite,mot.matprec)<mot.motscorethr2){
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
            outf << scoref(site,mot.matprec)-initscore << endl;

         }
         //         getchar();
      }
   }
   outf.close();
}

/*    void
 * fitscorepersite(Motif & mot)
 * {
 * 
 *    unsigned int numhamm=2;
 * 
 *    //number of sites for IC stat:
 *    unsigned int numforstat(50);
 * 
 *    double time_step=0.01; // for barrier
 *    vd dum(4,0.);
 *    vd w(dum),pdisp(dum);
 * 
 *    //distances to ref
 *    double distance;
 *    vd distances(nbspecies,0);
 *    if (species=="droso"){
 *       distances[0]=0;
 *       distances[1]=0.0968;
 *       distances[2]=0.1008;
 *       distances[3]=0.2334;
 *       distances[4]=0.2229;
 *       distances[5]=1.5141;
 *       distances[6]=1.2973;
 *       distances[7]=1.3057;
 *       distances[8]=1.5334;
 *       distances[9]=1.4741;
 *       distances[10]=1.5203;
 *       distances[11]=1.4661;
 *    }
 *    else if (species=="eutherian"){
 *       distances[0]=0;
 *       distances[1]=0.1587;
 *       distances[2]=0.4653;
 *       distances[3]=0.4662;
 *       distances[4]=0.4708;
 *       distances[5]=0.496;
 *       distances[6]=0.5512;
 *       distances[7]=0.5512;
 *       distances[8]=0.5397;
 *       distances[9]=0.4971;
 *    }
 * 
 * 
 *    evolutionary_model=2;
 *    inittreedist();
 *    vvd vdum(4,dum);
 *    vvd wdum(mot.motwidth,dum);
 *    vvvd M_h(mot.motwidth,vdum);
 *    vvvd M_f(mot.motwidth,vdum);
 * 
 *    for (unsigned int i=0;i<mot.motwidth;i++){
 *       M_h[i]=instant_rates_halpern (mot.matfreq[i], time_step);
 *       M_f[i]=instant_rates_felsen (mot.matfreq[i], time_step);
 *    }
 * 
 *    double initscore;
 * 
 *    // find best score for Felsen model
 *    cout << "Finding optimal score..." << endl;
 *    double scorebest,btotbest(1e6);
 *    ofstream chi2("chi2.dat");
 *    
 *    double maxinfo(0);
 *    int j(0);
 *    for (ivvd ivv=mot.matprec.begin();ivv!=mot.matprec.end();ivv++){
 *       int i(0);
 *       double maxcol(-10);
 *       for (ivd iv=ivv->begin();iv!=ivv->end();iv++){
 *          if (*iv>maxcol) maxcol=*iv;
 *          i++;
 *       }
 *       j++;
 *       maxinfo+=maxcol;
 *    }
 * 
 *    for (double score=3;score<.9*maxinfo;score+=0.05){
 * 
 *       mot.motscorethr2=score/10*mot.motwidth;
 *       mot.motscorethrcons=mot.motscorethr2;
 *       mot.matinithamming(mot.motscorethr2,numhamm);
 * 
 *       int index=0;
 *       double btot(0);
 *       for (unsigned int n=1;n<nbspecies;n++){
 *          vd distscore(80,0.); // between -4 and +4 with step 0.1
 *          vd distscore_f(80,0.); // between -4 and +4 with step 0.1
 *          //for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){}
 *          for (ivma ivm=mot.seqs.begin();ivm!=min(mot.seqs.begin()+numforstat,mot.seqs.end());ivm++){
 * 
 *             if (!ivm->matches[n]) continue;
 * 
 *             initscore=scoref(ivm->alignseq[0],mot.matprec);
 * 
 *             // previous m initialized to Id
 *             vvd pm(vdum);
 *             pm[0][0]=1.;
 *             pm[1][1]=1.;
 *             pm[2][2]=1.;
 *             pm[3][3]=1.;
 * 
 *             // initial probability initialized to TFBS
 *             vvd Pdisp_f(wdum);
 *             vvd Initprob(wdum);
 *             for (unsigned int i=0;i<mot.motwidth;i++) 
 *                Initprob[i][ivm->alignseq[0][i]]=1;
 *             Pdisp_f=Initprob;
 *             vint site(ivm->alignseq[0]);
 * 
 *             // FELSEN
 *             for (unsigned int i=0;i<mot.motwidth;i++){
 *                Pdisp_f[i]=evolvedist_felsen(Initprob[i],mot.matfreq[i],distances[n]);
 *             }
 * 
 *             // SCORES COMPUTATION
 *             // First, the "true" TFBS
 *             index=(int)(10*(4+(scoref(ivm->alignseq[n],mot.matprec)-initscore)));
 *             if (index>=0 && index<80) distscore[index]++;
 *             // draw numforstat sites for ICs computation
 *             for (unsigned int ks=0;ks<numforstat;ks++){
 * 
 *                vint dumsite(mot.motwidth,4);
 *                vint tempsite(dumsite);
 * 
 *                // FELSEN
 *                tempsite=dumsite;
 *                while (scorefhamming(tempsite,site)>numhamm){
 *                   for (unsigned int j=0;j<mot.motwidth;j++){
 *                      double pdist[4]={Pdisp_f[j][0],Pdisp_f[j][1],Pdisp_f[j][2],Pdisp_f[j][3]};
 *                      gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
 *                      unsigned int base=gsl_ran_discrete (gslran,gsldist);
 *                      gsl_ran_discrete_free (gsldist);
 *                      site[j]=base;
 *                   }
 *                   tempsite=site;
 *                }
 *                index=(int)(10*(4+(scoref(site,mot.matprec)-initscore)));
 *                if (index>=0 && index<80) distscore_f[index]+=1./numforstat;
 * 
 *             }
 *          }
 * 
 *          for (unsigned int id=0;id!=distscore.size();id++){
 *             btot+=pow(distscore[id]-distscore_f[id],2)/min(numforstat,mot.seqs.size());
 *          }
 *       }
 *       chi2 << score << "\t" << btot << endl;
 *       cout << score << "\t" << btot << endl;
 *       if (btot<btotbest && scorebest<10){
 *          scorebest=score;
 *          btotbest=btot;
 *       }
 *    }
 *    cout << "Best score is " << scorebest << endl;
 *    chi2.close();
 * 
 * 
 *    scorebest=10;
 * 
 *    mot.motscorethr2=scorebest/10*mot.motwidth;
 *    mot.motscorethrcons=mot.motscorethr2;
 *    mot.matinithamming(mot.motscorethr2,numhamm);
 * 
 *    if (mot.seqs.size()<10){
 *       cout << "No sequences! Lowering threshold..." << endl;
 *       while (mot.seqs.size()<10){
 *          scorebest-=0.05;
 *          mot.motscorethr2=scorebest/10*mot.motwidth;
 *          mot.motscorethrcons=mot.motscorethr2;
 *          mot.matinithamming(mot.motscorethr2,numhamm);
 *       }
 * 
 *    }
 *    cout << "Best score is " << scorebest << endl;
 *    cout << "There are " << mot.seqs.size() << " sequences " << endl;
 *    // mot.compprec();
 *    //
 * 
 *    for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){
 *       ivm->print();
 *    }
 * 
 *    cout << "Computing plots at best score..." << endl;
 *    ofstream outf("fitscore_site.dat");
 *    for (unsigned int n=1;n<nbspecies;n++){
 *       //for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){}
 *       for (ivma ivm=mot.seqs.begin();ivm!=min(mot.seqs.begin()+numforstat,mot.seqs.end());ivm++){
 * 
 *          if (!ivm->matches[n]) continue;
 * 
 *          initscore=scoref(ivm->alignseq[0],mot.matprec);
 * 
 *          // previous m initialized to Id
 *          vvd pm(vdum);
 *          pm[0][0]=1.;
 *          pm[1][1]=1.;
 *          pm[2][2]=1.;
 *          pm[3][3]=1.;
 * 
 *          // initial probability initialized to TFBS
 *          vvd Pdisp_f(wdum);
 *          vvd Pdisp_h(wdum);
 *          vvd Initprob(wdum);
 *          for (unsigned int i=0;i<mot.motwidth;i++) 
 *             Initprob[i][ivm->alignseq[0][i]]=1;
 *          Pdisp_f=Initprob;
 *          Pdisp_h=Initprob;
 * 
 * 
 *          // HALPERN
 *          for (unsigned int i=0;i<mot.motwidth;i++){
 *             Pdisp_h[i]=evolvedist_halpern(Initprob[i],mot.matfreq[i],distances[n]);
 *          }
 * 
 *          // FELSEN
 *          for (unsigned int i=0;i<mot.motwidth;i++){
 *             Pdisp_f[i]=evolvedist_felsen(Initprob[i],mot.matfreq[i],distances[n]);
 *          }
 * 
 * 
 *          // SCORES COMPUTATION
 *          //
 *          // First, the "true" TFBS
 *          //
 * 
 *          outf << numtospecies(n) << " ";
 *          outf << "Data" << " ";
 *          outf << scoref(ivm->alignseq[n],mot.matprec)-initscore << endl;
 * 
 * 
 *          // draw numforstat sites for ICs computation
 *          //cout << "data " << vinttostring(ivm->alignseq[n]) << endl;
 *          for (unsigned int ks=0;ks<numforstat;ks++){
 * 
 *             vint dumsite(mot.motwidth,4);
 *             vint site(ivm->alignseq[0]);
 * 
 *             // HALPERN
 *             vint tempsite(dumsite);
 *             while (scorefhamming(tempsite,site)>numhamm){
 *                for (unsigned int j=0;j<mot.motwidth;j++){
 *                   double pdist[4]={Pdisp_h[j][0],Pdisp_h[j][1],Pdisp_h[j][2],Pdisp_h[j][3]};
 *                   gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
 *                   unsigned int base=gsl_ran_discrete (gslran,gsldist);
 *                   gsl_ran_discrete_free (gsldist);
 *                   site[j]=base;
 *                }
 *                tempsite=site;
 *             }
 *             outf << numtospecies(n) << " ";
 *             outf << "Halpern" << " ";
 *             outf << scoref(site,mot.matprec)-initscore << endl;
 * 
 *             //cout << "halpern " << vinttostring(site) << endl;
 * 
 *             // FELSEN
 *             tempsite=dumsite;
 *             while (scorefhamming(tempsite,site)>numhamm){
 *                for (unsigned int j=0;j<mot.motwidth;j++){
 *                   double pdist[4]={Pdisp_f[j][0],Pdisp_f[j][1],Pdisp_f[j][2],Pdisp_f[j][3]};
 *                   gsl_ran_discrete_t * gsldist=gsl_ran_discrete_preproc (4,pdist);
 *                   unsigned int base=gsl_ran_discrete (gslran,gsldist);
 *                   gsl_ran_discrete_free (gsldist);
 *                   site[j]=base;
 *                }
 *                tempsite=site;
 *             }
 *             outf << numtospecies(n) << " ";
 *             outf << "Felsen" << " ";
 *             outf << scoref(site,mot.matprec)-initscore << endl;
 * 
 * 
 * 
 *          }
 *          //         getchar();
 *       }
 *    }
 *    outf.close();
 * }
 */


   void   // Currently used phylogenetic tree for drosophilae :  Heger and Pontig, 2007
inittreedist()
{
   pa=conca;
   pc=concc;

   treedist.clear();
   if (species=="droso"){
      treedist.push_back(noeud(1,2,12,0.017,0.021)); // 12
      treedist.push_back(noeud(0,12,13,0.0518,0.028)); // 13
      treedist.push_back(noeud(3,4,14,0.0885,0.0780)); // 14
      treedist.push_back(noeud(13,14,15,0.0623,0.0308)); // 15
      treedist.push_back(noeud(5,15,16,0.82,0.58)); // 16
      treedist.push_back(noeud(6,7,17,0.0056,0.0140)); // 17
      treedist.push_back(noeud(16,17,18,0.11,0.4876)); // 18
      treedist.push_back(noeud(9,10,19,0.3118,0.358)); // 19
      treedist.push_back(noeud(11,19,20,0.3598,0.056)); // 20
      treedist.push_back(noeud(8,18,21,0.6513,0.078)); // 21
      treedist.push_back(noeud(20,21,22,0.2092,0.0150)); // 22

      noemax=22;
   } 
   else if (species=="eutherian"){
      // Arbre de ensembl epo 10 eutharian  *** To be translated
      treedist.push_back(noeud(0,1,10,0.0770,0.0817)); // 10 mus and rat
      treedist.push_back(noeud(2,3,11,0.0067,0.0076)); // 11 hom and pan
      treedist.push_back(noeud(4,11,12,0.0220,0.0098)); // 12 pon and 11
      treedist.push_back(noeud(5,12,13,0.0593,0.0121)); // 13 mac and 12
      treedist.push_back(noeud(10,13,14,0.2526,0.1072)); // 14 10 and 13
      treedist.push_back(noeud(6,7,15,0.0796,0.0796)); // 15 bos and sus
      treedist.push_back(noeud(8,15,16,0.1477,0.0796)); // 16 can and 15
      treedist.push_back(noeud(9,16,17,0.1100,0.0049)); // 17 equ and 16
      treedist.push_back(noeud(14,17,18,0.0230,0.0345)); // 18 14 and 17

      noemax=18;
   }

   vtransi.clear();
   if (evolutionary_model==2){
      for (unsigned int i=0;i<noemax+1;i++){
         gsl_matrix * m=gsl_matrix_alloc(4,4);
         vtransi.push_back(m);
      }
      instrates = gsl_matrix_calloc(4,4);
      pnoe1out = gsl_vector_calloc(4);
      pnoe2out = gsl_vector_calloc(4);

      id = gsl_matrix_calloc(4,4);
      for (unsigned int i=0;i<4;i++){
         gsl_matrix_set(id,i,i,1);
      }

      m1=gsl_matrix_calloc(4,4);
      m2=gsl_matrix_calloc(4,4);
      m3=gsl_matrix_calloc(4,4);
      m4=gsl_matrix_calloc(4,4);
      pij=gsl_matrix_calloc(4,4);
      pijp=gsl_matrix_calloc(4,4);
   } 
   else if (evolutionary_model==1){
      proba2=gsl_vector_alloc(4);
   }


   vd distca(nbspecies,0);
   //here prox is defined as a distance
   for (int i=0;i<nbspecies;i++){
      int j=i;
      while (j!=noemax){
         for (ivnoe iv=treedist.begin();iv!=treedist.end();iv++){
            int esp1=iv->esp1;
            int esp2=iv->esp2;
            int noe=iv->noe;
            double prox1=iv->prox1;
            double prox2=iv->prox2;
            if (evolutionary_model==1){
               double correction=0.5+4.0*conca*concc;
               prox1=-log(prox1)*correction;
               prox2=-log(prox2)*correction;
            }
            if (j==esp1){
               distca[i]+=prox1;
               j=noe;
               break;
            }
            else if (j==esp2){
               distca[i]+=prox2;
               j=noe;
               break;
            }
         }
      }
   }

   dlca=distca;
}

   int
speciestonum(string name)//drosonum
{
   if (species=="droso"){
      if (name=="DroMel") return 0;
      else if (name=="DroSim") return 1;
      else if (name=="DroSec") return 2;
      else if (name=="DroYak") return 3;
      else if (name=="DroEre") return 4;
      else if (name=="DroAna") return 5;
      else if (name=="DroPse") return 6;
      else if (name=="DroPer") return 7;
      else if (name=="DroWil") return 8;
      else if (name=="DroVir") return 9;
      else if (name=="DroMoj") return 10;
      else if (name=="DroGri") return 11;
      else return -1;
   }
   else if (species=="eutherian"){
      if (name == "MusMus") return 0;
      else if (name == "RatNor") return 1;
      else if (name == "HomSap") return 2;
      else if (name == "PanTro") return 3;
      else if (name == "PonPyg") return 4;
      else if (name == "MacMul") return 5;
      else if (name == "BosTau") return 6;
      else if (name == "SusScr") return 7;
      else if (name == "CanFam") return 8;
      else if (name == "EquCab") return 9;
      else return -1;
   }
}

string numtospecies(int num)//numdroso
{
   if (species=="droso"){
      if (num==0) return "DroMel";
      else if (num==1) return "DroSim";
      else if (num==2) return "DroSec";
      else if (num==3) return "DroYak";
      else if (num==4) return "DroEre";
      else if (num==5) return "DroAna";
      else if (num==6) return "DroPse";
      else if (num==7) return "DroPer";
      else if (num==8) return "DroWil";
      else if (num==9) return "DroVir";
      else if (num==10) return "DroMoj";
      else if (num==11) return "DroGri";
      else return "No name";
   }
   else if (species=="eutherian"){
      if (num == 0) return "MusMus";
      else if (num == 1) return "RatNor";
      else if (num == 2) return "HomSap";
      else if (num == 3) return "PanTro";
      else if (num == 4) return "PonPyg";
      else if (num == 5) return "MacMul";
      else if (num == 6) return "BosTau";
      else if (num == 7) return "SusScr";
      else if (num == 8) return "CanFam";
      else if (num == 9) return "EquCab";
      else return "No name";
   }

}

   int
instant_rates (const gsl_vector * w, gsl_matrix * rates)
{
   double w0=gsl_vector_get(w,0);
   double w1=gsl_vector_get(w,1);
   double w2=gsl_vector_get(w,2);
   double w3=1-w0-w1-w2;
   if (w0<0 || w1<0 || w2<0 || w3<0){
      return -1;
   }

   double fat=proba_fixation_rel(w1/w0);
   double fac=proba_fixation_rel(pa*w2/pc/w0);
   double fag=proba_fixation_rel(pa*w3/pc/w0);

   double fta=proba_fixation_rel(w0/w1);
   double ftc=proba_fixation_rel(pa*w2/pc/w1);
   double ftg=proba_fixation_rel(pa*w3/pc/w1);

   double fca=proba_fixation_rel(pc*w0/pa/w2);
   double fct=proba_fixation_rel(pc*w1/pa/w2);
   double fcg=proba_fixation_rel(w3/w2);

   double fga=proba_fixation_rel(pc*w0/pa/w3);
   double fgt=proba_fixation_rel(pc*w1/pa/w3);
   double fgc=proba_fixation_rel(w2/w3);

   double prefact=1.0/(4*kappa*pa*pc+0.5);

   gsl_matrix_set(m1, 0, 0, -(pa*fat+pc*fac+pc*kappa*fag));
   gsl_matrix_set(m1, 0, 1, pa*fta);
   gsl_matrix_set(m1, 0, 2, pa*fca);
   gsl_matrix_set(m1, 0, 3, pa*kappa*fga);

   gsl_matrix_set(m1, 1, 0, pa*fat);
   gsl_matrix_set(m1, 1, 1, -(pa*fta+pc*kappa*ftc+pc*ftg));
   gsl_matrix_set(m1, 1, 2, pa*kappa*fct);
   gsl_matrix_set(m1, 1, 3, pa*fgt);

   gsl_matrix_set(m1, 2, 0, pc*fac);
   gsl_matrix_set(m1, 2, 1, pc*kappa*ftc);
   gsl_matrix_set(m1, 2, 2, -(pa*fca+pa*kappa*fct+pc*fcg));
   gsl_matrix_set(m1, 2, 3, pc*fgc);

   gsl_matrix_set(m1, 3, 0, pc*kappa*fag);
   gsl_matrix_set(m1, 3, 1, pc*ftg);
   gsl_matrix_set(m1, 3, 2, pc*fcg);
   gsl_matrix_set(m1, 3, 3, -(pa*kappa*fga+pa*fgt+pc*fgc));

   //   if (species==1) integr_step=0.01;
   //   else if (species==2) integr_step=0.001;
   gsl_matrix_scale(m1,prefact*integr_step);
   gsl_matrix * mattemp;

   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,0.5,m1,m1,0.0,m2);
   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1/3.0,m1,m2,0.0,m3);
   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1/4.0,m1,m3,0.0,m4);

   //Matrice d'évolution Runge Kutta 4
   // Mat(RG4) = Id + h*M + h^2/2*M^2 + h^3/6*M^3 + h^4/24*M^4

   gsl_matrix_memcpy(rates,id);
   gsl_matrix_add(rates,m1);
   gsl_matrix_add(rates,m2);
   gsl_matrix_add(rates,m3);
   gsl_matrix_add(rates,m4);

   return 0;
}

   int
transi (double d, gsl_matrix * rates, gsl_matrix * transitions)
{
   const gsl_odeiv_step_type * T 
      = gsl_odeiv_step_rk8pd;

   gsl_odeiv_step * s 
      = gsl_odeiv_step_alloc (T, 4);
   gsl_odeiv_control * c 
      = gsl_odeiv_control_y_new (1e-3, 0.0);
   gsl_odeiv_evolve * e 
      = gsl_odeiv_evolve_alloc (4);

   gsl_odeiv_system sys = {func, jac, 4, rates};

   double t = 0.0;
   double h = 1e-2;

   double betat=d/2.0/(pa*pa+pc*pc+8.0*pa*pc);


   double p[4] = { 1.0, 0.0 , 0.0, 0.0};
   while (t < betat)
   {
      int status = gsl_odeiv_evolve_apply (e, c, s,
            &sys, 
            &t, betat,
            &h, p);
      if (status != GSL_SUCCESS)
         break;
   }
   gsl_matrix_set(transitions,0,0,p[0]);
   gsl_matrix_set(transitions,1,0,p[1]);
   gsl_matrix_set(transitions,2,0,p[2]);
   gsl_matrix_set(transitions,3,0,p[3]);

   p[0]=0;
   p[1]=1;
   p[2]=0;
   p[3]=0;
   t=0;
   while (t < betat)
   {
      int status = gsl_odeiv_evolve_apply (e, c, s,
            &sys, 
            &t, betat,
            &h, p);
      if (status != GSL_SUCCESS)
         break;
   }
   gsl_matrix_set(transitions,0,1,p[0]);
   gsl_matrix_set(transitions,1,1,p[1]);
   gsl_matrix_set(transitions,2,1,p[2]);
   gsl_matrix_set(transitions,3,1,p[3]);

   p[0]=0;
   p[1]=0;
   p[2]=1;
   p[3]=0;
   t=0;
   while (t < betat)
   {
      int status = gsl_odeiv_evolve_apply (e, c, s,
            &sys, 
            &t, betat,
            &h, p);
      if (status != GSL_SUCCESS)
         break;
   }
   gsl_matrix_set(transitions,0,2,p[0]);
   gsl_matrix_set(transitions,1,2,p[1]);
   gsl_matrix_set(transitions,2,2,p[2]);
   gsl_matrix_set(transitions,3,2,p[3]);

   p[0]=0;
   p[1]=0;
   p[2]=0;
   p[3]=1;
   t=0;
   while (t < betat)
   {
      int status = gsl_odeiv_evolve_apply (e, c, s,
            &sys, 
            &t, betat,
            &h, p);
      if (status != GSL_SUCCESS)
         break;
   }
   gsl_matrix_set(transitions,0,3,p[0]);
   gsl_matrix_set(transitions,1,3,p[1]);
   gsl_matrix_set(transitions,2,3,p[2]);
   gsl_matrix_set(transitions,3,3,p[3]);

   gsl_odeiv_evolve_free (e);
   gsl_odeiv_control_free (c);
   gsl_odeiv_step_free (s);

   return 0;
}		/* -----  end of function transi  ----- */

   int
func (double t, const double y[], double f[],
      void *params)
{
   gsl_matrix * rates=(gsl_matrix *)params;
   gsl_vector_const_view yview=gsl_vector_const_view_array(y,4);
   gsl_vector_view fview=gsl_vector_view_array(f,4);

   gsl_blas_dgemv(CblasNoTrans,1.0,rates,&yview.vector,0.0,&fview.vector);

   return GSL_SUCCESS;
}

   int
jac (double t, const double y[], double *dfdy, 
      double dfdt[], void *params)
{
   gsl_matrix_view dfdy_mat 
      = gsl_matrix_view_array (dfdy, 4, 4);
   gsl_matrix * m = &dfdy_mat.matrix; 

   gsl_matrix * rates=(gsl_matrix *)params;
   gsl_matrix_memcpy(m,rates);

   dfdt[0] = 0.0;
   dfdt[1] = 0.0;
   dfdt[2] = 0.0;
   dfdt[3] = 0.0;
   return GSL_SUCCESS;
}

   double
proba_fixation_rel(double ratio)
{
   if (ratio==1.0) ratio=1.0+1e-6;
   return ratio*log(ratio)/(ratio-1.0);
}

   double
initprobatree(const unsigned int pos,Motalign & ma, gsl_matrix * probatree)
{
   for (unsigned int i=0;i<nbspecies;i++){ 
      if (ma.matches[i]){
         gsl_matrix_set(probatree,0,i,0);
         gsl_matrix_set(probatree,1,i,0);
         gsl_matrix_set(probatree,2,i,0);
         gsl_matrix_set(probatree,3,i,0);
         gsl_matrix_set(probatree,ma.alignseq[i][pos],i,1);
      } else {
         gsl_matrix_set(probatree,0,i,-1);
      }
   }
}


   double
loglikely_column(const unsigned int pos,Motalign & ma, vpgslmat & vtrans, const gsl_vector * w)
{
   gsl_matrix * probatree=gsl_matrix_alloc(4,noemax+1); 
   gsl_vector * wfull=gsl_vector_alloc(4);

   initprobatree(pos,ma,probatree);

   gsl_vector_set(wfull,0,gsl_vector_get(w,0));
   gsl_vector_set(wfull,1,gsl_vector_get(w,1));
   gsl_vector_set(wfull,2,gsl_vector_get(w,2));
   gsl_vector_set(wfull,3,1-gsl_vector_get(w,0)-gsl_vector_get(w,1)-gsl_vector_get(w,2));

   for (ivnoe iv=treedist.begin();iv!=treedist.end();iv++){
      int n1=iv->esp1;
      int n2=iv->esp2;
      if (evolutionary_model==2){
         if (gsl_matrix_get(probatree,0,n1)<-0.5 && gsl_matrix_get(probatree,0,n2)<-0.5){
            gsl_matrix_set(probatree,0,iv->noe,-1.0);
         } else if (gsl_matrix_get(probatree,0,n2)<-0.5){
            gsl_vector_view pnoe1=gsl_matrix_column(probatree,(*iv).esp1);
            gsl_vector_view pnoe=gsl_matrix_column(probatree,(*iv).noe);
            gsl_blas_dgemv(CblasTrans,1.0,vtrans[2*((*iv).noe-nbspecies)],&pnoe1.vector,0.0,&pnoe.vector);
         } else if (gsl_matrix_get(probatree,0,n1)<-0.5){
            gsl_vector_view pnoe2=gsl_matrix_column(probatree,(*iv).esp2);
            gsl_vector_view pnoe=gsl_matrix_column(probatree,(*iv).noe);
            gsl_blas_dgemv(CblasTrans,1.0,vtrans[2*((*iv).noe-nbspecies)+1],&pnoe2.vector,0.0,&pnoe.vector);
         } else {
            gsl_vector_view pnoe1=gsl_matrix_column(probatree,(*iv).esp1);
            gsl_vector_view pnoe2=gsl_matrix_column(probatree,(*iv).esp2);
            gsl_vector_view pnoe=gsl_matrix_column(probatree,(*iv).noe);
            gsl_blas_dgemv(CblasTrans,1.0,vtrans[2*((*iv).noe-nbspecies)],&pnoe1.vector,0.0,pnoe1out);
            gsl_blas_dgemv(CblasTrans,1.0,vtrans[2*((*iv).noe-nbspecies)+1],&pnoe2.vector,0.0,pnoe2out);
            gsl_vector_memcpy(&pnoe.vector,pnoe1out);
            gsl_vector_mul(&pnoe.vector,pnoe2out);
         }
      }
      else if (evolutionary_model=1){
         double prox1=iv->prox1;
         double prox2=iv->prox2;
         if (gsl_matrix_get(probatree,0,n1)<-0.5 && gsl_matrix_get(probatree,0,n2)<-0.5){
            gsl_matrix_set(probatree,0,iv->noe,-1.0);
         } else if (gsl_matrix_get(probatree,0,n2)<-0.5){
            gsl_vector_view pnoe1=gsl_matrix_column(probatree,(*iv).esp1);
            gsl_vector_view pnoe=gsl_matrix_column(probatree,(*iv).noe);
            double sum;
            gsl_blas_ddot(wfull,&pnoe1.vector,&sum);
            sum*=(1-prox1);
            gsl_vector_set_all(&pnoe.vector,sum);
            gsl_blas_daxpy(prox1,&pnoe1.vector,&pnoe.vector);
         } else if (gsl_matrix_get(probatree,0,n1)<-0.5){
            gsl_vector_view pnoe2=gsl_matrix_column(probatree,(*iv).esp2);
            gsl_vector_view pnoe=gsl_matrix_column(probatree,(*iv).noe);
            double sum;
            gsl_blas_ddot(wfull,&pnoe2.vector,&sum);
            sum*=(1-prox2);
            gsl_vector_set_all(&pnoe.vector,sum);
            gsl_blas_daxpy(prox2,&pnoe2.vector,&pnoe.vector);
         } else {
            gsl_vector_view pnoe1=gsl_matrix_column(probatree,iv->esp1);
            gsl_vector_view pnoe2=gsl_matrix_column(probatree,iv->esp2);
            gsl_vector_view pnoe=gsl_matrix_column(probatree,iv->noe);
            double sum1,sum2;
            gsl_blas_ddot(wfull,&pnoe1.vector,&sum1);
            sum1*=(1-prox1);
            gsl_blas_ddot(wfull,&pnoe2.vector,&sum2);
            sum2*=(1-prox2);
            gsl_vector_set_all(&pnoe.vector,sum1);
            gsl_vector_set_all(proba2,sum2);
            gsl_blas_daxpy(prox1,&pnoe1.vector,&pnoe.vector);
            gsl_blas_daxpy(prox2,&pnoe2.vector,proba2);
            gsl_vector_mul(&pnoe.vector,proba2);
         }
      }
   }

   if (gsl_matrix_get(probatree,0,noemax)<-0.5) cout << "error!!!\n"; // 16 instead of 22
   double w0=gsl_vector_get(w,0);
   double w1=gsl_vector_get(w,1);
   double w2=gsl_vector_get(w,2);
   double w3=1-w0-w1-w2;
   double res=log(w0*gsl_matrix_get(probatree,0,noemax)+ //!!
         w1*gsl_matrix_get(probatree,1,noemax)+
         w2*gsl_matrix_get(probatree,2,noemax)+
         w3*gsl_matrix_get(probatree,3,noemax));

   gsl_matrix_free(probatree);
   gsl_vector_free(wfull);

   return res;
}

   double
loglikely(const gsl_vector *w, void *params)
{
   void **par=(void **) params;
   Motif & mot=*((Motif *)(par[0]));
   const unsigned int pos =*((const unsigned int *)(par[1]));

   gsl_matrix * pmattemp;

   if (evolutionary_model==2){

      if (instant_rates(w,instrates)) return 1e10;

      gsl_matrix_memcpy(pij,id);

      if (species=="droso"){
         for (unsigned int i=1;i<117;i++){
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
            pmattemp=pij;
            pij=pijp;
            pijp=pmattemp;
            //      gsl_matrix_memcpy(pij,pijp);
            if (i==1){
               gsl_matrix_memcpy(vtransi[10],pij);
            } else if (i==2){
               gsl_matrix_memcpy(vtransi[0],pij);
               gsl_matrix_memcpy(vtransi[11],pij);
            } else if (i==3){
               gsl_matrix_memcpy(vtransi[1],pij);
               gsl_matrix_memcpy(vtransi[3],pij);
            } else if (i==4){
               gsl_matrix_memcpy(vtransi[7],pij);
            } else if (i==7){
               gsl_matrix_memcpy(vtransi[2],pij);
               gsl_matrix_memcpy(vtransi[6],pij);
            } else if (i==10){
               gsl_matrix_memcpy(vtransi[5],pij);
               gsl_matrix_memcpy(vtransi[17],pij);
            } else if (i==11){
               gsl_matrix_memcpy(vtransi[4],pij);
            } else if (i==14){
               gsl_matrix_memcpy(vtransi[20],pij);
            } else if (i==15){
               gsl_matrix_memcpy(vtransi[21],pij);
            } else if (i==35){
               gsl_matrix_memcpy(vtransi[12],pij);
               gsl_matrix_memcpy(vtransi[14],pij);
            } else if (i==42){
               gsl_matrix_memcpy(vtransi[19],pij);
            } else if (i==46){
               gsl_matrix_memcpy(vtransi[16],pij);
            } else if (i==49){
               gsl_matrix_memcpy(vtransi[15],pij);
            } else if (i==58){
               gsl_matrix_memcpy(vtransi[9],pij);
            } else if (i==68){
               gsl_matrix_memcpy(vtransi[13],pij);
            } else if (i==82){
               gsl_matrix_memcpy(vtransi[8],pij);
            } else if (i==116){
               gsl_matrix_memcpy(vtransi[18],pij);
            }
         }
      }
      else if (species=="eutherian"){
         //         for (unsigned int i=1;i<254;i++){
         //            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
         //            pmattemp=pij;
         //            pij=pijp;
         //            pijp=pmattemp;
         //            //      gsl_matrix_memcpy(pij,pijp);
         //            //      good approx, for integr_step=0.001
         //            if (i==5){
         //               gsl_matrix_memcpy(vtransi[15],pij);
         //            } else if (i==7){
         //               gsl_matrix_memcpy(vtransi[2],pij);
         //            } else if (i==8){
         //               gsl_matrix_memcpy(vtransi[3],pij);
         //            } else if (i==10){
         //               gsl_matrix_memcpy(vtransi[5],pij);
         //            } else if (i==12){
         //               gsl_matrix_memcpy(vtransi[7],pij);
         //            } else if (i==22){
         //               gsl_matrix_memcpy(vtransi[4],pij);
         //            } else if (i==23){
         //               gsl_matrix_memcpy(vtransi[16],pij);
         //            } else if (i==35){
         //               gsl_matrix_memcpy(vtransi[17],pij);
         //            } else if (i==59){
         //               gsl_matrix_memcpy(vtransi[6],pij);
         //            } else if (i==77){
         //               gsl_matrix_memcpy(vtransi[0],pij);
         //            } else if (i==80){
         //               gsl_matrix_memcpy(vtransi[10],pij);
         //               gsl_matrix_memcpy(vtransi[11],pij);
         //               gsl_matrix_memcpy(vtransi[13],pij);
         //            } else if (i==82){
         //               gsl_matrix_memcpy(vtransi[1],pij);
         //            } else if (i==107){
         //               gsl_matrix_memcpy(vtransi[9],pij);
         //            } else if (i==110){
         //               gsl_matrix_memcpy(vtransi[14],pij);
         //            } else if (i==148){
         //               gsl_matrix_memcpy(vtransi[12],pij);
         //            } else if (i==253){
         //               gsl_matrix_memcpy(vtransi[8],pij);
         //            }  		
         //            }
         // approx, integr_step=0.01
         for (unsigned int i=1;i<26;i++){
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
            pmattemp=pij;
            pij=pijp;
            pijp=pmattemp;
            //      gsl_matrix_memcpy(pij,pijp);
            if (i==1){
               gsl_matrix_memcpy(vtransi[15],pij);
               gsl_matrix_memcpy(vtransi[2],pij);
               gsl_matrix_memcpy(vtransi[3],pij);
               gsl_matrix_memcpy(vtransi[5],pij);
               gsl_matrix_memcpy(vtransi[7],pij);
            } else if (i==2){
               gsl_matrix_memcpy(vtransi[4],pij);
               gsl_matrix_memcpy(vtransi[16],pij);
            } else if (i==3){
               gsl_matrix_memcpy(vtransi[17],pij);
            } else if (i==6){
               gsl_matrix_memcpy(vtransi[6],pij);
            } else if (i==8){
               gsl_matrix_memcpy(vtransi[0],pij);
               gsl_matrix_memcpy(vtransi[10],pij);
               gsl_matrix_memcpy(vtransi[11],pij);
               gsl_matrix_memcpy(vtransi[13],pij);
               gsl_matrix_memcpy(vtransi[1],pij);
            } else if (i==11){
               gsl_matrix_memcpy(vtransi[9],pij);
               gsl_matrix_memcpy(vtransi[14],pij);
            } else if (i==15){
               gsl_matrix_memcpy(vtransi[12],pij);
            } else if (i==25){
               gsl_matrix_memcpy(vtransi[8],pij);
            }
         }
      }
   }
   //   for (unsigned int i=0;i<4;i++){
   //      for (unsigned int j=0;j<4;j++){
   //         cout << gsl_matrix_get(vtransi[8],i,j) << " ";
   //      }
   //      cout << endl;
   //   }
   //   cout << endl;

   double logli=0;
   for (ivma ima=mot.seqs.begin();ima!=mot.seqs.end();ima++){
      logli+=loglikely_column(pos,*ima,vtransi,w);
   }

   double w0=gsl_vector_get(w,0);
   double w1=gsl_vector_get(w,1);
   double w2=gsl_vector_get(w,2);
   double w3=1-w0-w1-w2;
   logli+=alpha*(log(w0)+log(w1))+beta*(log(w2)+log(w3));
   // logli+=(alpha-1.)*(log(w0)+log(w1))+(beta-1.)*(log(w2)+log(w3));
   

   return -logli;
}

   double
likelyhood(vd x, void *params)
{
   gsl_vector *w;
   w = gsl_vector_alloc (3);
   gsl_vector_set (w, 0, x[0]);
   gsl_vector_set (w, 1, x[1]);
   gsl_vector_set (w, 2, x[2]);
   void **par=(void **) params;
   Motif & mot=*((Motif *)(par[0]));
   const unsigned int pos =*((const unsigned int *)(par[1]));

   gsl_matrix * pmattemp;

   if (evolutionary_model==2){

      if (instant_rates(w,instrates)) return 1e10;

      gsl_matrix_memcpy(pij,id);

      if (species=="droso"){
         for (unsigned int i=1;i<117;i++){
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
            pmattemp=pij;
            pij=pijp;
            pijp=pmattemp;
            //      gsl_matrix_memcpy(pij,pijp);
            if (i==1){
               gsl_matrix_memcpy(vtransi[10],pij);
            } else if (i==2){
               gsl_matrix_memcpy(vtransi[0],pij);
               gsl_matrix_memcpy(vtransi[11],pij);
            } else if (i==3){
               gsl_matrix_memcpy(vtransi[1],pij);
               gsl_matrix_memcpy(vtransi[3],pij);
            } else if (i==4){
               gsl_matrix_memcpy(vtransi[7],pij);
            } else if (i==7){
               gsl_matrix_memcpy(vtransi[2],pij);
               gsl_matrix_memcpy(vtransi[6],pij);
            } else if (i==10){
               gsl_matrix_memcpy(vtransi[5],pij);
               gsl_matrix_memcpy(vtransi[17],pij);
            } else if (i==11){
               gsl_matrix_memcpy(vtransi[4],pij);
            } else if (i==14){
               gsl_matrix_memcpy(vtransi[20],pij);
            } else if (i==15){
               gsl_matrix_memcpy(vtransi[21],pij);
            } else if (i==35){
               gsl_matrix_memcpy(vtransi[12],pij);
               gsl_matrix_memcpy(vtransi[14],pij);
            } else if (i==42){
               gsl_matrix_memcpy(vtransi[19],pij);
            } else if (i==46){
               gsl_matrix_memcpy(vtransi[16],pij);
            } else if (i==49){
               gsl_matrix_memcpy(vtransi[15],pij);
            } else if (i==58){
               gsl_matrix_memcpy(vtransi[9],pij);
            } else if (i==68){
               gsl_matrix_memcpy(vtransi[13],pij);
            } else if (i==82){
               gsl_matrix_memcpy(vtransi[8],pij);
            } else if (i==116){
               gsl_matrix_memcpy(vtransi[18],pij);
            }
         }
      }
      else if (species=="eutherian"){
         //         for (unsigned int i=1;i<254;i++){
         //            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
         //            pmattemp=pij;
         //            pij=pijp;
         //            pijp=pmattemp;
         //            //      gsl_matrix_memcpy(pij,pijp);
         //            //      good approx, for integr_step=0.001
         //            if (i==5){
         //               gsl_matrix_memcpy(vtransi[15],pij);
         //            } else if (i==7){
         //               gsl_matrix_memcpy(vtransi[2],pij);
         //            } else if (i==8){
         //               gsl_matrix_memcpy(vtransi[3],pij);
         //            } else if (i==10){
         //               gsl_matrix_memcpy(vtransi[5],pij);
         //            } else if (i==12){
         //               gsl_matrix_memcpy(vtransi[7],pij);
         //            } else if (i==22){
         //               gsl_matrix_memcpy(vtransi[4],pij);
         //            } else if (i==23){
         //               gsl_matrix_memcpy(vtransi[16],pij);
         //            } else if (i==35){
         //               gsl_matrix_memcpy(vtransi[17],pij);
         //            } else if (i==59){
         //               gsl_matrix_memcpy(vtransi[6],pij);
         //            } else if (i==77){
         //               gsl_matrix_memcpy(vtransi[0],pij);
         //            } else if (i==80){
         //               gsl_matrix_memcpy(vtransi[10],pij);
         //               gsl_matrix_memcpy(vtransi[11],pij);
         //               gsl_matrix_memcpy(vtransi[13],pij);
         //            } else if (i==82){
         //               gsl_matrix_memcpy(vtransi[1],pij);
         //            } else if (i==107){
         //               gsl_matrix_memcpy(vtransi[9],pij);
         //            } else if (i==110){
         //               gsl_matrix_memcpy(vtransi[14],pij);
         //            } else if (i==148){
         //               gsl_matrix_memcpy(vtransi[12],pij);
         //            } else if (i==253){
         //               gsl_matrix_memcpy(vtransi[8],pij);
         //            }  		
         //            }
         // approx, integr_step=0.01
         for (unsigned int i=1;i<26;i++){
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
            pmattemp=pij;
            pij=pijp;
            pijp=pmattemp;
            //      gsl_matrix_memcpy(pij,pijp);
            if (i==1){
               gsl_matrix_memcpy(vtransi[15],pij);
               gsl_matrix_memcpy(vtransi[2],pij);
               gsl_matrix_memcpy(vtransi[3],pij);
               gsl_matrix_memcpy(vtransi[5],pij);
               gsl_matrix_memcpy(vtransi[7],pij);
            } else if (i==2){
               gsl_matrix_memcpy(vtransi[4],pij);
               gsl_matrix_memcpy(vtransi[16],pij);
            } else if (i==3){
               gsl_matrix_memcpy(vtransi[17],pij);
            } else if (i==6){
               gsl_matrix_memcpy(vtransi[6],pij);
            } else if (i==8){
               gsl_matrix_memcpy(vtransi[0],pij);
               gsl_matrix_memcpy(vtransi[10],pij);
               gsl_matrix_memcpy(vtransi[11],pij);
               gsl_matrix_memcpy(vtransi[13],pij);
               gsl_matrix_memcpy(vtransi[1],pij);
            } else if (i==11){
               gsl_matrix_memcpy(vtransi[9],pij);
               gsl_matrix_memcpy(vtransi[14],pij);
            } else if (i==15){
               gsl_matrix_memcpy(vtransi[12],pij);
            } else if (i==25){
               gsl_matrix_memcpy(vtransi[8],pij);
            }
         }
      }
   }
   //   for (unsigned int i=0;i<4;i++){
   //      for (unsigned int j=0;j<4;j++){
   //         cout << gsl_matrix_get(vtransi[8],i,j) << " ";
   //      }
   //      cout << endl;
   //   }
   //   cout << endl;

   double logli=0;
   for (ivma ima=mot.seqs.begin();ima!=mot.seqs.end();ima++){
      logli+=loglikely_column(pos,*ima,vtransi,w);
   }
   double w0=gsl_vector_get(w,0);
   double w1=gsl_vector_get(w,1);
   double w2=gsl_vector_get(w,2);
   double w3=1-w0-w1-w2;
   logli+=(alpha-1.)*(log(w0)+log(w1))+(beta-1.)*(log(w2)+log(w3));

   return exp(logli);
}

   double
loglikelyhood(vd x, void *params)
{
   gsl_vector *w;
   w = gsl_vector_alloc (3);
   gsl_vector_set (w, 0, x[0]);
   gsl_vector_set (w, 1, x[1]);
   gsl_vector_set (w, 2, x[2]);
   void **par=(void **) params;
   Motif & mot=*((Motif *)(par[0]));
   const unsigned int pos =*((const unsigned int *)(par[1]));

   gsl_matrix * pmattemp;

   if (evolutionary_model==2){

      if (instant_rates(w,instrates)) return 1e10;

      gsl_matrix_memcpy(pij,id);

      if (species=="droso"){
         for (unsigned int i=1;i<117;i++){
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
            pmattemp=pij;
            pij=pijp;
            pijp=pmattemp;
            //      gsl_matrix_memcpy(pij,pijp);
            if (i==1){
               gsl_matrix_memcpy(vtransi[10],pij);
            } else if (i==2){
               gsl_matrix_memcpy(vtransi[0],pij);
               gsl_matrix_memcpy(vtransi[11],pij);
            } else if (i==3){
               gsl_matrix_memcpy(vtransi[1],pij);
               gsl_matrix_memcpy(vtransi[3],pij);
            } else if (i==4){
               gsl_matrix_memcpy(vtransi[7],pij);
            } else if (i==7){
               gsl_matrix_memcpy(vtransi[2],pij);
               gsl_matrix_memcpy(vtransi[6],pij);
            } else if (i==10){
               gsl_matrix_memcpy(vtransi[5],pij);
               gsl_matrix_memcpy(vtransi[17],pij);
            } else if (i==11){
               gsl_matrix_memcpy(vtransi[4],pij);
            } else if (i==14){
               gsl_matrix_memcpy(vtransi[20],pij);
            } else if (i==15){
               gsl_matrix_memcpy(vtransi[21],pij);
            } else if (i==35){
               gsl_matrix_memcpy(vtransi[12],pij);
               gsl_matrix_memcpy(vtransi[14],pij);
            } else if (i==42){
               gsl_matrix_memcpy(vtransi[19],pij);
            } else if (i==46){
               gsl_matrix_memcpy(vtransi[16],pij);
            } else if (i==49){
               gsl_matrix_memcpy(vtransi[15],pij);
            } else if (i==58){
               gsl_matrix_memcpy(vtransi[9],pij);
            } else if (i==68){
               gsl_matrix_memcpy(vtransi[13],pij);
            } else if (i==82){
               gsl_matrix_memcpy(vtransi[8],pij);
            } else if (i==116){
               gsl_matrix_memcpy(vtransi[18],pij);
            }
         }
      }
      else if (species=="eutherian"){
         //         for (unsigned int i=1;i<254;i++){
         //            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
         //            pmattemp=pij;
         //            pij=pijp;
         //            pijp=pmattemp;
         //            //      gsl_matrix_memcpy(pij,pijp);
         //            //      good approx, for integr_step=0.001
         //            if (i==5){
         //               gsl_matrix_memcpy(vtransi[15],pij);
         //            } else if (i==7){
         //               gsl_matrix_memcpy(vtransi[2],pij);
         //            } else if (i==8){
         //               gsl_matrix_memcpy(vtransi[3],pij);
         //            } else if (i==10){
         //               gsl_matrix_memcpy(vtransi[5],pij);
         //            } else if (i==12){
         //               gsl_matrix_memcpy(vtransi[7],pij);
         //            } else if (i==22){
         //               gsl_matrix_memcpy(vtransi[4],pij);
         //            } else if (i==23){
         //               gsl_matrix_memcpy(vtransi[16],pij);
         //            } else if (i==35){
         //               gsl_matrix_memcpy(vtransi[17],pij);
         //            } else if (i==59){
         //               gsl_matrix_memcpy(vtransi[6],pij);
         //            } else if (i==77){
         //               gsl_matrix_memcpy(vtransi[0],pij);
         //            } else if (i==80){
         //               gsl_matrix_memcpy(vtransi[10],pij);
         //               gsl_matrix_memcpy(vtransi[11],pij);
         //               gsl_matrix_memcpy(vtransi[13],pij);
         //            } else if (i==82){
         //               gsl_matrix_memcpy(vtransi[1],pij);
         //            } else if (i==107){
         //               gsl_matrix_memcpy(vtransi[9],pij);
         //            } else if (i==110){
         //               gsl_matrix_memcpy(vtransi[14],pij);
         //            } else if (i==148){
         //               gsl_matrix_memcpy(vtransi[12],pij);
         //            } else if (i==253){
         //               gsl_matrix_memcpy(vtransi[8],pij);
         //            }  		
         //            }
         // approx, integr_step=0.01
         for (unsigned int i=1;i<26;i++){
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
            pmattemp=pij;
            pij=pijp;
            pijp=pmattemp;
            //      gsl_matrix_memcpy(pij,pijp);
            if (i==1){
               gsl_matrix_memcpy(vtransi[15],pij);
               gsl_matrix_memcpy(vtransi[2],pij);
               gsl_matrix_memcpy(vtransi[3],pij);
               gsl_matrix_memcpy(vtransi[5],pij);
               gsl_matrix_memcpy(vtransi[7],pij);
            } else if (i==2){
               gsl_matrix_memcpy(vtransi[4],pij);
               gsl_matrix_memcpy(vtransi[16],pij);
            } else if (i==3){
               gsl_matrix_memcpy(vtransi[17],pij);
            } else if (i==6){
               gsl_matrix_memcpy(vtransi[6],pij);
            } else if (i==8){
               gsl_matrix_memcpy(vtransi[0],pij);
               gsl_matrix_memcpy(vtransi[10],pij);
               gsl_matrix_memcpy(vtransi[11],pij);
               gsl_matrix_memcpy(vtransi[13],pij);
               gsl_matrix_memcpy(vtransi[1],pij);
            } else if (i==11){
               gsl_matrix_memcpy(vtransi[9],pij);
               gsl_matrix_memcpy(vtransi[14],pij);
            } else if (i==15){
               gsl_matrix_memcpy(vtransi[12],pij);
            } else if (i==25){
               gsl_matrix_memcpy(vtransi[8],pij);
            }
         }
      }
   }
   //   for (unsigned int i=0;i<4;i++){
   //      for (unsigned int j=0;j<4;j++){
   //         cout << gsl_matrix_get(vtransi[8],i,j) << " ";
   //      }
   //      cout << endl;
   //   }
   //   cout << endl;

   double logli=0;
   for (ivma ima=mot.seqs.begin();ima!=mot.seqs.end();ima++){
      logli+=loglikely_column(pos,*ima,vtransi,w);
   }
   double w0=gsl_vector_get(w,0);
   double w1=gsl_vector_get(w,1);
   double w2=gsl_vector_get(w,2);
   double w3=1-w0-w1-w2;
   logli+=(alpha-1.)*(log(w0)+log(w1))+(beta-1.)*(log(w2)+log(w3));

   return logli;
}
