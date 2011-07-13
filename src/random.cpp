#include <gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<iostream>
#include<fstream>

#include "vectortypes.hpp"
#include "const.hpp"
using namespace std;

gsl_rng * gslran;

   void
rnginit()
{
   const gsl_rng_type * T;
   gsl_rng_env_setup();
   T = gsl_rng_default;
   gslran = gsl_rng_alloc (T);
   long seed=time(NULL) * getpid();
   gsl_rng_set(gslran,seed);
}

void rngtest()
{
   vint ntimes;
   ntimes.push_back(10);
   ntimes.push_back(100);
   ntimes.push_back(1000);
   ntimes.push_back(10000);
   ntimes.push_back(100000);
   ntimes.push_back(1000000);
   ofstream outf("cv.dat");
   long seed=time(NULL) * getpid();
   // gsl_rng_set(gslran,seed);
   int nt=ntimes.size();
   vd moya(nt,0),moyt(nt,0),moyc(nt,0),moyg(nt,0);
   for (int imoy=0;imoy<50;imoy++){
      vvd diff;
      gsl_rng_set(gslran,seed);
      seed+=1;
      int ic=0;
      for (ivint in=ntimes.begin();in!=ntimes.end();in++){
         //   cout << "Time: " << time(NULL) << endl;
         //   cout << "Seed: " << seed << endl;
         double pa,pt,pc,pg,ptot;
         gsl_ran_discrete_t * g;
         unsigned int n[4];
         double p[4];
         p[0]=conca;
         p[1]=conct;
         p[2]=concc;
         p[3]=concg;
         ptot=0;
         g=gsl_ran_discrete_preproc (4,p);
         for (int j=0;j<*in;j++){
            vint bases;
            for (int i=0;i<10;i++) 
            { 
               unsigned int base=gsl_ran_discrete (gslran,g);
               if (base==0) pa+=1.;
               else if (base==1) pt+=1.;
               else if (base==2) pc+=1.;
               else if (base==3) pg+=1.;
               ptot+=1;
               bases.push_back(base);
            }
            //cout << bases << endl;
            //if (j<5) cout << bases << endl;
         }
         //cout << endl;
         pa/=ptot;
         pt/=ptot;
         pc/=ptot;
         pg/=ptot;

         moya[ic]+=abs(pa-gsl_ran_discrete_pdf(0,g))/50;
         moyt[ic]+=abs(pt-gsl_ran_discrete_pdf(1,g))/50;
         moyc[ic]+=abs(pc-gsl_ran_discrete_pdf(2,g))/50;
         moyg[ic]+=abs(pg-gsl_ran_discrete_pdf(3,g))/50;
         //cout << "pa=" << pa << " pt=" << pt << " pp=" << pc << " pg=" << pg << endl;
         //cout << "pa=" << gsl_ran_discrete_pdf(0,g) << " pt=" << gsl_ran_discrete_pdf(1,g) << " pp=" << gsl_ran_discrete_pdf(2,g) << " pg=" << gsl_ran_discrete_pdf(3,g) << endl;

         gsl_ran_discrete_free (g);
         ic++;
      }
   }

   for (int im=0;im!=nt;im++){
      outf << moya[im] << " " << moyt[im] << " " << moyc[im] << " " << moyg[im] << endl;
   }
   outf.close();
}


