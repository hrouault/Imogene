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

   void
evolvesite(Motif & mot)
{

   //number of sites for evolution stat:
   unsigned int numforstat(100);

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
//   for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){
//      ivm->print();
//   }

   for (unsigned int n=1;n<nbspecies;n++){
      cout << numtospecies(n) << endl;
      for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){

         if (!ivm->matches[n]) continue;

         initscore=scoref(ivm->alignseq[0],mot.matprec);

         vint site(ivm->alignseq[0]);

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
// !!!! To be restored
//            Pdisp_h[i]=evolvedist_halpern(Initprob[i],mot.matfreq[i],distances[n]);
         }

         // FELSEN
         for (unsigned int i=0;i<mot.motwidth;i++){
// !!!! To be restored
//            Pdisp_f[i]=evolvedist_felsen(Initprob[i],mot.matfreq[i],distances[n]);
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
         outf << scoref(ivm->alignseq[n],mot.matprec)-initscore << endl;


         // draw numforstat sites for evolved sites scores computation
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

         }
         //         getchar();
      }
   }

   outf.close();
}

   void
test_args_init()
{
   if (!strcmp(test_args.species_arg, "droso")) {
      species = "droso";
      nbspecies = 12;
   } else if (!strcmp(test_args.species_arg, "eutherian")) {
      species = "eutherian";
      nbspecies = 12;
   }
   initconc();
   width = test_args.width_arg;
   scorethr2 = test_args.threshold_arg * log(2);
   scorethr = scorethr2 * (1 - 2.0 / width);
   scorethrcons = scorethr2 * (1 - 1.0 / width);
   evolutionary_model = test_args.evolutionary_model_arg;
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
   
   compalpha();

   Motif mot;
   mot = mots[0];
   mot.setscorethr2meaninfo();
   mot.matinit(mot.motscorethr2);

   evolvesite(mot);

   gsl_rng_free(gslran);
   test_cmdline_parser_free(&test_args);
   return EXIT_SUCCESS;
}		/* -----  end of function test  ----- */


