/*    
 * Copyright (C) 2003-2011 Hervé Rouault
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


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <string>
#include <cstring>

#include <algorithm> // used by sort

#include "genmot_cmdline.h"


// *** See if the following are required...
#include<cmath>
#include<iomanip>
#include<vector>
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

#include "config.h"

#include "const.hpp"
#include "vectortypes.hpp"
#include "random.hpp"
#include "sequence.hpp"
#include "motif.hpp"
#include "tree.hpp"
//#include "montecarlo.hpp"

#include "imogene.hpp"

#include "extract.hpp"
#include "genmot.hpp"
#include "help.hpp"
#include "scangen.hpp"

#include "distinfo.hpp"

genmot_args_info args_info;
vginst potregs;
vvginst groupedinst;
vstring phenos;
vstring gbacks;
vchrom chromints;

unsigned int sizepos,sizeneg;
unsigned int cutoff_for_combination=3;


const char usage_string[] =
	"imogene [--version] [--help]\n"
	"           <command> [<args>]";

const char more_info_string[] =
	"See 'imogene help <command>' for more information on a specific command.";


static int handle_options(char ***argv, int *argc, int *envchanged)
{
	char **orig_argv = *argv;

	while (*argc > 0) {
		char *cmd = (*argv)[0];
		if (cmd[0] != '-')
			break;

		if (!strcmp(cmd, "--help") || !strcmp(cmd, "--version")){
			break;
		} else {
			fprintf(stderr, "Unknown option: %s\n", cmd);
			fprintf(stderr, "Usage : %s\n", usage_string);
		}

		(*argv)++;
		(*argc)--;
	}
	return (*argv) - orig_argv;
}

struct cmd_struct {
	const char *cmd;
	int (*fn)(int, char **);
   const char *help;
};

static int run_builtin(struct cmd_struct *p, int argc, char **argv)
{
	int status, help;
//	struct stat st;
	const char *prefix;

	prefix = NULL;
	help = argc == 2 && !strcmp(argv[1], "-h");


	status = p->fn(argc, argv);
	if (status)
		return status;

	return 0;
}

static struct cmd_struct commands[] = {
   { "distinfo", cmd_distinfo, "Computes the distance between two motifs." },
   { "extract", cmd_extract, "Extracts alignments from coordinates." },
   { "genmot", cmd_genmot, "generate motifs" },
   { "help", cmd_help, "Help message" },
   { "scangen", cmd_scangen, "infere CRMs" },
   { "version", cmd_version, "Print Imogene version" }
};


static void
handle_command(int argc, char **argv)
{
	char *cmd = argv[0];
   for (int i = 0; i < ARRAY_SIZE(commands); i++) {
      struct cmd_struct *p = commands+i;
      if (strcmp(p->cmd, cmd))
         continue;
      exit(run_builtin(p, argc, argv));
   }
   cout << "error!" << endl;
}

static inline void mput_char(char c, unsigned int num)
{
	while(num--)
		putchar(c);
}


void list_cmds_help(void)
{
	int i, longest = 0;

	for (i = 0; i < ARRAY_SIZE(commands); i++) {
		if (longest < strlen(commands[i].cmd))
			longest = strlen(commands[i].cmd);
	}

	puts("The available Imogene commands are:");
	for (i = 0; i < ARRAY_SIZE(commands); i++) {
		printf("   %s   ", commands[i].cmd);
		mput_char(' ', longest - strlen(commands[i].cmd));
		puts(commands[i].help);
	}
}


   void
findnearestgene_sides(vcoord & vgenes, vcoord & vpeaks)
{

   double peak;
   for (ivcoord ivc=vpeaks.begin();ivc!=vpeaks.end();ivc++){
      cout << ivc->name << " -> ";
      peak=(ivc->stop+ivc->start)/2;
      double distmin(1e9);
      ivcoord bestcoord;
      for (ivcoord ivg=vgenes.begin();ivg!=vgenes.end();ivg++){
         if (ivg->chrom==ivc->chrom){
            double dist;
            dist=fabs(peak-ivg->start);
            if (fabs(peak-ivg->stop)<dist) dist=fabs(peak-ivg->stop);
            if (dist<distmin){
               ivc->name=ivg->name;
               bestcoord=ivg;
               distmin=dist;
            }
         }
      }

      if (distmin==1e6){
         ivc->neargenes.push_back("None");
         ivc->neargenes.push_back("None");
         ivc->neargenes.push_back("None");
         ivc->neargenes.push_back("None");
         continue;
      }

      //TSS in 3', or 5'
      if (peak-bestcoord->start>0){
         for (ivcoord ivg=bestcoord-1;ivg!=bestcoord+3;ivg++){
            if (ivg>=vgenes.begin() && ivg<vgenes.end() && ivg->chrom==bestcoord->chrom){
               ivc->neargenes.push_back(ivg->name);
            }
            else {
               ivc->neargenes.push_back("None");
            } 
         }
      }
      else{
         for (ivcoord ivg=bestcoord-2;ivg!=bestcoord+2;ivg++){
            if (ivg>=vgenes.begin() && ivg<vgenes.end() && ivg->chrom==bestcoord->chrom){
               ivc->neargenes.push_back(ivg->name);
            }
            else {
               ivc->neargenes.push_back("None");
            } 
         }
      }

      for (ivstring ivs=ivc->neargenes.begin();ivs!=ivc->neargenes.end();ivs++){
         cout << *ivs << "\t";
      }
      cout << endl;

   }
   return;
}

   void
findnearestgene(vcoord & vgenes, vcoord & vpeaks)
{

   //   ofstream outf("peak_dist.dat");
   double peak;
   for (ivcoord ivc=vpeaks.begin();ivc!=vpeaks.end();ivc++){
      //cout << ivc->name << " -> ";
      peak=(ivc->stop+ivc->start)/2;
      double distmin(1e12);
      double distminsigned(1e12);
      for (ivcoord ivg=vgenes.begin();ivg!=vgenes.end();ivg++){
         if (ivg->chrom==ivc->chrom){
            double dist;
            // distance to TSS, depends on strand
            //if (ivg->strand==-1) dist=peak-ivg->stop;
            //else dist=peak-ivg->start;
            dist=peak-ivg->start;
            //        if (fabs(peak-ivg->stop)<dist) dist=fabs(peak-ivg->stop);
            if (fabs(dist)<distmin){
               ivc->name=ivg->name;
               distmin=fabs(dist);
               if (ivg->strand==-1) distminsigned=-dist;
               else distminsigned=dist;
            }
         }
      }
      //      if (fabs(distminsigned)!=1e6) outf << distminsigned << " " << ivc->name << endl;
      //cout << ivc->name << endl;
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

      // *** CREATE A CPP FOR DISPLAY
      // *** OPTION FOR MASKING FOUND SITES?

//      if (args_info.disp_svg_given){
//         nm=(*im).nbmatchnmaskforsvg(seq,moti);
//      } else {
         //nm=(*im).nbmatchnmask(seq,moti);
         nm=(*im).nbmatchwomask(seq,moti);
//      }
      nbmot=nm;

      // *** ENSURE MOTIFS ARE WELL WEIGHTED (eg CONSEERVED INSTANCES ON BACKGROUND)
 //     if (args_info.weightmots_given){
         nmcorr+=nm*log((*im).lambdatrain/(*im).lambda);
  //    }
 //     else{
 //        nmcorr+=nm;
 //     }
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
print_reportbugs()
{
   cout << endl;
   cout << "Maintained by Hervé Rouault rouault@lps.ens.fr and" << endl;
   cout << "Marc Santolini santolin@lps.ens.fr" << endl;
   cout << "Report bugs to one of us." << endl;
}

   void
print_copyright ()
{

   cout << " Copyright (C) 2003-2011 Hervé Rouault" << endl;
   cout << endl;
   cout << " Imogene is free software: you can redistribute it and/or modify" << endl;
   cout << " it under the terms of the GNU General Public License as published by" << endl;
   cout << " the Free Software Foundation, either version 3 of the License, or" << endl;
   cout << " (at your option) any later version." << endl;
   cout << endl;
   cout << " Imogene is distributed in the hope that it will be useful," << endl;
   cout << " but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
   cout << " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
   cout << " GNU General Public License for more details." << endl;
   cout << endl;
   cout << " You should have received a copy of the GNU General Public License" << endl;
   cout << " along with Imogene.  If not, see <http://www.gnu.org/licenses/>." << endl;
   cout << endl;
   cout << " Written by Hervé Rouault and Marc Santolini." << endl;

}

#define is_dir_sep(c) ((c) == '/' || (c) == '\\')

char *
extract_argv0_cmd(char *argv0)
{
	char *slash;

	if (!argv0 || !*argv0)
		return NULL;
	slash = argv0 + strlen(argv0);

	while (argv0 <= slash && !is_dir_sep(*slash))
		slash--;

	if (slash >= argv0) {
		return slash + 1;
	}

	return argv0;
}

int prefixcmp(char *str, const char *prefix)
{
	for (; ; str++, prefix++)
		if (!*prefix)
			return 0;
		else if (*str != *prefix)
			return (unsigned char)*prefix - (unsigned char)*str;
}


/** 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  imogene-genmot main file, responsible for the de novo inference of
 *  motifs
 * =====================================================================================
 */
   int
main(int argc, char **argv)
{

   char *cmd;

	cmd = extract_argv0_cmd(argv[0]);
	if (!cmd)
		strncpy(cmd,"imogene-help",15);

   if (!prefixcmp(cmd, "imogene-")) {
		cmd += 8;
		argv[0] = cmd;
		handle_command(argc, argv);
      return EXIT_FAILURE;
	}
   

   argv++;
   argc--;
   handle_options(&argv, &argc, NULL);
   if (argc > 0){
      if (!strcmp(argv[0], "--help") || !strcmp(argv[0], "--version")){
         argv[0]+=2;
      }
   } else {
      cout <<  "Usage : " << usage_string << endl;
      list_common_cmds_help();
      cout << "\n" << more_info_string << endl;
      return EXIT_FAILURE;
   }
   cmd = argv[0];

   if (!prefixcmp(cmd, "imogene-")) {
		cmd += 8;
		argv[0] = cmd;
		handle_command(argc, argv);
      return EXIT_FAILURE;
	}

   cmd = argv[0];

   handle_command(argc, argv);

   fprintf(stderr, "Failed to run command '%s': %s\n",
		cmd, strerror(errno));

   return EXIT_FAILURE;
}

