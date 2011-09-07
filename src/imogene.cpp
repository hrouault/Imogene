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

#include <algorithm> // used by sort

#include "genmot_cmdline.h"


// *** See if the following are required...
#include<cmath>
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
vmot motsdef;
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

//	/* Somebody closed stdout? */
//	if (fstat(fileno(stdout), &st))
//		return 0;
//	/* Ignore write errors for pipes and sockets.. */
//	if (S_ISFIFO(st.st_mode) || S_ISSOCK(st.st_mode))
//		return 0;

//	/* Check for ENOSPC and EIO errors.. */
//	if (fflush(stdout))
//		die_errno("write failure on standard output");
//	if (ferror(stdout))
//		die("unknown write failure on standard output");
//	if (fclose(stdout))
//		die_errno("close failed on standard output");
	return 0;
}

static struct cmd_struct commands[] = {
   { "distinfo", cmd_distinfo, "Computes the distance between two motifs." },
   { "extract", cmd_extract, "Extracts alignments from coordinates." },
   { "genmot", cmd_genmot, "generate motifs" },
   { "help", cmd_help, "Help message" },
   { "scangen", cmd_scangen, "infere CRMs" }
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

//   void
//printbackreg(vcoord & vcds,string folder)
//{
//   //ofstream outf;
//   ofstream outf("randpeaks.dat");
//   unsigned int numback(0);
//   unsigned int i(0);
//   unsigned int backsize;
//   if (args_info.size_given) backsize=args_info.size_arg;
//   else backsize=2000;
//   while (numback<10000){
//      if (i>=vcds.size()) i=0;
//      int range;
//      range=vcds[i].stop-backsize-vcds[i].start;
//      //cout << vcds[i];
//      if (range<0){
//         i++;
//         //     cout << "range < 0" << endl;
//         continue;
//      }
//      unsigned int newstart;
//      unsigned int newpos;
//      newpos= gsl_rng_uniform_int (gslran,range+1);
//      newstart=vcds[i].start+newpos;
//      Coordinate cdtmp;
//      cdtmp.start=newstart;
//      cdtmp.stop=newstart+backsize-1;
//      cdtmp.chrom=vcds[i].chrom;
//      cdtmp.name=vcds[i].name;
//      // FOR PEAKS
//      // 
//      outf << "peaks_" << numback+1 << "\t";
//      outf << chromfromint(cdtmp.chrom) << "\t";
//      outf << cdtmp.start << "\n";
//      numback++;
//
//      //FOR SEQS
//      //
//      //         Sequence & s=seq;
//      //         for (int k=0;k<nbspecies;k++){
//      //            if (s.species[k]){
//      //               if (k==0) outf << ">" << numtospecies(k) << " " <<
//      //                  "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
//      //               else  outf << ">" << numtospecies(k) << endl;
//      //               outf << s.seqsrealigned[k] << endl;
//      //            }; 
//      //         }
//      //         outf.close();
//      //         numback++;
//      //      }
//      //      Sequence seq=coordtoseq(cdtmp);
//      //      int nbfr=0;
//      //      bool bc=0;
//      //      if (species==1){
//      //         if (seq.species[5]) nbfr++; 
//      //         if (seq.species[6] || seq.species[7]) nbfr++;
//      //         if (seq.species[8]) nbfr++;
//      //         if (seq.species[9] || seq.species[10] || seq.species[11]) nbfr++;
//      //         if (nbfr>1) bc=1;
//      //      }
//      //      else if (species==2){
//      //         if (seq.species[2] || seq.species[3] || seq.species[4] || seq.species[5]) nbfr++;
//      //         if (seq.species[6] || seq.species[7]) nbfr++;
//      //         if (seq.species[8]) nbfr++;
//      //         if (seq.species[9]) nbfr++;
//      //         if (nbfr>1) bc=1;
//      //      }
//      //      if (seq.species[0] && bc && seq.nbtb>backsize/2){
//      //         cout << "->" << numback << ". " <<cdtmp;
//      //         stringstream os;
//      //         os << numback;
//      //         os >> seq.name;
//      //         string filename=folder;
//      //         filename+=seq.name;
//      //         filename.append(".fa");
//      //         outf.open(filename.c_str());
//      //         Sequence & s=seq;
//      //         for (int k=0;k<nbspecies;k++){
//      //            if (s.species[k]){
//      //               if (k==0) outf << ">" << numtospecies(k) << " " <<
//      //                  "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
//      //               else  outf << ">" << numtospecies(k) << endl;
//      //               outf << s.seqsrealigned[k] << endl;
//      //            }; 
//      //         }
//      //         outf.close();
//      //         numback++;
//      //      }
//      //      else {
//      //         //cout << "sequence not conserved or too much masked" << endl;
//      //      }
//      cout.flush();
//      i++;
//}
//outf.close();
//}

//void
//printbackregwcoords(vcoord & vcds,string folder){
//
//   //vcds is a vector of random shuffled intergenic regions
//   //regints are our interest sequences
//   //
//   ofstream outf;
//   unsigned int i(0);
//   unsigned int numback(1);
//
//   unsigned int repeat=10;
//   if (args_info.numrepeat_given) repeat=args_info.numrepeat_arg;
//
//   for (int j=1;j<=repeat;j++){
//      for (ivseq ivs=regints.begin();ivs!=regints.end();ivs++){
//
//         //size of our interest sequence
//         unsigned int backsize=ivs->stop-ivs->start+1;
//         unsigned int state=0;
//
//         //if we find a corresponding background region, then state=1;
//         while (state==0){
//
//            //if we have searched all intergenic sequences we loop on the first
//            if (i==vcds.size()-1) i=0;
//
//            int range;
//            range=vcds[i].stop-backsize-vcds[i].start+1;
//
//            //intergenic size has to be superior to interest sequence
//            if (range<1){
//               i++;
//               continue;
//            }
//
//            //randomly choose a region of a same size than interest in the intergenic region
//            unsigned int newstart;
//            unsigned int newpos;
//            newpos= gsl_rng_uniform_int (gslran,range);
//            newstart=vcds[i].start+newpos;
//            Coordinate cdtmp;
//            cdtmp.start=newstart;
//            cdtmp.stop=newstart+backsize-1;
//            cdtmp.chrom=vcds[i].chrom;
//            cdtmp.name=vcds[i].name;
//
//            //extract the alignment
//            Sequence seq=coordtoseq(cdtmp);
//            int nbfr=0;
//            bool bc=0;
//
//            if (species=="droso"){
//               if (seq.species[5]) nbfr++; 
//               if (seq.species[6] || seq.species[7]) nbfr++;
//               if (seq.species[8]) nbfr++;
//               if (seq.species[9] || seq.species[10] || seq.species[11]) nbfr++;
//               if (nbfr>1) bc=1;
//            }
//            else if (species=="eutherian"){
//               if (seq.species[2] || seq.species[3] || seq.species[4] || seq.species[5]) nbfr++;
//               if (seq.species[6] || seq.species[7]) nbfr++;
//               if (seq.species[8]) nbfr++;
//               if (seq.species[9]) nbfr++;
//               if (nbfr>1) bc=1;
//            }
//
//
//            //we want mus to be present, the sequence to be conserved, and at least half of bases to be unmasked
//            if (seq.species[0] && bc && seq.nbtb>backsize/2){
//               cout << "->" << numback << ". " <<cdtmp;
//               cout.flush();
//               stringstream os;
//               os << numback;
//               os >> seq.name;
//               string filename=folder;
//               filename+=seq.name;
//               filename.append(".fa");
//               outf.open(filename.c_str());
//               Sequence & s=seq;
//               for (int k=0;k<nbspecies;k++){
//                  if (s.species[k]){
//                     if (k==0) outf << ">" << numtospecies(k) << " " <<
//                        "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
//                     else  outf << ">" << numtospecies(k) << endl;
//                     outf << s.seqsrealigned[k] << endl;
//                  }; 
//               }
//               outf.close();
//               numback++;
//               state=1;
//            }
//            i++;
//         }
//
//      }
//   }
//
//   outf.close();
//}

   void
disptexinit(ofstream & outf)
{
   outf << "\\documentclass[11pt,twoside,reqno,a4paper]{article}\n"<<
      "\\usepackage[french]{babel}\n"<<
      "\\usepackage{color}\n"<<
      "\\usepackage[utf8]{inputenc}\n"<<
      "\\usepackage{geometry}\n"<<
      "\\geometry{a4paper}\n"<<          
      "\\usepackage{graphicx}\n"<<
      "\\usepackage{array,multirow}\n" <<
      "\\begin{document}\n"<<
      "\\noindent\n";
}

   void
dispinit(ofstream & outf)
{
   outf << "\\documentclass[11pt,twoside,reqno,a4paper]{article}\n"<<
      "\\usepackage[french]{babel}\n"<<
      "\\usepackage{color}\n"<<
      "\\usepackage[utf8]{inputenc}\n"<<
      "\\usepackage{geometry}\n"<<
      "\\geometry{a4paper}\n"<<          
      "\\usepackage{graphicx}\n"<<
      "\\usepackage{array,multirow,longtable}\n" <<
      "\\setlongtables\n" <<
      "\\begin{document}\n"<<
      "\\noindent\n"<<
      "\\begin{longtable}[!t]\n" <<
      "{|c|c|c|c|}\n" <<
      "\\hline\n" <<
      "Nom & \\# cons & Score & Pos/Neg " <<
      "\\endfirsthead\n" <<
      "\\hline\n" <<
      "Nom & \\# cons & Score & Pos/Neg\\\\\n" <<
      "\\hline\n" <<
      "\\endhead\n" <<
      "\\endlastfoot \\hline\n";
}

   void
dispinit(ofstream & outf,vmot & vmoti)
{
   outf << "\\documentclass[11pt,twoside,reqno,a4paper]{article}\n"<<
      "\\usepackage[french]{babel}\n"<<
      "\\usepackage{color}\n"<<
      "\\usepackage[utf8]{inputenc}\n"<<
      "\\usepackage{geometry}\n"<<
      "\\geometry{a4paper}\n"<<          
      "\\usepackage{graphicx}\n"<<
      "\\usepackage{array,multirow,longtable}\n" <<
      "\\setlongtables\n" <<
      "\\begin{document}\n"<<
      "\\noindent\n"<<
      "\\begin{longtable}[!t]\n" <<
      "{|c|";
   for (int i=0;i<vmoti.size();i++){
      outf << "c|";
   }
   outf <<"c|c|}\n" <<
      "\\hline\n" <<
      "Nom &";
   for (ivmot ivm=vmoti.begin();ivm!=vmoti.end();ivm++){
      outf << "\\#" << ivm->index+1 << " &";
   }
   outf <<" Score & Pos/Neg " <<
      "\\endfirsthead\n" <<
      "\\hline\n" <<
      "Nom &";
   for (ivmot ivm=vmoti.begin();ivm!=vmoti.end();ivm++){
      outf << "\\#" << ivm->index+1 << " &";
   }
   outf <<" Score & Pos/Neg " <<
      "\\hline\n" <<
      "\\endhead\n" <<
      "\\endlastfoot \\hline\n";
}

   void
dispinitforrank(ofstream & outf,vmot & vmoti)
{
   outf << "\\documentclass[11pt,landscape,twoside,reqno,a4paper]{article}\n"<<
      "\\usepackage[french]{babel}\n"<<
      "\\usepackage{color}\n"<<
      "\\usepackage[utf8]{inputenc}\n"<<
      "\\usepackage{geometry}\n"<<
      "\\geometry{a4paper}\n"<<          
      "\\usepackage{graphicx}\n"<<
      "\\usepackage{array,multirow,longtable}\n" <<
      "\\setlongtables\n" <<
      "\\begin{document}\n"<<
      "\\noindent\n"<<
      "\\begin{longtable}[!t]\n" <<
      "{|c|c|c|c|";
   for (int i=0;i<vmoti.size();i++){
      outf << "c|";
   }
   outf <<"c|}\n" <<
      "\\hline\n" <<
      "Nom & Chrom & Start & Stop &";
   for (ivmot ivm=vmoti.begin();ivm!=vmoti.end();ivm++){
      if (ivm->name!=""){
         outf << ivm->name << " &";
      }
      else{
         outf << "\\#" << ivm->index+1 << " &";
      }
   }
   outf <<" Score" <<
      "\\endfirsthead\n" <<
      "\\hline\n" <<
      "Nom & Chrom & Start & Stop &";
   for (ivmot ivm=vmoti.begin();ivm!=vmoti.end();ivm++){
      if (ivm->name!=""){
         outf <<  ivm->name << " &";
      }
      else{
         outf << "\\#" << ivm->index+1 << " &";
      }
   }
   outf <<" Score" <<
      "\\endhead\n" <<
      "\\endlastfoot \\hline\n";
}

   void
dispclose(ofstream & outf)
{
   outf <<	"\\end{document}";
}

   void
dispscore(vseq & vscore, Motif & mot)
{
   ofstream disp;
   ostringstream disps;
   disps << "score/dispscore/dispscore_" << mot.index << ".tex";
   disp.open(disps.str().c_str());
   dispinit(disp);
   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
      string sign;
      if ((*ivs).sign==1) sign="+";
      else sign="-";
      disp << (*ivs).name << " & " << (*ivs).motis[mot.index] << " & "
         << setprecision(3) <<  (*ivs).score << " & " << sign << " \\\\\n" <<
         "\\hline\n";

   }
   disp  << "\\caption{\\small Motif " << mot.index << 
      " at threshold " << scorethr2 << "}\n"<<
      "\\end{longtable}\n";
   dispclose(disp);
   disp.close();
}

   void
dispscore(vseq & vscore, vmot & vmoti)
{
   ofstream disp;
   ostringstream disps;
   disps << "score/dispscore/dispscore_";
   for (ivmot mot=vmoti.begin();mot!=vmoti.end();mot++){
      if (mot!=vmoti.begin()) disps << "+";
      disps << mot->index;
   }
   disps << ".tex";
   disp.open(disps.str().c_str());
   dispinit(disp,vmoti);
   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
      string sign;
      if ((*ivs).sign==1) sign="+";
      else sign="-";
      disp << (*ivs).name << " & "; 
      for (ivmot ivm=vmoti.begin();ivm!=vmoti.end();ivm++){
         disp << (*ivs).motis[ivm->index] << " & ";
      }
      disp << setprecision(3) <<  (*ivs).score << " & " << sign << " \\\\\n" <<
         "\\hline\n";

   }
   disp  << "\\caption{\\small Motifs ";

   for (ivmot ivm=vmoti.begin();ivm!=vmoti.end();ivm++){
      disp << ivm->index << " ";
   }
   disp <<" at threshold " << scorethr2 << "}\n"<<
      "\\end{longtable}\n";
   dispclose(disp);
   disp.close();
}

   void
dispscoreforrank(vseq & vscore, vmot & vmoti,ofstream & disp)
{
   dispinitforrank(disp,vmoti);
   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
      disp << (*ivs).name << " & "; 
      disp << chromfromint(ivs->chrom) << " & "; 
      disp << ivs->start << " & "; 
      disp << ivs->stop << " & "; 
      for (ivmot ivm=vmoti.begin();ivm!=vmoti.end();ivm++){
         disp << (*ivs).motis[ivm->index] << " & ";
      }
      disp << setprecision(3) <<  (*ivs).score << " \\\\\n" <<
         "\\hline\n";

   }
   disp  << "\\caption{\\small Motifs ";

   for (ivmot ivm=vmoti.begin();ivm!=vmoti.end();ivm++){
      if (ivm->name!=""){
         disp <<  ivm->name << " ";
      }
      else{
         disp << "\\#" << ivm->index+1 << " &";
      }
   }
   disp <<" at threshold " << scorethr2 << "}\n"<<
      "\\end{longtable}\n";
   dispclose(disp);
}

   void
dispseqwmots (Sequence & seq, vmot & mots, ofstream & outf,double scorethrinit)
{
   for (unsigned int i=0;i<1;i++){//nbspecies;i++){}
      unsigned int texpos=0;
      if(seq.species[i]){
         unsigned int pos=0;
         string name;
         name=numtospecies(i);
         int found=seq.name.find("_");
         string texname=seq.name;
         if(found!=string::npos) texname.insert(found,"\\");
         if (i==0){
            outf << "$>$" << name <<  "\t" << texname << "\tchr" << chromfromint(seq.chrom) << 
               "\t" << seq.start << "\t" << seq.stop << " \\\\\n";
         } else{
            outf << "$>$" << name <<  "\\\\\n";
         }
         outf << "\\texttt{";

         outf << pos+seq.start << "\t";

         double widthinit=10;

         int minwidth=mots[0].bsinit.size();
         if (mots.size()>1 && mots[1].bsinit.size()<minwidth) minwidth=mots[1].bsinit.size();
         if (mots.size()>2 && mots[2].bsinit.size()<minwidth) minwidth=mots[2].bsinit.size();
         if (mots.size()>3 && mots[3].bsinit.size()<minwidth) minwidth=mots[3].bsinit.size();


         civint istrf;
         for (civint istr=seq.iseqs[i].begin();istr!=seq.iseqs[i].end()-minwidth+1;istr++){
            double score=0;
            width=mots[0].bsinit.size();
            if (istr>seq.iseqs[i].end()-width) break;
            int state=mots[0].statemot(seq,pos,i,score);
            //if (score>100) exit(9);
            string color="red";
            if (state==0 && mots.size()>1){
               width=mots[1].bsinit.size();
               if (istr>seq.iseqs[i].end()-width) break;
               state=mots[1].statemot(seq,pos,i,score);
               color="blue";
            }
            if (state==0 && mots.size()>2){
               width=mots[2].bsinit.size();
               if (istr>seq.iseqs[i].end()-width) break;
               state=mots[2].statemot(seq,pos,i,score);
               color="green";
            }
            if (state==0 && mots.size()>3){
               width=mots[3].bsinit.size();
               if (istr>seq.iseqs[i].end()-width) break;
               state=mots[3].statemot(seq,pos,i,score);
               color="yellow";
            }
            if (state==0){
               vint oneb;
               oneb.push_back(*istr);
               outf << vinttostring(oneb);
               pos++;
               texpos++;
            }
            else if (state==1){
               outf << "\\textcolor{" << color << "}{";
               for (unsigned int j=0;j<width;j++){
                  vint oneb;
                  oneb.push_back(*(istr+j));
                  outf << vinttostring(oneb);
                  pos++;
                  texpos++;
                  //                  if (texpos>0 && texpos%10==0 && texpos%60!=0){
                  //                     outf << "\t";
                  //                  }
                  if (texpos>0 && texpos%60==0){
                     outf << "}";
                     outf << "\\\\\n" << pos+seq.start << "\t";
                     outf << "\\textcolor{" << color << "}{";
                  }
               }
               for (int ipos=0;ipos<6;ipos++){
                  texpos++;
                  if (texpos>0 && (texpos)%60==0){
                     outf << "}";
                     outf << "\\\\\n" << pos+seq.start << "\t";
                     outf << "\\textcolor{" << color << "}{";
                  }
               }
               outf << " ($" << setprecision(2) << score << "$)";
               outf << "} ";
               istr+=width-1;
            }
            else if (state==2){
               outf << "\\textit{\\textcolor{" << color << "}{";
               for (unsigned int j=0;j<width;j++){
                  vint oneb;
                  oneb.push_back(*(istr+j));
                  outf << vinttostring(oneb);
                  pos++;
                  texpos++;
                  //                  if (texpos>0 && texpos%10==0 && texpos%60!=0){
                  //                     outf << "\t";
                  //                  }
                  if (texpos>0 && texpos%60==0){
                     outf << "}}";
                     outf << "\\\\\n" << pos+seq.start << "\t";
                     outf << "\\textit{\\textcolor{" << color << "}{";
                  }
               }
               outf << "}} ";
               for (int ipos=0;ipos<6;ipos++){
                  texpos+=1;
                  if (texpos>0 && (texpos)%60==0){
                     outf << "\\\\\n" << pos+seq.start << "\t";
                  }
               }
               outf << "\\textcolor{" << color << "}{";
               outf << " ($" << setprecision(2) << score << "$)";
               outf << "} ";
               istr+=width-1;
            }
            //            if (texpos>0 && texpos%10==0 && texpos%60!=0){
            //               outf << "\t";
            //            }
            if (texpos>0 && texpos%60==0){
               outf << "\\\\\n" << pos+seq.start << "\t";
            }

            istrf=istr;
         }
         if (istrf<seq.iseqs[i].end()-1){
            for (civint istr=istrf;istr!=seq.iseqs[i].end();istr++){
               vint oneb;
               oneb.push_back(*istr);
               outf << vinttostring(oneb);
               pos++;
               texpos++;
               //            if (texpos>0 && texpos%10==0 && texpos%60!=0){
               //               outf << "\t";
               //            }
               if (texpos>0 && texpos%60==0){
                  outf << "\\\\\n" << pos+seq.start << "\t";
               }
            }
         }
         outf << "}\n";
         outf << "\\\\\n";
      }
      //      i++;
      //      cout << "ok" << endl;
   }
   outf << "\\\\\n";
}

   string
colfromint(int i)
{
   if (i==0) return "red";
   if (i==1) return "blue";
   if (i==2) return "green";
   if (i==3) return "yellow";

   return "black";

}

   void
dispseqwmotswgaps (Sequence & seq, ofstream & outf)
{
   //HEADER
   string name=numtospecies(0);
   int found=seq.name.find("_");
   string texname=seq.name;
   if(found!=string::npos) texname.insert(found,"\\");
   outf << "$>$" << name <<  "\t";
   outf << texname << "\t" ;
   outf << "chr" << chromfromint(seq.chrom) << "\t";
   outf << seq.start << "\t" ;
   outf << seq.stop << " \\\\\n";
   outf << " \\\\\n";

   //DEFINING STATES, COLORS AND SCORES   
   unsigned int sizeseq=seq.seqsrealigned[0].size();
   vint vdum(sizeseq,0);
   vvint vvstate(nbspecies,vdum);// 0 nothing, 1 normal tfbs, 2 cons tfbs
   vvint vvcol(nbspecies,vdum); // 0 red, 1 blue, 2 green, 3 yellow
   vd vs(sizeseq,0.);
   vvd scores(nbspecies,vs);

   //NOT CONS
   for (ivinstseq ivs=seq.instances.begin();ivs!=seq.instances.end();ivs++){
      int spe=ivs->species;
      int motwidth=motsdef[ivs->motindex].motwidth;
      int ipos=seq.imaps[spe][ivs->pos];
      int stop=seq.imapsinv[spe][ipos+motwidth-1];

      for (int i=ivs->pos;i<stop+1;i++){
         vvstate[spe][i]=1;
         vvcol[spe][i]=ivs->motindex;
         scores[spe][i]=ivs->score;
         //if (seq.seqsrealigned[spe][i]=='-') vvcol[spe][i]=-1;
      }
   }

   //CONS
   for (ivvinstseq ivvs=seq.instancescons.begin();ivvs!=seq.instancescons.end();ivvs++){
      for (ivinstseq ivs=(*ivvs).begin();ivs!=(*ivvs).end();ivs++){
         int spe=ivs->species;
         int motwidth=motsdef[ivs->motindex].motwidth;
         int ipos=seq.imaps[spe][ivs->pos];
         int stop=seq.imapsinv[spe][ipos+motwidth-1];

         for (int i=ivs->pos;i<stop+1;i++){
            vvstate[spe][i]=2;
            scores[spe][i]=ivs->score;
         }
      }
   }

   //DISPLAY
   int start=0;
   int stop=min(60,sizeseq);
   outf << "\\texttt{";
   while (start<sizeseq){
      for (unsigned int spe=0;spe<nbspecies;spe++){
         if(seq.species[spe]){
            outf << numtospecies(spe) << "\t";
            for (unsigned int i=start;i<stop;i++){
               if (vvstate[spe][i]==0){
                  outf << seq.seqsrealigned[spe][i];
               } else if (vvstate[spe][i]==1){
                  outf << "\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                  outf << seq.seqsrealigned[spe][i];
                  outf << "}";

               } else if (vvstate[spe][i]==2){
                  outf << "\\textit{\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                  outf << seq.seqsrealigned[spe][i];
                  outf << "}}";

               }
            }
            outf << "\\\\\n";
            outf << "\\textcolor{white}{";
            outf << numtospecies(spe) << "\t";
            outf << "}";
            int pstate=0;
            for (unsigned int i=start;i<stop;i++){
               if (vvstate[spe][i]>0  && vvstate[spe][i]!=pstate){
                  if (vvstate[spe][i]==2) outf << "\\textit{";
                  outf << "\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                  outf << setprecision(1) << fixed <<  scores[spe][i] ;
                  outf << "}";
                  if (vvstate[spe][i]==2) outf << "}";
                  i+=2;
               } else{
                  outf << "\\ ";
               }
               pstate=vvstate[spe][i];
            }
            outf << "\\\\\n";
         }
      }
      outf << "\\\\\n";
      start=stop;
      stop=min(stop+60,sizeseq);
   }

   outf << "}\n";
}

   void
disptexclose(ofstream & outf)
{
   outf << "\\end{document}";
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

//   void
//scanseqs(vstring & regs)
//{
//   vcoord pcoords;
//   if (args_info.masktrain_given){
//
//      cout << "(Masking training set)" << endl;
//      ifstream inf;
//      inf.open(args_info.masktrain_arg);
//      back_insert_iterator<vcoord> pdest(pcoords);
//      copy(iiscoord(inf),iiscoord(),pdest);
//      inf.close();
//
//   }
//   
//   string pchrom("");
//   for (ivstring is=regs.begin();is!=regs.end();is++){
//      //cout << *is << endl;
//      Sequence seq;
//      ifstream fseq;
//      fseq.open((*is).c_str());
//      string fseqline;
//      getline(fseq,fseqline);
//      seq.name=fseqline;
//      stringstream firstline(fseqline);
//      string speciesname;
//      firstline >> speciesname;
//      string chrom;
//      firstline >> chrom;
//      if (chrom!=pchrom){
//         cout << chrom << endl;
//         pchrom=chrom;
//      }
//      seq.chrom=intfromchrom(chrom.substr(3));
//      if (seq.chrom==-1) continue;
//      //cout << "chromosome \"" << seq.chrom << "\"\n";
//      firstline >> seq.start;
//      //cout << "start \"" << seq.start << "\"\n";
//      firstline >> seq.stop;
//      seq.finame=*is;
//      string dum;
//      seq.species=vint(nbspecies,0);
//      vint dumi;
//      seq.iseqs=vvint(nbspecies,dumi);
//      seq.imaps=vvint(nbspecies,dumi);
//      seq.imapsinv=vvint(nbspecies,dumi);
//      while (!fseq.eof()){
//         //                cout << fseqline << endl;
//         int seqnum=speciestonum(fseqline.substr(1,6));
//         getline(fseq,fseqline);
//         seq.imaps[seqnum]=alignedtomap(fseqline);
//         seq.imapsinv[seqnum]=alignedtorevmap(fseqline);
//         string seqwogap=remgaps(fseqline);
//         if (seqwogap.size()>30){
//            seq.species[seqnum]=1;
//         }
//         seq.iseqs[seqnum]=stringtoint(seqwogap);
//         getline(fseq,fseqline);
//      }
//      seq.nbN=compN(seq.iseqs[0]);
//      seq.nbtb=seq.iseqs[0].size()-seq.nbN;
//      
//      for (ivcoord ivc=pcoords.begin();ivc!=pcoords.end();ivc++){
//         if (ivc->chrom==seq.chrom){ 
//            if ((ivc->start>=seq.start && ivc->start<=seq.stop) ||
//                  (ivc->stop>=seq.start && ivc->stop <=seq.stop)){
//                  int maskstart=max(ivc->start,seq.start);
//                  int maskstop=min(ivc->stop,seq.stop);
//      //            cout << ivc->name << endl;
//      //            cout << vinttostring(seq.iseqs[0]) << endl;
//                  for (unsigned int base=maskstart-seq.start;base<=maskstop-seq.start;base++){
//                     seq.iseqs[0][base]=4;
//                  }
//      //            cout << vinttostring(seq.iseqs[0]) << endl;
//            }
//         }
//      }
//
//      //check that ref species contains at least 30bp
//      if (seq.species[0]){
//         scanseq(seq,motsdef);
//      }
//
//      fseq.close();
//   }
//   unsigned int i=0;
//
//   cout << "Finding Instances" << endl;
//   for (ivmot im=motsdef.begin();im!=motsdef.begin()+nbmots_for_score;im++){
//      for (unsigned int j=0;j<nbchrom;j++){
//         //           cout << "Instances size : " << (*im).instances.size() << "\n";
//         sort((*im).instances[j].begin(),(*im).instances[j].end());
//         for (ivinst iins=(*im).instances[j].begin();iins!=(*im).instances[j].end();iins++){
//            Instance & curinst=*iins;
//            int start=0;
//            if ((curinst.coord-scanwidth)/scanstep+2>0) start=(curinst.coord-scanwidth)/scanstep+1;
//            for (unsigned int k=start;k<curinst.coord/scanstep+1;k++){
//               groupedinst[curinst.chrom][k].nbmots[i]++;
//               groupedinst[curinst.chrom][k].totmots++;
//               groupedinst[curinst.chrom][k].instances.push_back(curinst);
//            }
//         }
//      }
//      i++;
//   }
//
//   for (unsigned int j=0;j<nbchrom;j++){
//      for (ivginst ivg=groupedinst[j].begin();ivg!=groupedinst[j].end();ivg++){
//         if ((*ivg).totmots==0){
//            (*ivg).discarded=1;
//         }
//      }
//   }
//   for (unsigned int j=0;j<nbchrom;j++){
//      for (ivginst ivg=groupedinst[j].begin();ivg!=groupedinst[j].end();ivg++){
//         if (!(*ivg).discarded){
//            potregs.push_back(*ivg);
//         }
//      }
//   }
//
//   cout << "output results" << endl;
//   if (args_info.all_mots_only_given){
//         outputresults(nbmots_for_score);
//   }
//   else{
//      for (unsigned int i=1;i<nbmots_for_score+1;i++){
//         cout << "Motif " << i << endl;
//         outputresults(i);
//      }
//   }
//
//}

//   void
//scanseqs(ifstream & list,vcoord & coords)
//{
//   
//   vstring regs;
//   back_insert_iterator<vstring> dest(regs);
//   copy(iisstring(list),iisstring(),dest);
//
//   string pchrom("");
//   for (ivstring is=regs.begin();is!=regs.end();is++){
//      Sequence seq;
//      ifstream fseq;
//      fseq.open((*is).c_str());
//      string fseqline;
//      getline(fseq,fseqline);
//      seq.name=fseqline;
//      stringstream firstline(fseqline);
//      string speciesname;
//      firstline >> speciesname;
//      string chrom;
//      firstline >> chrom;
//      if (chrom!=pchrom){
//         cout << chrom << endl;
//         pchrom=chrom;
//      }
//      seq.chrom=intfromchrom(chrom.substr(3));
//      if (seq.chrom==-1) continue;
//      firstline >> seq.start;
//      firstline >> seq.stop;
//      seq.finame=*is;
//      string dum;
//      seq.species=vint(nbspecies,0);
//      vint dumi;
//      seq.iseqs=vvint(nbspecies,dumi);
//      seq.imaps=vvint(nbspecies,dumi);
//      seq.imapsinv=vvint(nbspecies,dumi);
//      while (!fseq.eof()){
//         //                cout << fseqline << endl;
//         int seqnum=speciestonum(fseqline.substr(1,6));
//         getline(fseq,fseqline);
//         seq.imaps[seqnum]=alignedtomap(fseqline);
//         seq.imapsinv[seqnum]=alignedtorevmap(fseqline);
//         string seqwogap=remgaps(fseqline);
//         if (seqwogap.size()>30){
//            seq.species[seqnum]=1;
//         }
//         seq.iseqs[seqnum]=stringtoint(seqwogap);
//         getline(fseq,fseqline);
//      }
//      seq.nbN=compN(seq.iseqs[0]);
//      seq.nbtb=seq.iseqs[0].size()-seq.nbN;
//      
//      scanseq(seq,motsdef);
//
//      fseq.close();
//   }
//   unsigned int i=0;
//
//   cout << "Finding Instances" << endl;
//   for (ivmot im=motsdef.begin();im!=motsdef.begin()+nbmots_for_score;im++){
//      for (unsigned int j=0;j<nbchrom;j++){
//         //           cout << "Instances size : " << (*im).instances.size() << "\n";
//         sort((*im).instances[j].begin(),(*im).instances[j].end());
//         for (ivinst iins=(*im).instances[j].begin();iins!=(*im).instances[j].end();iins++){
//            Instance & curinst=*iins;
//            int start=0;
//            if ((curinst.coord-scanwidth)/scanstep+2>0) start=(curinst.coord-scanwidth)/scanstep+1;
//            for (unsigned int k=start;k<curinst.coord/scanstep+1;k++){
//               groupedinst[curinst.chrom][k].nbmots[i]++;
//               groupedinst[curinst.chrom][k].totmots++;
//               groupedinst[curinst.chrom][k].instances.push_back(curinst);
//            }
//         }
//      }
//      i++;
//   }
//
//   cout << "Discarding empty Instances" << endl;
//   for (unsigned int j=0;j<nbchrom;j++){
//      for (ivginst ivg=groupedinst[j].begin();ivg!=groupedinst[j].end();ivg++){
//         if ((*ivg).totmots==0){
//            (*ivg).discarded=1;
//         }
//      }
//   }
//   for (unsigned int j=0;j<nbchrom;j++){
//      for (ivginst ivg=groupedinst[j].begin();ivg!=groupedinst[j].end();ivg++){
//         if (!(*ivg).discarded){
//            potregs.push_back(*ivg);
//         }
//      }
//   }
//
//   cout << "output results" << endl;
//   if (args_info.all_mots_only_given){
//         outputresults(nbmots_for_score);
//   }
//   else{
//      for (unsigned int i=1;i<nbmots_for_score+1;i++){
//         cout << "Motif " << i << endl;
//         outputresults(i);
//      }
//   }
//}
//   
//   void
//scanmots()
//{
//   ifstream potregs;
//   if (species=="droso") potregs.open("/home/santolin/these/files/droso/align/all/align-files.dat");
//   else if (species=="eutherian") potregs.open("/home/santolin/these/files/mus/epo/align-files.dat");
//   //else if (species==2) potregs.open("/home/santolin/these/files/transfac/matrices/align-test.dat");
//   vstring regs;
//   back_insert_iterator<vstring> dest(regs);
//   copy(iisstring(potregs),iisstring(),dest);
//   potregs.close();
//
//   string pchrom("");
//   system("if ! test -d scanmots;then mkdir scanmots;fi;");      
//   unsigned int totlen(0),totlentb(0);
//   //random_shuffle(regs.begin(),regs.end());
//   for (ivstring is=regs.begin();is!=regs.end();is++){
//      Sequence seq;
//      ifstream fseq;
//      fseq.open((*is).c_str());
//      string fseqline;
//      getline(fseq,fseqline);
//      seq.name=fseqline;
////cout << seq.name << endl;
//      stringstream firstline(fseqline);
//      string speciesname;
//      firstline >> speciesname;
//      string chrom;
//      firstline >> chrom;
//      if (chrom!=pchrom){
//         cout << chrom << endl;
//         pchrom=chrom;
//      }
//      seq.chrom=intfromchrom(chrom.substr(3));
//      if (seq.chrom==-1) continue;
//      firstline >> seq.start;
//      firstline >> seq.stop;
//      seq.finame=*is;
//      string dum;
//      seq.species=vint(nbspecies,0);
//      vint dumi;
//      seq.iseqs=vvint(nbspecies,dumi);
//      seq.imaps=vvint(nbspecies,dumi);
//      seq.imapsinv=vvint(nbspecies,dumi);
//      while (!fseq.eof()){
////cout << fseqline << endl;
//         int seqnum=speciestonum(fseqline.substr(1,6));
//         getline(fseq,fseqline);
//         seq.imaps[seqnum]=alignedtomap(fseqline);
//         seq.imapsinv[seqnum]=alignedtorevmap(fseqline);
//         string seqwogap=remgaps(fseqline);
//         if (seqwogap.size()>30){
//            seq.species[seqnum]=1;
//         }
//         seq.iseqs[seqnum]=stringtoint(seqwogap);
//         getline(fseq,fseqline);
//      }
//      seq.nbN=compN(seq.iseqs[0]);
//      seq.nbtb=seq.iseqs[0].size()-seq.nbN;
//
//      scanseqforconsinstances(seq,motsdef);
//
//      fseq.close();
//      
//      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
//         //sort(ivm->refinstances_short.begin(),ivm->refinstances_short.end());
//         ostringstream outfs;
//         outfs << "scanmots/" << ivm->name << "_" << ivm->motscorethr2 << ".dat";
//         ofstream outf;
//         if (totlen==0) outf.open(outfs.str().c_str());
//         else outf.open(outfs.str().c_str(),ios::app);
//
//         for (ivinst ivi=ivm->refinstances_short.begin();ivi!=ivm->refinstances_short.end();ivi++){
//            outf << ivi->site << "\t";
//            outf << chromfromint(ivi->chrom) << "\t";
//            outf << ivi->coord << "\t";
//            outf << ivi->sens << "\t";
//            outf << ivi->score << endl;
//         }
//         outf.close();
//         ivm->refinstances_short.clear();
//      }
//
//      totlen+=seq.nbtb+seq.nbN;
//      totlentb+=seq.nbtb;
//   }
//
//   cout << "Total length of the alignement: " << totlen << " bp, including " << totlentb << " unmasked bp" << endl;
//
//   return;
//}
//
//bool operator<(const Sequence & seqscr1,const Sequence & seqscr2)
//{
//   if (seqscr1.score == seqscr2.score)
//   {
//      return seqscr1.sign > seqscr2.sign; //for identical scores, we want the negatives first (for ROC consistency)
//   }
//   else
//   {
//      return seqscr1.score > seqscr2.score; 
//   };
//}
//
//bool operator<(const Motif & mot1,const Motif & mot2)
//{
//   if (mot1.optauc == mot2.optauc)
//   {
//      return mot1.index < mot2.index;//the lower index first 
//   }
//   else
//   {
//      return mot1.optauc > mot2.optauc;//the best auc first
//   };
//}
//
//   void
//calcscore(Motif & im,Sequence & seqscr)//(vmot & mots)
//{
//   unsigned int ncons=im.nbmatchcons(seqscr);//seqscr.motis[im->index];
//   //   unsigned int ncons=im.nbmatchmat(seqscr.seq);//seqscr.motis[im->index];
//   //   unsigned int ncons=seqscr.motis[im.index];
//   double lseq=seqscr.nbtb;
//   double score=0;
//   //cout << ncons << " ";
//   //score= ncons;//*log(im.lambdatrain/im.lambda) + lseq*(im.lambda-im.lambdatrain);
//   double lnj=0;
//   for (int i=0;i<ncons;i++){
//      lnj+=log(i+1);
//   }
//   seqscr.motis[im.index]=ncons;
//   score=ncons;
//   //score=ncons*log(im.lambdatrain/im.lambda);
//   //score= lnj-ncons*log(im.lambda*lseq) + lseq*im.lambda;
//   //   score= (double)ncons/lseq;//*log(im.lambdatrain/im.lambda) + lseq*(im.lambda-im.lambdatrain);
//   seqscr.score=score;
//   return;
//   //cout << seqscr.score << endl;
//}
//
//   void
//calcscore(vmot & vm,Sequence & seqscr)//(vmot & mots)
//{
//
//   Sequence seqtmp=seqscr;
//   double lseq=seqscr.nbtb;
//   double score=0;
//   for (ivmot im=vm.begin();im!=min(vm.begin()+nbmots_for_score,vm.end());im++){
//      width=im->bsinit.size();
//      unsigned int ncons;
//      if (args_info.maskforscore_given){
//         ncons=im->nbmatchconsnmask(seqtmp);//seqscr.motis[im->index];
//      } else{
//         ncons=im->nbmatchcons(seqtmp);//seqscr.motis[im->index];
//      }
//      //   score+=ncons;
//      //      double lnj=0;
//      //      for (int i=0;i<ncons;i++){
//      //         lnj+=log(i+1);
//      //      }
//      //   //   score+=lnj-ncons*log(im->lambda*lseq) + lseq*im->lambda;
//      //      if (args_info.byaffinity_given){
//      //         score+=im->scorematchcons(seqtmp);
//      //      } else {
//      //         score+=ncons;
//      //      }
//      if (args_info.weightmots_given){
//
//         score+= ncons*log(im->lambdatrain/im->lambda);// + lseq*(im.lambda-im.lambdatrain);
//      }
//      else {
//         score+= ncons;//*log(im->lambdatrain/im->lambda);// + lseq*(im.lambda-im.lambdatrain);
//      }
//      seqscr.motis[im->index]=ncons;
//   }
//   seqscr.score=score;
//   return;
//   //cout << seqscr.score << endl;
//}
//
//   void
//calcscore(vmot & vm,vseq & vs)//(vmot & mots)
//{
//   for (ivseq ivs=vs.begin();ivs!=vs.end();ivs++){
//      calcscore(vm,*ivs);
//   }
//   return;
//   //cout << seqscr.score << endl;
//}
//
//// allncons is a vint where each int is a number of conserved motifs
//   double
//calcscore(vmot & vm,vint & allncons,int lseq)
//{
//
//   double score=0;
//   for (ivmot im=vm.begin();im!=min(vm.begin()+nbmots_for_score,vm.end());im++){
//      unsigned int ncons=allncons[im->index];
//      double lnj=0;
//      for (int i=0;i<ncons;i++){
//         lnj+=log(i+1);
//      }
//      //score+=lnj-ncons*log(im->lambda*lseq) + lseq*im->lambda;
//      //    score+=allncons[im->index];
//      score+= allncons[im->index];//*log(im->lambdatrain/im->lambda);// + lseq*(im.lambda-im.lambdatrain);
//   }
//   return score;
//}
//
//
//svg::svg()
//{
//   xsize=800;
//   ysize=600;
//   xoffset=140;
//   yoffset=65;
//   pos=0;
//};
//
//   void
//svginit(ofstream & svgfile, svg s)
//{
//   svgfile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
//   svgfile << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n";
//   svgfile << "<svg version=\"1.0\" id=\"Calque_1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" x=\"0px\" y=\"0px\"\n";
//   svgfile << "	 width=\""<< s.xsize <<"px\" height=\""<< s.ysize <<"px\" viewBox=\"0 0 "<< s.xsize <<" "<< s.ysize <<
//      "\" enable-background=\"new 0 0 "<< s.xsize <<" "<< s.ysize <<"\" xml:space=\"preserve\">\n";
//
//   svgfile << "<line fill=\"none\" stroke=\"" << "black" << "\" stroke-width=\"" << 3 << "\" x1=\"" << 0 << "\" y1=\"" << 0 << "\" x2=\"" << 4 << "\" y2=\"" << 0 << "\"/>\n";
//}
//
//   void
//svgdisplay(ofstream & svgfile,Sequence & seq, svg & s)
//{
//   double ytext=s.yoffset+300*s.pos;
//   double xbegin=s.xoffset;
//   double xend=xbegin+0.4*seq.imaps[0].size();
//
//   svgfile << "<text transform=\"matrix(1 0 0 1 55.5 " << ytext << ")\" font-size=\"12\">" << seq.finame << "</text>\n";
//
//   for (unsigned int i=0;i<nbspecies;i++){
//      if (seq.species[i]){
//
//         double yline=s.yoffset+20+300*s.pos+20*i;
//
//         string dro;
//         dro = seq.speciesname[i];
//         svgfile << "<text transform=\"matrix(1 0 0 1 40 " << yline << ")\" font-size=\"12\">" << dro << "</text>\n";
//         svgfile << "<line fill=\"none\" stroke=\"#000000\" stroke-width=\"0.5\" x1=\"" << xbegin << "\" y1=\"" << yline << "\" x2=\"" << xend << "\" y2=\"" << yline << "\"/>\n";
//
//         vvint coords;
//         vint curcoord;
//         int seqorali=0;
//         unsigned int p=0;
//         for (istring is=seq.seqsrealigned[i].begin();is!=seq.seqsrealigned[i].end();is++){
//            if (*is=='-'){
//               if (seqorali==1){
//                  curcoord.push_back(p);
//                  coords.push_back(curcoord);
//                  curcoord.clear();
//                  seqorali=0;
//               }
//            } else {
//               if (seqorali==0){
//                  curcoord.push_back(p);
//                  seqorali=1;
//               }
//               if (seqorali==1 && is==seq.seqsrealigned[i].end()-1){
//                  curcoord.push_back(p);
//                  coords.push_back(curcoord);
//                  curcoord.clear();
//               }
//            }
//            p++;
//         }
//
//         for (ivvint ivv=coords.begin();ivv!=coords.end();ivv++){
//            svgfile << "<line fill=\"none\" stroke=\"#000000\" stroke-width=\"3\" x1=\"" << xbegin+0.4*(*ivv)[0] << 
//               "\" y1=\"" << yline << "\" x2=\"" << xbegin+0.4*(*ivv)[1] << "\" y2=\"" << yline << "\"/>\n";
//         }
//      }
//   }
//
//   for (ivinstseq ivi=seq.instances.begin();ivi!=seq.instances.end();ivi++){
//      int moti=(*ivi).motindex;
//      if (moti<8){
//         string color;
//         if (moti==0){
//            color="red";
//         } else if (moti==1){
//            color="blue";
//         } else if (moti==2){
//            color="green";
//         } else if (moti==3){
//            color="purple";
//         } else if (moti==4){
//            color="gray";
//         } else if (moti==5){
//            color="orange";
//         } else if (moti==6){
//            color="brown";
//         } else if (moti==7){
//            color="gold";
//         }
//         double xmot=xbegin+0.4*(*ivi).pos;
//         double yline=s.yoffset+20+300*s.pos+20*(*ivi).species;
//         string width;
//         if ((*ivi).score>scorethr2) width="3";
//         else width="1";
//         svgfile << "<line fill=\"none\" stroke=\"" << color << "\" stroke-width=\"" << width << "\" x1=\"" << xmot << "\" y1=\"" << yline-5 << "\" x2=\"" << xmot << "\" y2=\"" << yline+5 << "\"/>\n";
//         svgfile << "<text transform=\"matrix(1 0 0 1 " << xmot-2 << " " << yline-8 << ")\" font-size=\"8\">" << fixed << setprecision(1) << (*ivi).score << "</text>\n";
//      }
//   }
//   s.pos++;
//}
//
//
//   void
//svgclose(ofstream & svgfile)
//{
//   svgfile << "</svg>\n";
//}
//
//   void
//scanseqsforsvg(vseq & align)
//{
//
//   system("if ! test -d display;then mkdir display;fi;");      
//   system("if ! test -d display/svg;then mkdir display/svg;fi;");      
//
//   for (ivseq is=align.begin();is!=align.end();is++){
//
//      int xsize(0);
//      int ysize(0);
//      svg s;
//      string filename("display/svg/");
//      filename+=(*is).name;
//      filename+=".svg";
//      ofstream svgfile(filename.c_str());
//
//      scanseq(*is,motsdef);
//
//      //we set the size for the svg file
//      xsize=s.xoffset+(int)(0.4*(*is).imaps[0].size());
//      if (xsize>s.xsize) s.xsize=xsize;
//      s.xsize+=s.xoffset;
//      ysize=s.yoffset+340; 
//      s.ysize=ysize+s.yoffset;
//
//      svginit(svgfile,s);
//      svgdisplay(svgfile,*is,s);
//      svgclose(svgfile);
//
//      //cout << s.xsize << " " << s.ysize << endl; 
//   }
//}
//
//   void
//outputresults(unsigned int nbmots_score)
//{
//   cout << "Shuffling" << endl;
//   random_shuffle(potregs.begin(),potregs.end());
//   for (ivginst ivg=potregs.begin();ivg!=potregs.end();ivg++){
//      if (args_info.weightmots_given){
//         (*ivg).compscoreweight(motsdef,nbmots_score);
//      }
//      else{
//         (*ivg).compscore(motsdef,nbmots_score);
//      }
//   }
//   cout << "Sorting" << endl;
//   sort(potregs.begin(),potregs.end());
//   cout << "Discarding" << endl;
//   vstring gnames;
//   vginst potregs_def;
//
//   // Keep best scoring enhancer per gene
//   // OR
//   // Keep all enhancers per gene, removing overlapping low score ones
//   if (args_info.discard_on_gene_names_given){
//      cout << "discard on gene" << endl;
//      int test=0;
//      for (ivginst ivg=potregs.begin();ivg!=potregs.end();ivg++){
//         test=0;
//         for (ivstring ivs=gnames.begin();ivs!=gnames.end();ivs++){
//            if ((*ivg).besttss.gene==*ivs){
//               (*ivg).discarded=1;
//               test=1;
//               break;
//            }
//         }
//         if (test==0){
//            gnames.push_back((*ivg).besttss.gene);
//         }
//      }
//   } 
//   else {
//      for (ivginst ivg=potregs.begin();ivg!=potregs.end();ivg++){
//         int test=0;
//         for (ivginst ivg2=potregs_def.begin();ivg2!=potregs_def.end();ivg2++){
//            if ((*ivg).distance(*ivg2)<scanwidth){
//               test=1;
//               break;
//            }
//         }
//         if (test==0){
//            potregs_def.push_back(*ivg);
//         }
//      }
//   }
//   cout << "Displaying" << endl;
//   stringstream filename;
//   filename << "result" << nbmots_score << ".dat";
//   ofstream res(filename.str().c_str());
//   for (ivginst ivg=potregs_def.begin();ivg!=potregs_def.end();ivg++){
//      if (!(*ivg).discarded){
//         res << (*ivg).score << " " << chromfromint((*ivg).chrom) << ":" << (*ivg).start << ".." << (*ivg).stop << " ";
//         res << (*ivg).besttss.gene << " ";
//         for (ivTSS ivt=(*ivg).TSSs.begin();ivt!=(*ivg).TSSs.end();ivt++){
//            res << (*ivt).gene << ";";
//         }
//         res << " ";
//         for (ivint ivi=(*ivg).nbmots.begin();ivi!=(*ivg).nbmots.end();ivi++){
//            res << *ivi << " ";
//         }
//         res << "\n";
//         //         for (ivinst ivi=(*ivg).instances.begin();ivi!=(*ivg).instances.end();ivi++){
//         //            cout << "\t" << chromfromint((*ivi).chrom) << ":" << (*ivi).coord << "\n";
//         //         }
//      }
//   }
//   res << "\n";
//   res.close();
//
//   //	stringstream filenameseqs;
//   //	filenameseqs << "seqs-" << nbmots_score << ".dat";
//   //	ofstream seqs(filenameseqs.str().c_str());
//   //	unsigned int count=0;
//   //	for (ivginst ivg=potregs_def.begin();ivg!=potregs_def.end();ivg++){
//   //		if (!(*ivg).discarded){
//   //			unsigned int mini=scanwidth;
//   //			unsigned int max=0;
//   //			seqs << "chrom : " << (*ivg).chrom << "\n";
//   //			seqs << "start : " << (*ivg).start << "\n";
//   //			seqs << "stop : " << (*ivg).stop << "\n";
//   //			for (ivinst ivi=(*ivg).instances.begin();ivi!=(*ivg).instances.end();ivi++){
//   //				if ((*ivi).coord<mini) mini=(*ivi).coord;
//   //				if ((*ivi).coord>max) max=(*ivi).coord;
//   //			}
//   //			seqs << "mini : " << mini << "\n";
//   //			seqs << "max : " << max << "\n";
//   //			unsigned int start=(mini+max)/2-scanwidth/2;
//   //			unsigned int stop=(mini+max)/2+scanwidth/2;
//   //			for (ivinst ivi=(*ivg).instances.begin();ivi!=(*ivg).instances.end();ivi++){
//   //				seqs << "motif " << (*ivi).motindex << " at position " << (*ivi).coord-start << "\n";;
//   //			}
//   //			seqs << ">seq_" << count << " annotated " << (*ivg).besttss.gene;
//   //			seqs << " ";
//   //			seqs << "\n";
//   //			seqs << chromints[(*ivg).chrom].seq.substr(start,stop-start) << "\n";
//   //			count++;
//   //			if (count>100) break;
//   //		}
//   //	}
//   //	seqs << "\n";
//   //	seqs.close();
//
//   if (args_info.phenotype_given){
//      stringstream filehist;
//      filehist << "hist" << nbmots_score << ".dat";
//      ofstream histo(filehist.str().c_str());
//      displayhist(potregs_def,histo);
//      histo.close();
//      if (args_info.print_histo_sets_given){
//         stringstream filehist_back;
//         filehist_back << "hist-back" << nbmots_score << ".dat";
//         ofstream histo_back(filehist_back.str().c_str());
//         displayhist_set(potregs_def,gbacks,histo_back);
//         histo_back.close();
//
//         stringstream filehist_interest;
//         filehist_interest << "hist-interest" << nbmots_score << ".dat";
//         ofstream histo_interest(filehist_interest.str().c_str());
//         displayhist_set(potregs_def,phenos,histo_interest);
//         histo_interest.close();
//      }
//   }
//
//}
//
//
//
//
//vseq regtests;
//vseq regints;
//
//   void
//extracttofastawfullname(string folder)
//{
//   ofstream outf;
//   for (ivseq iv=regints.begin();iv!=regints.end();iv++){
//      Sequence seq=*iv;
//      if (seq.species[0] && seq.nbtb>0){ 
//         stringstream file;
//         file << folder;
//         file << seq.name << "_";
//         file << chromfromint(seq.chrom) << "_";
//         file << seq.start << "_";
//         file << seq.stop << ".fa";
//         outf.open(file.str().c_str());
//         Sequence & s=seq;
//         for (int i=0;i<nbspecies;i++){
//            if (s.species[i]){
//               if (i==0) outf << ">" << numtospecies(i) << " " <<
//                  "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
//               else  outf << ">" << numtospecies(i) << endl;
//               outf << s.seqsrealigned[i] << endl;
//            }; 
//         }
//         outf.close();
//      }
//   }
//}
//
//
//
//   double
//distcv(vvd& mat1, vvd& mat2)
//{
//   double max(0);
//   int col(0);
//   for (ivvd iv=mat1.begin();iv!=mat1.end();iv++){
//      int line(0);
//      for (ivd il=(*iv).begin();il!=(*iv).end();il++){
//         double diff;
//         diff=fabs((*il)-mat2[col][line]);
//         if (diff>max) max=diff;
//         line++;
//      }
//      col++;
//   }
//   return max;
//}
//
//   void
//matflank(Motif & mot,int numflank)
//{
//   vd dum(4,0);
//   vvd flank_left(numflank,dum);
//   vvd flank_right(numflank,dum);
//
//   for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){
//      vint vdum;
//      for (int i=1;i<=numflank;i++){
//         if ((*ivm).strand==1){
//            ivint iv_left=(*ivm).matchespos[0]-i;
//            if (iv_left>(*ivm).seq_start) flank_left[numflank-i][(*iv_left)]+=1;
//            ivint iv_right=(*ivm).matchespos[0]+width+i;
//            if (iv_right<(*ivm).seq_stop) flank_right[i-1][(*iv_right)]+=1;
//         }
//         else if ((*ivm).strand==-1){
//            ivint iv_left=(*ivm).matchespos[0]+width+i;
//            if (iv_left<(*ivm).seq_stop){
//               int brc;
//               if (*iv_left==0) brc=1;
//               if (*iv_left==1) brc=0;
//               if (*iv_left==2) brc=3;
//               if (*iv_left==3) brc=2;
//               flank_left[numflank-i][brc]+=1;
//            }  
//            ivint iv_right=(*ivm).matchespos[0]-i;
//            if (iv_right>(*ivm).seq_start){ 
//               int brc;
//               if (*iv_right==0) brc=1;
//               if (*iv_right==1) brc=0;
//               if (*iv_right==2) brc=3;
//               if (*iv_right==3) brc=2;
//               flank_right[i-1][brc]+=1;
//            }  
//         }
//      }
//
//      for (ivint iv=(*ivm).matchespos[0];iv!=(*ivm).matchespos[0]+width;iv++){
//         vdum.push_back(*iv);
//      }
//      // cout << flank_left << "\t" << vinttostring(vdum) << "\t" << flank_right << endl; 
//   }
//
//   alpha=conca;
//   beta=concc;
//   countfreq(flank_left);
//   freqtolog(flank_left);
//   countfreq(flank_right);
//   freqtolog(flank_right);
//   //   cout << flank_left <<  "\t" << flank_right << endl; 
//
//   vvd finmat;
//   for (ivvd iv=flank_left.begin();iv!=flank_left.end();iv++){
//      finmat.push_back(*iv);
//   }
//   for (ivvd iv=mot.matprec.begin();iv!=mot.matprec.end();iv++){
//      finmat.push_back(*iv);
//   }
//   for (ivvd iv=flank_right.begin();iv!=flank_right.end();iv++){
//      finmat.push_back(*iv);
//   }
//
//   width+=2*numflank;
//   mot.matprec=finmat;
//   mot.matprecrevcomp=reversecomp(finmat);
//
//   return;
//}
//
//
//   void
//seqanalysis(Sequence & currseq,vmot & genmots)
//{
//   unsigned int i=0;
////      for (int j=0;j<nbspecies;j++){
////      cout << currseq.iseqs[j] << endl;
////      }
//   for (vint::iterator istr=currseq.iseqs[0].begin();istr!=currseq.iseqs[0].end()-width+1;istr++){
//      //cout << "\r" << i+1 << "/" << currseq.iseqs[0].size()-width+1 ; 
//      vint bs(istr,istr+width);
//     //cout << i << " " << bs << endl;
//      if (compN(bs)>0) continue;
//      Motif currmot;
//      currmot.bsinit=vinttostring(bs);
//      currmot.seqinit=currseq.name;
//      currmot.pos=i;
//      motiftomat(bs,currmot);
//      currmot.matricerevcomp=reversecomp(currmot.matrice);
//      currmot.matprec=currmot.matrice;
//      currmot.matprecrevcomp=currmot.matricerevcomp;
//      vvd pmat=currmot.matprec;
//      unsigned int nbconv(0);
//      for (int nb=1;nb<=nbiter;nb++){
//         double max=0.01;
//         int iter(0);
//         while(max>0){
//            if (nb>2) currmot.matinit(scorethr2);
//            else currmot.matinit(scorethr);
//            if (currmot.nbmot<1) break;
//            
//            currmot.compprec();
//            max=distcv(currmot.matprec,pmat);
//            pmat=currmot.matprec;
//            iter++;
//            if (iter==20) break;
//            nbconv++;
//         }
//         if (nb==1){
//            currmot.corrprec();
//            pmat=currmot.matprec;
//            currmot.matprecrevcomp=reversecomp(currmot.matprec);
//         }
//      }
//      //cout << currmot.nbmot << " " ; 
//      //cout.flush();
//
//      currmot.matinit(scorethr2);
//      if (currmot.nbmot>2){
//         currmot.pvaluecomp();
//         //currmot.display(streamfile);
//         //		for (ivma ima=currmot.seqs.begin();ima!=currmot.seqs.end();ima++){
//         //  cout << (*ima).alignseq[0] << endl;
//         //		}
//         //cout << currmot.matprec << endl;
//         genmots.push_back(currmot);
//      }
//      i++;
//   }
//   cout << endl;
//}
//   
//vTSS TSSall;
//
//   void
//loadannots()
//{
//   ifstream annots;
//   annots.open(args_info.phenotype_arg);
//
//   back_insert_iterator<vstring> dest(phenos);
//   copy(iisstring(annots),iisstring(),dest);
//
//   annots.close();
//
//   ifstream glist;
//   //!!!! FIRST scangens were done with genelist-wosensory
//   //   if (species==1) glist.open("/home/santolin/these/files/droso/plaza/phenos/neg-ovo-pheno.dat");
//   //   else if (species==2) glist.open("/home/santolin/these/files/mus/affymetrix/e10/genelist-wo-pos-down-uniq.dat");
//   if (args_info.phenoback_given){
//      glist.open(args_info.phenoback_arg);
//   }
//   else {
//      if (species=="droso") glist.open("/home/rouault/these/sequence/genomes/genelist.dat");
//      else if (species=="eutherian") glist.open("/home/santolin/these/files/mus/biomart/genelist-protein-coding+miRNA.dat");
//   }
//
//   back_insert_iterator<vstring> destg(gbacks);
//   copy(iisstring(glist),iisstring(),destg);
//
//   glist.close();
//
//   for (ivstring ivs=phenos.begin();ivs!=phenos.end();ivs++){
//      for (ivstring ivs2=gbacks.begin();ivs2!=gbacks.end();ivs2++){
//         if (*ivs2==*ivs){
//            gbacks.erase(ivs2);
//            break;
//         }
//      }
//   }
//
//   //   for (ivstring ivs=phenos.begin();ivs!=phenos.end();ivs++){
//   //      cout << *ivs << endl;
//   //   }
//}
//
//   vcoord
//loadcoords()
//{
//   ifstream inf;
//   inf.open(args_info.coord_file_arg);
//   vcoord coords;
//   back_insert_iterator<vcoord> dest(coords);
//   copy(iiscoord(inf),iiscoord(),dest);
//   return coords;
//}
//
//
//   void
//loadchroms()
//{
//   ifstream chromsf;
//   if (species=="droso"){
//      chromsf.open("/home/rouault/these/sequence/genomes/melano-only/dmel-all-chromosome-r4.3.fasta");
//   } else if (species=="eutherian"){
//      chromsf.open("/home/santolin/these/files/mus/genome/fasta/all/all.chromosome_no_MT.fasta");
//      //chromsf.open("/home/santolin/these/files/mus/genome/fasta/all/chr4.fa");
//      //chromsf.open("/home/santolin/these/files/mus/genome/fasta/repeat_masked/all.chromosome_no_MT.fa");
//   }
//   back_insert_iterator<vchrom> dest(chromints);
//   copy(iischrom(chromsf),iischrom(),dest);
//   sort(chromints.begin(),chromints.end());
//
//   for (int i=0;i<chromints.size();i++){
//      cout << "chrom " << chromints[i].name << " ";
//      cout << "of length " << chromints[i].seq.size() << " : ";
//      cout << chromints[i].seq.substr(0,50) << " ... ";
//      cout << chromints[i].seq.substr(chromints[i].seq.size()-50);
//      cout << "\n";
//   }
//
//
//   chromsf.close();
//}
//
//   vint 
//loadlengthchrom()
//{
//
//   vint lchr(nbchrom,0);
//
//   if (species=="droso"){ // *** to be corrected to be easily updated
//
//      lchr[0]=22410834;//chr2L
//      lchr[1]=20769785;//chr2R
//      lchr[2]=23774897;//chr3L
//      lchr[3]=27908053;//chr3R
//      lchr[4]=1284640;//chr4
//      lchr[5]=22227390;//chrX
//
//
//   } else if (species=="eutherian"){
//
//      lchr[0]=197195432;
//      lchr[1]=181748087;
//      lchr[2]=159599783;
//      lchr[3]=155630120;
//      lchr[4]=152537259;
//      lchr[5]=149517037;
//      lchr[6]=152524553;
//      lchr[7]=131738871;
//      lchr[8]=124076172;
//      lchr[9]=129993255;
//      lchr[10]=121843856;
//      lchr[11]=121257530;
//      lchr[12]=120284312;
//      lchr[13]=125194864;
//      lchr[14]=103494974;
//      lchr[15]=98319150;
//      lchr[16]=95272651;
//      lchr[17]=90772031;
//      lchr[18]=61342430;
//      lchr[19]=166650296;//X
//      lchr[20]=15902555;//Y
//   }
//
//   return lchr;
//
//}
//
//   void
//initgroupedinst()
//{
//   cout << "Loading length chroms..." << endl;
//   lengthchrom=loadlengthchrom();
//
//   vginst dumvginst;
//   groupedinst=vvginst(nbchrom,dumvginst);
//   for (unsigned int i=0;i<nbchrom;i++){
//      for (unsigned int j=0;j<lengthchrom[i]/scanstep+1;j++){
//         groupedinst[i].push_back(GroupInstance(j*scanstep,j*scanstep+scanwidth,i));
//      }
//   }
//   ifstream annots;
//   if (species=="droso"){
//      annots.open("/home/rouault/these/sequence/genomes/regres-wellform-all.dat");
//   } else if (species=="eutherian"){
//      annots.open("/home/santolin/these/files/mus/biomart/genes-n-strand-protein-coding+miRNA-no-MT-n-NT-TSS.dat");
//   }
//   importTSS(TSSall,annots);
//   annots.close();
//   // assign CRMs to genes in an annotextent region
//   cout << "Assign CRMs to genes in annotextend..." << endl;
//   for (ivTSS ivt=TSSall.begin();ivt!=TSSall.end();ivt++){
//      int start=0;
//      if (((*ivt).coord-annotextent-scanwidth/2)/scanstep+1>0) start=((*ivt).coord-annotextent-scanwidth/2)/scanstep;
//      int stop=((*ivt).coord+annotextent-scanwidth/2)/scanstep+1;
//      if (stop>lengthchrom[(*ivt).chrom]/scanstep+1) stop=lengthchrom[(*ivt).chrom]/scanstep+1;
//      for (unsigned int i=start;i<stop;i++){
//         groupedinst[(*ivt).chrom][i].TSSs.push_back(*ivt);
//         //       cout << (*ivt).gene << endl;
//      }
//   }
//   // assign nearest TSS to genes
//   cout << "Assign CRMs to nearest gene..." << endl;
//   for (unsigned int i=0;i<nbchrom;i++){
//      for (ivginst ivg=groupedinst[i].begin();ivg!=groupedinst[i].end();ivg++){
//         (*ivg).compbestannot();
//      }
//   }
//
//   // attribute phenotype to CRMs
//   cout << "Attribute phenotype to CRM..." << endl;
//   for (unsigned int i=0;i<nbchrom;i++){
//      for (ivginst ivg=groupedinst[i].begin();ivg!=groupedinst[i].end();ivg++){
//         (*ivg).isdiscarded();
//         if (!(*ivg).discarded){
//            for (ivstring ivs=phenos.begin();ivs!=phenos.end();ivs++){
//               if ((*ivg).besttss.gene==*ivs){
//                  (*ivg).goodpheno=1;
//                  break;
//               } else {
//                  (*ivg).goodpheno=0;
//               }
//            }
//         } else {
//            (*ivg).goodpheno=0;
//         }
//      }
//   }
//
//}
//
//   void
//getsites(ifstream & file, Sequence & bds)
//{
//   string dum;
//   string dumtest;
//   int i=0;
//   cout << "seq " << endl;
//   getline(file,dum);
//   while(!file.eof()){
//      bds.sitenames.push_back(dum.substr(1));
//      getline(file,dum);
//      if (i>0 && dum.length()!=dumtest.length()) {
//         cout << dum << endl;
//         cout << "Error in site length" << endl;
//         exit(2);
//      }
//      dumtest=dum;
//      i++;
//      bds.seqs.push_back(dum);
//      cout << dum << endl;
//      bds.iseqs.push_back(stringtoint(dum));
//      getline(file,dum);
//   }
//}
//
//   void
//disptex(vseq & seqs)
//{
//   system("if ! test -d display;then mkdir display;fi;");      
//   //   system("if ! test -d display/tex;then mkdir display/tex;fi;");      
//
//   double scoreinit=scorethr2;
//   string filename("display/");
//   filename+="results.tex";
//   ofstream outf(filename.c_str());
//   disptexinit(outf);
//   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){ 
//      cout << "Scanning " << (*ivs).name << endl;
//      dispseqwmots(*ivs,motsdef,outf,scoreinit);
//   }
//   disptexclose(outf);
//
//   outf.close();
//}
//
//   void
//disptexwgaps(vseq & align)
//{
//   system("if ! test -d display;then mkdir display;fi;");      
//
//   cout << "Scanning sequences for instances..." << endl;
//   scanseqsforinstancesnmask(align,motsdef);
//
//   cout << "Defining conserved instances..." << endl;
//   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
//      ivs->instances2instancescons();
//      //cout << ivs->name << "\n" << ivs->instancescons;
//   }
//
//   string folder("display/");
//   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){ 
//      stringstream file;
//      file << folder;
//      file << ivs->name << "_";
//      file << chromfromint(ivs->chrom) << "_";
//      file << ivs->start << "_";
//      file << ivs->stop;
//      file << ".tex";
//      ofstream outf(file.str().c_str());
//      disptexinit(outf);
//      cout << "Scanning " << (*ivs).name << endl;
//      dispseqwmotswgaps(*ivs,outf);
//      disptexclose(outf);
//      outf.close();
//   }
//
//}
//
//   Sequence
//findbestseq(vmot & vm, Sequence & seq)
//{
//
//   Sequence bestseq=seq;
//
//   //Find the best 1kb piece (if possible!) maximizing the sites score function
//   bestseq.score=-1.e6;
//   vint bestmotis(nbmots_for_score,0);
//
//   //we define the 1kb piece with best score
//   for (ivvinstseq iv=seq.instancescons.begin();iv!=seq.instancescons.end();iv++){
//
//      vinstseq vic=*iv;
//      vint bestmotistmp(nbmots_for_score,0);
//      int start=seq.imaps[0][vic[0].pos]-neighbext;
//      if (start<0) start=0;
//      int stop=start+scanwidth-1;
//      if (stop>seq.stop-seq.start){
//         stop=seq.stop-seq.start;
//         start=stop-scanwidth+1;
//         if (start<0) start=0;
//      }
//      for (ivvinstseq iv1=seq.instancescons.begin();iv1!=seq.instancescons.end();iv1++){
//         vinstseq vic1=*iv1;
//         int start1=seq.imaps[0][vic1[0].pos];
//         if (start1>=start && start1<stop-width+1){
//            bestmotistmp[vic1[0].motindex]++;
//         }
//      }
//
//      double score;
//      int lseq=stop-start+1;
//      score=calcscore(vm,bestmotistmp,lseq);
//
//      if (score>bestseq.score){
//         bestseq.score=score;
//         bestseq.start=seq.start+start;
//         bestseq.stop=seq.start+stop;
//         bestmotis=bestmotistmp;
//      }
//   }
//
//   if (seq.instancescons.size()==0) bestseq.score=calcscore(vm,bestmotis,seq.nbtb);
//
//   //   cout << bestseq.name << " " << bestseq.score  <<" " << bestseq.start << " "  <<bestseq.stop << endl;
//   //   cout << seq.name << " "  << seq.start << " "  <<seq.stop << endl;
//
//   //cut original sequence into best piece
//   int rawstart,rawstop;
//   // Coorects minor bug due to extraction (incorrect stop due to gap etc...)
//   if (bestseq.stop-seq.start+1>=seq.imapsinv[0].size()) bestseq.stop=seq.start+seq.imapsinv[0].size()-1;
//   rawstart=seq.imapsinv[0][bestseq.start-seq.start];
//   rawstop=seq.imapsinv[0][bestseq.stop-seq.start];
//   for (int i=0;i<nbspecies;i++){
//      if (seq.species[i]){
//         civint istart=seq.iseqs[i].begin()+seq.imaps[i][rawstart];
//         civint istop=seq.iseqs[i].begin()+seq.imaps[i][rawstop]+1;
//         if (istop>=seq.iseqs[i].end()) istop=seq.iseqs[i].end();
//         vint itmp(istart,istop);
//         bestseq.iseqs[i]=itmp;
//         bestseq.seqs[i]=vinttostring(itmp);
//         bestseq.seqsrealigned[i]=seq.seqsrealigned[i].substr(rawstart,rawstop-rawstart+1);
//         bestseq.imaps[i]=alignedtomap(bestseq.seqsrealigned[i]);
//         bestseq.imapsinv[i]=alignedtorevmap(bestseq.seqsrealigned[i]);
//         //         cout << bestseq.iseqs[i] << endl;
//         //         cout << bestseq.seqs[i] << endl;
//         //         cout << bestseq.seqsrealigned[i] << endl;
//         int nbtb=bestseq.iseqs[i].size()-compN(bestseq.iseqs[i]);
//         if (nbtb>width+neighbext){//>5){}
//            bestseq.species[i]=1;
//         } else {
//            bestseq.species[i]=0;
//         }
//      }
//   }
//   bestseq.nbN=compN(bestseq.iseqs[0]);
//   bestseq.nbtb=bestseq.iseqs[0].size()-bestseq.nbN;
//   //   cout << bestseq.nbN << endl;
//   //   cout << bestseq.nbtb << endl;
//
//   return bestseq;
//
//}
//
//   vseq 
//findbestseqs(vmot & vm, vseq & seqs1)
//{
//   // INIT
//   vseq seqs=seqs1;
//
//   // We need different enhancers of the same gene to have same name
//
//   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
//      if ((*ivs).name.find('-')!=string::npos)
//      {
//         (*ivs).name.erase((*ivs).name.begin()+(*ivs).name.find('-'),(*ivs).name.end());
//      }
//   }
//
//   // SCANNING OLD SEQS
//   //
//   //   cout << "Scanning sequences for instances..." << endl;
//   scanseqsforinstancesnmask(seqs,vm);
//
//   // KEEPING BEST SEQS
//   //
//   //   cout << "Defining conserved instances..." << endl;
//   //For each sequence find best part and cut it into vtmp
//   vseq vtmp;
//   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
//      ivs->instances2instancescons();
//      vtmp.push_back(findbestseq(vm,*ivs));
//   }
//
//   // UNICITY: ONE SEQUENCE PER GENE
//   //
//   vseq vbest;
//   // only keep the best seq among all surrounding / intronic regions of a given gene
//   if (vtmp.size()>1){
//      Sequence bestseq=*(vtmp.begin());
//      for (ivseq ivs=vtmp.begin()+1;ivs!=vtmp.end();ivs++){ 
//         if (ivs->name==bestseq.name){ 
//            if (ivs->score>bestseq.score){
//               bestseq=*ivs;
//            }
//         } else{
//            vbest.push_back(bestseq);
//            bestseq=*ivs;
//            if (ivs==vtmp.end()-1) vbest.push_back(*ivs);
//         }
//      }
//   } else {
//      vbest=vtmp;
//   }
//
//   // RE-SCANNING BEST SEQUENCES
//   //
//   for (ivseq iv=vbest.begin();iv!=vbest.end();iv++){
//      iv->instances.clear();
//      iv->instancescons.clear();
//      vint dum(nbmots_for_score,0);
//      iv->motis=dum;
//   }
//   //   cout << "Scanning best 1kbs for instances..." << endl;
//   //scanseqsforinstancesnmask(vbest,vm);
//   scanseqsforinstances(vbest,vm);
//
//   //   cout << "Defining best 1kbs conserved instances..." << endl;
//   for (ivseq ivs=vbest.begin();ivs!=vbest.end();ivs++){
//      ivs->instances2instancescons();
//      //   cout << ivs->name << endl;
//      //  cout << ivs->instancescons << endl;
//   }
//
//
//   // RE-DEFINE SIZEPOS AND SIZENEG
//   vseq vscorepos;
//   sizepos=0;
//   for (ivseq ivs=vbest.begin();ivs!=vbest.end();ivs++){
//      if (ivs->sign==1){
//         sizepos++;
//      }
//   }
//   sizeneg=vbest.size()-sizepos;
//
//   return vbest;
//}
//
//   void
//dispbestseqs(vseq & seqs, ofstream & outf)
//{
//   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
//      outf << chromfromint(ivs->chrom) << "\t";
//      outf << ivs->start << "\t";
//      outf << ivs->stop << "\t";
//      outf << ivs->name << "\n";
//   }
//}
//
//
//   void
//loadseqsforscore(vseq & vscore)
//{
//
//   if (args_info.posalign_given){ 
//      ifstream inf;
//      inf.open(args_info.posalign_arg);
//      vseq seqs;
//      seqs=loadsequencesconserv(inf);
//      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
//         if (ivs->species[0] && ivs->nbtb>0){
//            ivs->sign=1;
//            vscore.push_back(*ivs);
//         }
//      }
//      sizepos=vscore.size();
//      cout << "Positive set: " << sizepos << "\n";
//      inf.close();
//   }
//
//   if (args_info.poscoords_given){
//      //Working on coordinate files
//      ifstream align;
//      if (species=="droso") align.open("/home/santolin/these/files/droso/align/all/align-masked.dat"); // *** To be changed...
//      else if (species=="eutherian") align.open("/home/santolin/these/files/mus/epo/align-masked.dat");
//
//      //POSITIVES
//      ifstream inf;
//      inf.open(args_info.poscoords_arg);
//      alignscoord=loadcoordconserv(align);
//      align.close();
//      vcoord pcoords;
//      back_insert_iterator<vcoord> pdest(pcoords);
//      copy(iiscoord(inf),iiscoord(),pdest);
//      for (ivcoord ivc=pcoords.begin();ivc!=pcoords.end();ivc++){
//         Sequence seqtoimport=coordtoseq(*ivc);
//         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
//            seqtoimport.sign=1;
//            vscore.push_back(seqtoimport);
//         }
//      }
//      inf.close();
//      sizepos=vscore.size();
//      cout << "Positive set: " << sizepos << "\n";
//   }
//
//   if (args_info.negalign_given){
//      ifstream inf;
//      inf.open(args_info.negalign_arg);
//      vseq seqs;
//      seqs=loadsequencesconserv(inf);
//      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
//         if (ivs->species[0] && ivs->nbtb>0){
//            ivs->sign=-1;
//            vscore.push_back(*ivs);
//         }
//      }
//      sizeneg=vscore.size()-sizepos;
//      cout << "Negative set: " << sizeneg << "\n";
//      inf.close();
//   }
//
//   if (args_info.negcoords_given){
//      ifstream align;
//      if (species=="droso") align.open("/droso/align/all/align-masked.dat");
//      else if (species=="eutherian") align.open("/mus/epo/align-masked.dat");
//
//      //NEGATIVES
//      ifstream inf;
//      inf.open(args_info.negcoords_arg);
//      alignscoord=loadcoordconserv(align);
//      align.close();
//      vcoord ncoords;
//      back_insert_iterator<vcoord> ndest(ncoords);
//      copy(iiscoord(inf),iiscoord(),ndest);
//      for (ivcoord ivc=ncoords.begin();ivc!=ncoords.end();ivc++){
//         Sequence seqtoimport=coordtoseq(*ivc);
//         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
//            seqtoimport.sign=-1;
//            vscore.push_back(seqtoimport);
//         }
//      }
//      inf.close();
//      sizeneg=vscore.size()-sizepos;
//      cout << "Negative set: " << sizeneg << "\n";
//   }
//
//   if (!(args_info.negcoords_given||args_info.negalign_given)){
//      cout << "Using background as negative" << endl;
//      ifstream inf;
//      if (species=="droso") inf.open("/droso/backreg/regs2000-align-1000.dat"); 
//      else if (species=="eutherian") inf.open("/mus/backreg/regs2000-align-1000.dat");//regs2000.fa");
//      vseq seqs;
//      seqs=loadsequencesconserv(inf);
//      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
//         if (ivs->species[0] && ivs->nbtb>0){
//            ivs->sign=-1;
//            vscore.push_back(*ivs);
//         }
//      }
//      sizeneg=vscore.size()-sizepos;
//      cout << "Negative set: " << sizeneg << "\n";
//      inf.close();
//   }
//
//   if (!(args_info.poscoords_given || args_info.posalign_given)) {
//      cout << "Please give positive file" << endl;
//      exit(1);
//   }
//
//   //   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
//   //      if (!(ivs->sign==-1 && ivs->name.find("reg")!=string::npos)){
//   ////         if ((*ivs).name.find('_')!=string::npos)
//   ////         {
//   ////            if (args_info.display_given || args_info.byscore_given){
//   ////               (*ivs).name.insert((*ivs).name.find('_'),"\\");
//   ////            } else{
//   ////               (*ivs).name.erase((*ivs).name.begin()+(*ivs).name.find('_'),(*ivs).name.end());
//   ////            }
//   ////         }
//   //         //      if ((*ivs).name.find('-')!=string::npos)
//   //         //      {
//   //         //         (*ivs).name.erase((*ivs).name.begin()+(*ivs).name.find('-'),(*ivs).name.end());
//   //         //      }
//   //      }
//   //   }
//
//   return;
//}
//
//   void
//loadseqsforscorenotmasked(vseq & vscore)
//{
//
//   if (args_info.posalign_given){ 
//      ifstream inf;
//      inf.open(args_info.posalign_arg);
//      vseq seqs;
//      seqs=loadsequencesconserv(inf);
//      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
//         if (ivs->species[0] && ivs->nbtb>0){
//            ivs->sign=1;
//            vscore.push_back(*ivs);
//         }
//      }
//      sizepos=vscore.size();
//      cout << "Positive set: " << sizepos << "\n";
//      inf.close();
//   }
//
//   if (args_info.poscoords_given){
//      //Working on coordinate files
//      ifstream align;
//      if (species=="droso") align.open("/home/santolin/these/files/droso/align/all/align.dat");
//      else if (species=="eutherian") align.open("/home/santolin/these/files/mus/epo/align.dat");
//
//      //POSITIVES
//      ifstream inf;
//      inf.open(args_info.poscoords_arg);
//      alignscoord=loadcoordconserv(align);
//      align.close();
//      vcoord pcoords;
//      back_insert_iterator<vcoord> pdest(pcoords);
//      copy(iiscoord(inf),iiscoord(),pdest);
//      for (ivcoord ivc=pcoords.begin();ivc!=pcoords.end();ivc++){
//         Sequence seqtoimport=coordtoseq(*ivc);
//         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
//            seqtoimport.sign=1;
//            vscore.push_back(seqtoimport);
//         }
//      }
//      inf.close();
//      sizepos=vscore.size();
//      cout << "Positive set: " << sizepos << "\n";
//   }
//
//   if (args_info.negalign_given){
//      ifstream inf;
//      inf.open(args_info.negalign_arg);
//      vseq seqs;
//      seqs=loadsequencesconserv(inf);
//      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
//         if (ivs->species[0] && ivs->nbtb>0){
//            ivs->sign=-1;
//            vscore.push_back(*ivs);
//         }
//      }
//      sizeneg=vscore.size()-sizepos;
//      cout << "Negative set: " << sizeneg << "\n";
//      inf.close();
//   }
//
//   if (args_info.negcoords_given){
//      ifstream align;
//      if (species=="droso") align.open("/droso/align/all/align.dat");
//      else if (species=="eutherian") align.open("/mus/epo/align.dat");
//
//      //NEGATIVES
//      ifstream inf;
//      inf.open(args_info.negcoords_arg);
//      alignscoord=loadcoordconserv(align);
//      align.close();
//      vcoord ncoords;
//      back_insert_iterator<vcoord> ndest(ncoords);
//      copy(iiscoord(inf),iiscoord(),ndest);
//      for (ivcoord ivc=ncoords.begin();ivc!=ncoords.end();ivc++){
//         Sequence seqtoimport=coordtoseq(*ivc);
//         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
//            seqtoimport.sign=-1;
//            vscore.push_back(seqtoimport);
//         }
//      }
//      inf.close();
//      sizeneg=vscore.size()-sizepos;
//      cout << "Negative set: " << sizeneg << "\n";
//   }
//
//   if (!(args_info.negcoords_given||args_info.negalign_given)){
//      cout << "Using background as negative" << endl;
//      ifstream inf;
//      if (species=="droso") inf.open("/droso/backreg/regs2000-align-1000.dat"); 
//      else if (species=="eutherian") inf.open("/mus/backreg/regs2000-align-1000.dat");
//      vseq seqs;
//      seqs=loadsequencesconserv(inf);
//      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
//         if (ivs->species[0] && ivs->nbtb>0){
//            ivs->sign=-1;
//            vscore.push_back(*ivs);
//         }
//      }
//      sizeneg=vscore.size()-sizepos;
//      cout << "Negative set: " << sizeneg << "\n";
//      inf.close();
//   }
//
//   if (!(args_info.poscoords_given || args_info.posalign_given)) {
//      cout << "Please give positive file" << endl;
//      exit(1);
//   }
//
//   //   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
//   //      if (!(ivs->sign==-1 && ivs->name.find("reg")!=string::npos)){
//   ////         if ((*ivs).name.find('_')!=string::npos)
//   ////         {
//   ////            if (args_info.display_given || args_info.byscore_given){
//   ////               (*ivs).name.insert((*ivs).name.find('_'),"\\");
//   ////            } else{
//   ////               (*ivs).name.erase((*ivs).name.begin()+(*ivs).name.find('_'),(*ivs).name.end());
//   ////            }
//   ////         }
//   //         //      if ((*ivs).name.find('-')!=string::npos)
//   //         //      {
//   //         //         (*ivs).name.erase((*ivs).name.begin()+(*ivs).name.find('-'),(*ivs).name.end());
//   //         //      }
//   //      }
//   //   }
//
//   return;
//}
//
////Loads align-file (fasta) or coord-file (name/chrom/start/stop)
//   vseq
//loadseqs()
//{
//   vseq seqs;
//
//   if (args_info.align_file_given){
//
//      ifstream inf;
//      inf.open(args_info.align_file_arg);
//      seqs=loadsequencesconserv(inf);
//      inf.close();
//   }
//   else if (args_info.coord_file_given){
//
//      ifstream inf;
//      inf.open(args_info.coord_file_arg);
//
//      ifstream align;
//      if (species=="droso") align.open("/droso/align/all/align-masked.dat");
//      else if (species=="eutherian") align.open("/mus/epo/align-masked.dat");
//
//      alignscoord=loadcoordconserv(align);
//
//      align.close();
//
//      vcoord coords;
//      back_insert_iterator<vcoord> dest(coords);
//      copy(iiscoord(inf),iiscoord(),dest);
//      for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
//         Sequence seqtoimport=coordtoseq(*ivc);
//         //         cout << (*ivc).name << " " << (*ivc).start << endl;
//         //         cout << seqtoimport.name << " " << seqtoimport.start << endl;
//         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
//            seqs.push_back(seqtoimport);
//         }
//         //         for (int i=0;i<nbspecies;i++){
//         //         cout << seqtoimport.iseqs[i] << endl; 
//         //         }
//         //         exit(9);
//      }
//      inf.close();
//   }
//   else {
//      cout << "Error in loadseqs: please give a coord/alignment file. Exiting..." << endl;
//      exit(1);
//   }
//   return seqs;
//}
//
//   void
//args_init()
//{
//   width=args_info.width_arg;
//   scorethr2=args_info.threshold_arg;
//   species=args_info.species_arg;
//   if (species=="droso"){ // *** Shouldn't be hardcoded??
//      conca=0.3; 
//      concc=0.5-conca;
//      nbspecies=12;
//      cout << "Species: droso" << endl;
//   }
//   else if (species=="eutherian"){
//      conca=0.263; 
//      concc=0.5-conca;
//      nbspecies=10;
//      cout << "Species: mus" << endl;
//   }
//   conct=conca;
//   concg=concc;
//   evolutionary_model=args_info.evolutionary_model_arg;
//
//   scorethrcons=scorethr2-1.0;
//   scorethr=scorethr2-1.0;
//
//   neighbext=args_info.neighbext_arg;
//   nbmots_for_score=args_info.nbmots_arg;
//   scanwidth=args_info.scanwidth_arg;
//
//   if (species=="droso") nbchrom=6;
//   else if (species=="eutherian") nbchrom=21;
//
//}
//
//   void
//args_init_scangen()
//{
//   scanwidth=args_info.scanwidth_arg;
//   scanstep=args_info.scanstep_arg;
//
//   nbmots_for_score=args_info.nbmots_arg;
//}
//

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

const char package_name[] = PACKAGE_NAME;
const char package_version[] = PACKAGE_VERSION;

int
cmd_version(int argc, char **argv)
{
	printf("%s version %s\n", package_name, package_version);
	return 0;
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

