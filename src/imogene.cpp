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

