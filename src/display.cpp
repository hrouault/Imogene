#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

#include "display_cmdline.h"
#include "const.hpp"
#include "motif.hpp"
#include "sequence.hpp"
#include "tree.hpp"

vmot motsdef;

// general tex header
// TODO: more portable (french etc.)
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

// tex header for a table of conserved motifs per enhancer
// TODO: portability + translation

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

// more generality than previous.
// TODO: keep this one and modify code where needed to use it.
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

// for long tables
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

// display nb of conserved motifs in a tex file, for pos/neg seqs
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

// display nb of conserved motifs (for several motifs) in a tex file, for pos/neg seqs
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

// general tex display for nb of conserved motifs
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

// display TFBS on reference sequences
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

// for TFBS color
   string
colfromint(int i)
{
   if (i==0) return "red";
   if (i==1) return "blue";
   if (i==2) return "green";
   if (i==3) return "yellow";

   return "black";

}

// display TFBS on aligned sequences
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


display_args_info display_args;


/** 
 * ===  FUNCTION  ======================================================================
 *         Name:  cmd_display
 *  Description:  Display results
 * =====================================================================================
 */
int
cmd_display(int argc, char **argv)
{
   if ( display_cmdline_parser(argc, argv, & display_args)!=0)
      exit(1);
   
   return 1;
}		/* -----  end of function extract  ----- */
