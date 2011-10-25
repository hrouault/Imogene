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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <sys/stat.h>

using namespace std;

#include "display_cmdline.h"
#include "display.hpp"
#include "const.hpp"
#include "motif.hpp"
#include "sequence.hpp"
#include "tree.hpp"

display_args_info display_args;
   
void
dispweblogo(vmot& mots)
{

   unsigned int index=1;
   for ( ivmot ivm=mots.begin();ivm!=mots.end();ivm++ ) {
      if ( ivm->check ) {
         stringstream ss;
         ss << "python " << PYTHON_PATH"/weblogo-display.py ";
         ss << "Motif";
         ss << index << " ";
         ss << concc << " ";
         for (ivvd ivv=ivm->matfreq.begin();ivv!=ivm->matfreq.end();ivv++){
            for (ivd iv=ivv->begin();iv!=ivv->end()-1;iv++){
               ss << *iv << ",";
            }
            ss << *(ivv->end()-1) << " ";
         }
         system(ss.str().c_str());
         index++;
      }
   }
   return;

}

svg::svg()
{
   xsize=800;
   ysize=600;
   xoffset=140;
   yoffset=65;
   pos=0;
};

   void
texify (string & str)
{
   int found=str.find("_");
   unsigned int deca=0;
   while (found!=string::npos){
      str.insert(found,"\\");
      found=str.find("_",found+3);
   }
   found=str.find("#");
   deca=0;
   while (found!=string::npos){
      str.insert(found,"\\");
      found=str.find("#",found+3);
   }
   return;
}

// for TFBS color
   string
colfromint(int i)
{
   if (i==0) return "red";
   if (i==1) return "blue";
   if (i==2) return "green";
   if (i==3) return "yellow";
   if (i==4) return "brown";

   return "black";

}

   void
svginit(ofstream & svgfile, svg s)
{
   svgfile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
   svgfile << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n";
   svgfile << "<svg version=\"1.0\" id=\"Calque_1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" x=\"0px\" y=\"0px\"\n";
   svgfile << "	 width=\""<< s.xsize <<"px\" height=\""<< s.ysize <<"px\" viewBox=\"0 0 "<< s.xsize <<" "<< s.ysize <<
      "\" enable-background=\"new 0 0 "<< s.xsize <<" "<< s.ysize <<"\" xml:space=\"preserve\">\n";

   svgfile << "<line fill=\"none\" stroke=\"" << "black" << "\" stroke-width=\"" << 3 << "\" x1=\"" << 0 << "\" y1=\"" << 0 << "\" x2=\"" << 4 << "\" y2=\"" << 0 << "\"/>\n";
}

   void
svgdisplay(ofstream & svgfile,Sequence & seq, svg & s)
{
   double ytext=s.yoffset+300*s.pos;
   double xbegin=s.xoffset;
   double xend=xbegin+0.4*seq.imaps[0].size();

   svgfile << "<text transform=\"matrix(1 0 0 1 55.5 " << ytext << ")\" font-size=\"12\">" << seq.finame << "</text>\n";

   for (unsigned int i=0;i<nbspecies;i++){
      if (seq.species[i]){

         double yline=s.yoffset+20+300*s.pos+20*i;

         string dro;
         dro = numtospecies(i);
         svgfile << "<text transform=\"matrix(1 0 0 1 40 " << yline << ")\" font-size=\"12\">" << dro << "</text>\n";
         svgfile << "<line fill=\"none\" stroke=\"#000000\" stroke-width=\"0.5\" x1=\"" << xbegin << "\" y1=\"" << yline << "\" x2=\"" << xend << "\" y2=\"" << yline << "\"/>\n";

         vvint coords;
         vint curcoord;
         int seqorali=0;
         unsigned int p=0;
         for (istring is=seq.seqsrealigned[i].begin();is!=seq.seqsrealigned[i].end();is++){
            if (*is=='-'){
               if (seqorali==1){
                  curcoord.push_back(p);
                  coords.push_back(curcoord);
                  curcoord.clear();
                  seqorali=0;
               }
            } else {
               if (seqorali==0){
                  curcoord.push_back(p);
                  seqorali=1;
               }
               if (seqorali==1 && is==seq.seqsrealigned[i].end()-1){
                  curcoord.push_back(p);
                  coords.push_back(curcoord);
                  curcoord.clear();
               }
            }
            p++;
         }

         for (ivvint ivv=coords.begin();ivv!=coords.end();ivv++){
            svgfile << "<line fill=\"none\" stroke=\"#000000\" stroke-width=\"3\" x1=\"" << xbegin+0.4*(*ivv)[0] << 
               "\" y1=\"" << yline << "\" x2=\"" << xbegin+0.4*(*ivv)[1] << "\" y2=\"" << yline << "\"/>\n";
         }
      }
   }

   for (ivinstseq ivi=seq.instances.begin();ivi!=seq.instances.end();ivi++){
      int moti=(*ivi).motindex;
      if (moti<8){
         string color;
         if (moti==0){
            color="red";
         } else if (moti==1){
            color="blue";
         } else if (moti==2){
            color="green";
         } else if (moti==3){
            color="purple";
         } else if (moti==4){
            color="gray";
         } else if (moti==5){
            color="orange";
         } else if (moti==6){
            color="brown";
         } else if (moti==7){
            color="gold";
         }
         double xmot=xbegin+0.4*(*ivi).pos;
         double yline=s.yoffset+20+300*s.pos+20*(*ivi).species;
         string width;
         if ((*ivi).score>scorethr2) width="3";
         else width="1";
         svgfile << "<line fill=\"none\" stroke=\"" << color << "\" stroke-width=\"" << width << "\" x1=\"" << xmot << "\" y1=\"" << yline-5 << "\" x2=\"" << xmot << "\" y2=\"" << yline+5 << "\"/>\n";
         svgfile << "<text transform=\"matrix(1 0 0 1 " << xmot-2 << " " << yline-8 << ")\" font-size=\"6\">" << fixed << setprecision(1) << (*ivi).score << "</text>\n";
      }
   }
   s.pos++;
}

   void
svgclose(ofstream & svgfile)
{
   svgfile << "</svg>\n";
}

   void
scanseqsforsvg(vseq & align,vmot & mots)
{

   for (ivseq is=align.begin();is!=align.end();is++){

      int xsize(0);
      int ysize(0);
      svg s;
      string filename("display/");
      filename+=(*is).name;
      filename+=".svg";
      ofstream svgfile(filename.c_str());
      if (svgfile.fail()){
         cerr << "Cannot open file for svg recording: " << strerror(errno) << endl;
         exit(-1);
      }

      //we set the size for the svg file
      xsize=s.xoffset+(int)(0.4*is->iseqs[0].size());
      if (xsize>s.xsize) s.xsize=xsize;
      s.xsize+=s.xoffset;
      ysize=s.yoffset+340; 
      s.ysize=ysize+s.yoffset;

      svginit(svgfile,s);
      svgdisplay(svgfile,*is,s);
      svgclose(svgfile);
   }
}

// display TFBS on reference sequences
   void
dispseqwmots (Sequence & seq, vmot & mots, ofstream & outf)
{
   //HEADER
   string name=numtospecies(0);
   string texname=seq.name;
   outf << "$>$" << name <<  "\t";
   outf << texname << "\t" ;
   outf << chromfromint(seq.chrom) << "\t";
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
      int motwidth=mots[ivs->motindex].motwidth;
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
         int motwidth=mots[ivs->motindex].motwidth;
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
      for (unsigned int spe=0;spe<1;spe++){
         if(seq.species[spe]){
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
               if ( i>0 && i%10 == 0 ) {
                  outf << "\t";
               }
            }
            outf << "\\\\\n";
            int pstate=0;
            int ppos=start;
            for (unsigned int i=start;i<stop;i++){
               if (vvstate[spe][i]>0  && vvstate[spe][i]!=pstate){
                  int deca;
                  deca=i-ppos;
                  outf << "\\hspace*{" << deca << "\\charwidth}";
                  if (vvstate[spe][i]==2) outf << "\\textit{";
                  outf << "\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                  outf << setprecision(1) << fixed <<  scores[spe][i] ;
                  outf << "}";
                  if (vvstate[spe][i]==2) outf << "}";
                  // to count the number of digits:
                  stringstream score;
                  score << setprecision(1) << fixed <<  scores[spe][i] ;
                  for (int ii=1;ii<=score.str().size();ii++){
                     if ( i>0 && i%10 == 0 ) {
                        outf << "\\hspace*{1\\charwidth}";
                     }
                     i++;
                  }
                  ppos=i;
               }
               pstate=vvstate[spe][i];
               if ( i>0 && i%10 == 0 ) {
                  outf << "\\hspace*{1\\charwidth}";
               }
            }
            outf << "\\\\\n";
         }
      }
      start=stop;
      stop=min(stop+60,sizeseq);
   }
   outf << "}\n";
   outf << "\\\\\n";
}

// display TFBS on reference sequences
   void
dispseqwmots_html (Sequence & seq, vmot & mots, ofstream & outf)
{
   //HEADER
   outf << "<h3>";
   string name=numtospecies(0);
   outf << ">" << name <<  " ";
   outf << seq.name << " " ;
   outf << chromfromint(seq.chrom) << " ";
   outf << seq.start << " " ;
   outf << seq.stop << "</h3>" << endl;

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
      int motwidth=mots[ivs->motindex].motwidth;
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
         int motwidth=mots[ivs->motindex].motwidth;
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
   outf << "<pre>" << endl;
   while (start<sizeseq){
      for (unsigned int spe=0;spe<1;spe++){
         if(seq.species[spe]){
            for (unsigned int i=start;i<stop;i++){
               if (vvstate[spe][i]==0){
                  outf << seq.seqsrealigned[spe][i];
               } else if (vvstate[spe][i]==1){
                  outf << "<span class=\"mot" << vvcol[spe][i] << "\">";
                  outf << seq.seqsrealigned[spe][i];
                  outf << "</span>";

               } else if (vvstate[spe][i]==2){
                  outf << "\\textit{\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                  outf << seq.seqsrealigned[spe][i];
                  outf << "}}";
               }
               if ( i>0 && i%10 == 0 ) {
                  outf << "\t";
               }
            }
            outf << "\\\\\n";
            int pstate=0;
            int ppos=start;
            for (unsigned int i=start;i<stop;i++){
               if (vvstate[spe][i]>0  && vvstate[spe][i]!=pstate){
                  int deca;
                  deca=i-ppos;
                  outf << "\\hspace*{" << deca << "\\charwidth}";
                  if (vvstate[spe][i]==2) outf << "\\textit{";
                  outf << "\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                  outf << setprecision(1) << fixed <<  scores[spe][i] ;
                  outf << "}";
                  if (vvstate[spe][i]==2) outf << "}";
                  // to count the number of digits:
                  stringstream score;
                  score << setprecision(1) << fixed <<  scores[spe][i] ;
                  for (int ii=1;ii<=score.str().size();ii++){
                     if ( i>0 && i%10 == 0 ) {
                        outf << "\\hspace*{1\\charwidth}";
                     }
                     i++;
                  }
                  ppos=i;
               }
               pstate=vvstate[spe][i];
               if ( i>0 && i%10 == 0 ) {
                  outf << "\\hspace*{1\\charwidth}";
               }
            }
            outf << "\\\\\n";
         }
      }
      start=stop;
      stop=min(stop+60,sizeseq);
   }
   outf << "</pre>" << endl;
}

// display TFBS on aligned sequences
   void
dispseqwmotswgaps (Sequence & seq, vmot & mots, ofstream & outf)
{
   //HEADER
   string name=numtospecies(0);
   string texname=seq.name;
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
      int motwidth=mots[ivs->motindex].motwidth;
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
         int motwidth=mots[ivs->motindex].motwidth;
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
               if ( i>0 && i%10 == 0 ) {
                  outf << "\t";
               }
            }
            outf << "\\\\\n";
            int pstate=0;
            int ppos=start;
            outf << "\\hspace*{" << 7 << "\\charwidth}";// species name + space = 7 characters
            for (unsigned int i=start;i<stop;i++){

               if (vvstate[spe][i]>0  && vvstate[spe][i]!=pstate){
                  int deca;
                  deca=i-ppos;
                  outf << "\\hspace*{" << deca << "\\charwidth}";
                  if (vvstate[spe][i]==2) outf << "\\textit{";
                  outf << "\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                  outf << setprecision(1) << fixed <<  scores[spe][i] ;
                  outf << "}";
                  if (vvstate[spe][i]==2) outf << "}";
                  // to count the number of digits:
                  stringstream score;
                  score << setprecision(1) << fixed <<  scores[spe][i] ;
                  for (int ii=1;ii<=score.str().size();ii++){
                     if ( i>0 && i%10 == 0 ) {
                        outf << "\\hspace*{1\\charwidth}";
                     }
                     i++;
                  }
                  ppos=i;
               }
               pstate=vvstate[spe][i];
               if ( i>0 && i%10 == 0 ) {
                  outf << "\\hspace*{1\\charwidth}";
               }
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

// general tex header
// TODO: more portable for international languages?
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
      "\\noindent\n"<<
      "\\newlength{\\charwidth}"<<
      "\\settowidth{\\charwidth}{\\texttt{A}}";
}

// general html header
// TODO: more portable for international languages?
   void
disphtmlinit(ofstream & outf)
{
   outf << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\" ";
   outf << "\"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\">" << endl;
   outf << endl;
   outf << "<html xmlns=\"http://www.w3.org/1999/xhtml\">" << endl;
   outf << endl;
   outf << "<head>" << endl;
   outf << "<title>Reference sequence annotation</title>" << endl;
   outf << "</head>" << endl;
   outf << endl;
   outf << "<body>" << endl;
}

   void
disptexclose(ofstream & outf)
{
   outf << "\\end{document}";
}

   void
disphtmlclose(ofstream & outf)
{
   outf << "</body>" << endl;
   outf << endl;
   outf << "</html>" << endl;
}

   void
disptex(vseq & seqs,vmot & mots)
{
   string filename;
   filename ="display/results.tex";
   ofstream outf(filename.c_str());
   if (outf.fail()){
      cerr << "Cannot open file for tex recording: " << strerror(errno) << endl;
      exit(-1);
   }

   disptexinit(outf);
   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){ 
      cout << "Scanning " << ivs->name << endl;
      dispseqwmots(*ivs,mots,outf);
   }
   disptexclose(outf);
   outf.close();
}

   void
disphtml(vseq & seqs,vmot & mots)
{
   string filename;
   filename ="display/results.html";
   ofstream outf(filename.c_str());
   if (outf.fail()){
      cerr << "Cannot open file for html recording: " << strerror(errno) << endl;
      exit(-1);
   }

   disphtmlinit(outf);
   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){ 
      cout << "Scanning " << ivs->name << endl;
      dispseqwmots_html(*ivs,mots,outf);
   }
   disphtmlclose(outf);
   outf.close();
}

   void
disptexwgaps(vseq & align,vmot & mots)
{
   string folder("display/");
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){ 
      stringstream file;
      file << folder;
      file << ivs->name;
      file << ".tex";
      ofstream outf(file.str().c_str());
      disptexinit(outf);
      cout << "Scanning " << (*ivs).name << endl;
      dispseqwmotswgaps(*ivs,mots,outf);
      disptexclose(outf);
      outf.close();
   }
}

   void
display_args_init()
{
   if (!strcmp(display_args.species_arg,"droso")){
      species="droso";
      nbspecies=12;
      conca=0.3; 
   } else if (!strcmp(display_args.species_arg,"eutherian")){
      species="eutherian";
      nbspecies=12;
      conca=0.263; 
   }
   concc=0.5-conca;
   conct=conca;
   concg=concc;

   // *** It would be nice to set the threshold by bp, in bits.
   scorethr2=display_args.threshold_arg;
   scorethr=scorethr2-1.0;
   scorethrcons=scorethr2-1.0;

   nbmots_for_score=display_args.nbmots_arg;

}

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

   display_args_init();

   cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;

   cout << "Loading alignments " << endl;

   vseq align;
   align=loadseqs(display_args.align_arg);
   cout << "Nb sequences to scan: " << align.size() << endl;

   cout << "Loading Motifs" << endl;

   vmot mots;
   loadmots(display_args.motifs_arg,mots); 
   if (nbmots_for_score<mots.size()) mots.erase(mots.begin()+nbmots_for_score,mots.end());

   for (ivmot iv=mots.begin();iv!=mots.end();iv++){

      // use same score on all species for detection
      iv->motscorethrcons=iv->motscorethr2;

      // avoid problems with _ and # characters for latex
      texify(iv->name);
   }
   cout << "Loaded " << mots.size() << " motifs." << endl;

   int mkdir_res = mkdir("display",S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH); 
   if (mkdir && errno != EEXIST){
      cerr << "Cannot create display directory: " << strerror(errno) << endl;
      exit(-1);
   }

   cout << "Scanning sequences for instances..." << endl;
   scanseqsforinstancesnmask(align,mots);

   cout << "Defining conserved instances..." << endl;
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
      ivs->instances2instancescons();
      //cout << ivs->name << "\n" << ivs->instancescons;
   }

   if (display_args.tex_ref_given){
      cout << "Creating tex file for reference species... " << endl;
      disptex(align,mots);
   }
   else if (display_args.html_ref_given){
      cout << "Creating html file for reference species... " << endl;
      disphtml(align,mots);
   }
   else if (display_args.tex_align_given){
      cout << "Creating fasta/tex files... " << endl;
      disptexwgaps(align,mots);
   }
   else if (display_args.svg_given){
      scanseqsforsvg(align,mots);
   }
   else{
      cout << "No mode was given. Exiting." << endl;
   }

   return 1;
}		/* -----  end of function extract  ----- */
