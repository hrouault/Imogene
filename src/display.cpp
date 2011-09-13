#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
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
      if (display_args.svg_given){
         nm=(*im).nbmatchnmaskforsvg(seq,moti);
      } else {
         //nm=(*im).nbmatchnmask(seq,moti);
         nm=(*im).nbmatchwomask(seq,moti);
      }
      nbmot=nm;

//      if (args_info.weightmots_given){
//         nmcorr+=nm*log((*im).lambdatrain/(*im).lambda);
//      }
//      else{
//         nmcorr+=nm;
//      }
//      nbcorr.push_back(nmcorr);
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
         string texname=seq.name;
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
         if (mots.size()>1){
            for (ivmot ivm=mots.begin()+1;ivm!=mots.end();ivm++){
               if (ivm->bsinit.size()<minwidth) minwidth=ivm->bsinit.size();
            }
         }

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
dispseqwmotswgaps (Sequence & seq, vmot & mots, ofstream & outf)
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
disptex(vseq & seqs,vmot & mots)
{

   double scoreinit=scorethr2;
   string filename("display/");
   filename+="results.tex";
   ofstream outf(filename.c_str());
   disptexinit(outf);
   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){ 
      cout << "Scanning " << ivs->name << endl;
      dispseqwmots(*ivs,mots,outf,scoreinit);
   }
   disptexclose(outf);

   outf.close();
}

   void
disptexwgaps(vseq & align,vmot & mots)
{

   cout << "Scanning sequences for instances..." << endl;
   scanseqsforinstancesnmask(align,mots);

   cout << "Defining conserved instances..." << endl;
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
      ivs->instances2instancescons();
      //cout << ivs->name << "\n" << ivs->instancescons;
   }

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


svg::svg()
{
   xsize=800;
   ysize=600;
   xoffset=140;
   yoffset=65;
   pos=0;
};

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
         svgfile << "<text transform=\"matrix(1 0 0 1 " << xmot-2 << " " << yline-8 << ")\" font-size=\"8\">" << fixed << setprecision(1) << (*ivi).score << "</text>\n";
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

      scanseq(*is,mots);

      //we set the size for the svg file
      xsize=s.xoffset+(int)(0.4*(*is).imaps[0].size());
      if (xsize>s.xsize) s.xsize=xsize;
      s.xsize+=s.xoffset;
      ysize=s.yoffset+340; 
      s.ysize=ysize+s.yoffset;

      svginit(svgfile,s);
      svgdisplay(svgfile,*is,s);
      svgclose(svgfile);

      //cout << s.xsize << " " << s.ysize << endl; 
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
   cout << "Loaded " << mots.size() << " motifs." << endl;
   for (ivmot iv=mots.begin();iv!=mots.end();iv++){
      iv->motscorethrcons=iv->motscorethr2;
      texify(iv->name);
   }

   mkdir("display",S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH); 

   if (display_args.tex_ref_given){

      cout << "Creating fasta/tex files... " << endl;
      disptex(align,mots);
   
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
