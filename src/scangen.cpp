#include<cmath>
#include<iostream>
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

#include "const.hpp"
#include "vectortypes.hpp"
#include "random.hpp"
#include "sequence.hpp"
#include "motif.hpp"
#include "tree.hpp"
#include "montecarlo.hpp"
#include "scangen.hpp"
#include "cmdline.h"

gengetopt_args_info args_info;
vmot motsdef;
vginst potregs;
vvginst groupedinst;
vstring phenos;
vstring gbacks;
vchrom chromints;

unsigned int sizepos,sizeneg;
unsigned int cutoff_for_combination=3;

   void
printsvmheader(ofstream & svm)
{
   svm << "## Support Vector Machine training set (chip VS background)" << endl;
   svm << "# crm_id:name of the nearest gene" << endl;
   svm << "# chr, start, stop: genomic coordinates." << endl;
   svm << "# {ChIP, Background}: Does the CRM belong (1) or not (0) to the class {'ChIP', 'Background'}" << endl;
   svm << "# mot_1,...,mot_n: number of sites for each motif" << endl;
   svm << "crm_id" << "\t";
   svm << "chr" << "\t";
   svm << "start" << "\t";
   svm << "stop" << "\t";
   svm << "ChIP" << "\t";
   svm << "Background";
   for (ivmot ivm=motsdef.begin();ivm!=min(motsdef.end(),motsdef.begin()+nbmots_for_score);ivm++){
      svm << "\t" << ivm->name;
   }
   svm << "\n";
}

   void
printforsvm(vseq & vscore)
{
   ofstream svm;
   svm.open("score/svm.dat");
   printsvmheader(svm);
   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
      svm << ivs->name << "\t";
      svm << chromfromint(ivs->chrom) << "\t";
      svm << ivs->start << "\t";
      svm << ivs->stop << "\t";
      if (ivs->sign==1){
         svm << 1 << "\t";
         svm << 0;
      } else {
         svm << 0 << "\t";
         svm << 1;
      }
      for (ivint iv=ivs->motis.begin();iv!=ivs->motis.end();iv++){
         svm << "\t" << *iv;
      }
      svm << "\n";
   }
   svm.close();
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
printbackreg(vcoord & vcds,string folder)
{
   //ofstream outf;
   ofstream outf("randpeaks.dat");
   unsigned int numback(0);
   unsigned int i(0);
   unsigned int backsize;
   if (args_info.size_given) backsize=args_info.size_arg;
   else backsize=2000;
   while (numback<10000){
      if (i>=vcds.size()) i=0;
      int range;
      range=vcds[i].stop-backsize-vcds[i].start;
      //cout << vcds[i];
      if (range<0){
         i++;
         //     cout << "range < 0" << endl;
         continue;
      }
      unsigned int newstart;
      unsigned int newpos;
      newpos= gsl_rng_uniform_int (gslran,range+1);
      newstart=vcds[i].start+newpos;
      Coordinate cdtmp;
      cdtmp.start=newstart;
      cdtmp.stop=newstart+backsize-1;
      cdtmp.chrom=vcds[i].chrom;
      cdtmp.name=vcds[i].name;
      // FOR PEAKS
      // 
      outf << "peaks_" << numback+1 << "\t";
      outf << chromfromint(cdtmp.chrom) << "\t";
      outf << cdtmp.start << "\n";
      numback++;

      //FOR SEQS
      //
      //         Sequence & s=seq;
      //         for (int k=0;k<nbspecies;k++){
      //            if (s.species[k]){
      //               if (k==0) outf << ">" << numtospecies(k) << " " <<
      //                  "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
      //               else  outf << ">" << numtospecies(k) << endl;
      //               outf << s.seqsrealigned[k] << endl;
      //            }; 
      //         }
      //         outf.close();
      //         numback++;
      //      }
      //      Sequence seq=coordtoseq(cdtmp);
      //      int nbfr=0;
      //      bool bc=0;
      //      if (species==1){
      //         if (seq.species[5]) nbfr++; 
      //         if (seq.species[6] || seq.species[7]) nbfr++;
      //         if (seq.species[8]) nbfr++;
      //         if (seq.species[9] || seq.species[10] || seq.species[11]) nbfr++;
      //         if (nbfr>1) bc=1;
      //      }
      //      else if (species==2){
      //         if (seq.species[2] || seq.species[3] || seq.species[4] || seq.species[5]) nbfr++;
      //         if (seq.species[6] || seq.species[7]) nbfr++;
      //         if (seq.species[8]) nbfr++;
      //         if (seq.species[9]) nbfr++;
      //         if (nbfr>1) bc=1;
      //      }
      //      if (seq.species[0] && bc && seq.nbtb>backsize/2){
      //         cout << "->" << numback << ". " <<cdtmp;
      //         stringstream os;
      //         os << numback;
      //         os >> seq.name;
      //         string filename=folder;
      //         filename+=seq.name;
      //         filename.append(".fa");
      //         outf.open(filename.c_str());
      //         Sequence & s=seq;
      //         for (int k=0;k<nbspecies;k++){
      //            if (s.species[k]){
      //               if (k==0) outf << ">" << numtospecies(k) << " " <<
      //                  "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
      //               else  outf << ">" << numtospecies(k) << endl;
      //               outf << s.seqsrealigned[k] << endl;
      //            }; 
      //         }
      //         outf.close();
      //         numback++;
      //      }
      //      else {
      //         //cout << "sequence not conserved or too much masked" << endl;
      //      }
      cout.flush();
      i++;
}
outf.close();
}

void
printbackregwcoords(vcoord & vcds,string folder){

   //vcds is a vector of random shuffled intergenic regions
   //regints are our interest sequences
   //
   ofstream outf;
   unsigned int i(0);
   unsigned int numback(1);

   unsigned int repeat=10;
   if (args_info.numrepeat_given) repeat=args_info.numrepeat_arg;

   for (int j=1;j<=repeat;j++){
      for (ivseq ivs=regints.begin();ivs!=regints.end();ivs++){

         //size of our interest sequence
         unsigned int backsize=ivs->stop-ivs->start+1;
         unsigned int state=0;

         //if we find a corresponding background region, then state=1;
         while (state==0){

            //if we have searched all intergenic sequences we loop on the first
            if (i==vcds.size()-1) i=0;

            int range;
            range=vcds[i].stop-backsize-vcds[i].start+1;

            //intergenic size has to be superior to interest sequence
            if (range<1){
               i++;
               continue;
            }

            //randomly choose a region of a same size than interest in the intergenic region
            unsigned int newstart;
            unsigned int newpos;
            newpos= gsl_rng_uniform_int (gslran,range);
            newstart=vcds[i].start+newpos;
            Coordinate cdtmp;
            cdtmp.start=newstart;
            cdtmp.stop=newstart+backsize-1;
            cdtmp.chrom=vcds[i].chrom;
            cdtmp.name=vcds[i].name;

            //extract the alignment
            Sequence seq=coordtoseq(cdtmp);
            int nbfr=0;
            bool bc=0;

            if (species==1){
               if (seq.species[5]) nbfr++; 
               if (seq.species[6] || seq.species[7]) nbfr++;
               if (seq.species[8]) nbfr++;
               if (seq.species[9] || seq.species[10] || seq.species[11]) nbfr++;
               if (nbfr>1) bc=1;
            }
            else if (species==2){
               if (seq.species[2] || seq.species[3] || seq.species[4] || seq.species[5]) nbfr++;
               if (seq.species[6] || seq.species[7]) nbfr++;
               if (seq.species[8]) nbfr++;
               if (seq.species[9]) nbfr++;
               if (nbfr>1) bc=1;
            }


            //we want mus to be present, the sequence to be conserved, and at least half of bases to be unmasked
            if (seq.species[0] && bc && seq.nbtb>backsize/2){
               cout << "->" << numback << ". " <<cdtmp;
               cout.flush();
               stringstream os;
               os << numback;
               os >> seq.name;
               string filename=folder;
               filename+=seq.name;
               filename.append(".fa");
               outf.open(filename.c_str());
               Sequence & s=seq;
               for (int k=0;k<nbspecies;k++){
                  if (s.species[k]){
                     if (k==0) outf << ">" << numtospecies(k) << " " <<
                        "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
                     else  outf << ">" << numtospecies(k) << endl;
                     outf << s.seqsrealigned[k] << endl;
                  }; 
               }
               outf.close();
               numback++;
               state=1;
            }
            i++;
         }

      }
   }

   outf.close();
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
      if (args_info.disp_svg_given){
         nm=(*im).nbmatchnmaskforsvg(seq,moti);
      } else {
         //nm=(*im).nbmatchnmask(seq,moti);
         nm=(*im).nbmatchwomask(seq,moti);
      }
      nbmot=nm;

      if (args_info.weightmots_given){
         nmcorr+=nm*log((*im).lambdatrain/(*im).lambda);
      }
      else{
         nmcorr+=nm;
      }
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
scanseqs(vstring & regs)
{
   vcoord pcoords;
   if (args_info.masktrain_given){

      cout << "(Masking training set)" << endl;
      ifstream inf;
      inf.open(args_info.masktrain_arg);
      back_insert_iterator<vcoord> pdest(pcoords);
      copy(iiscoord(inf),iiscoord(),pdest);
      inf.close();

   }
   
   string pchrom("");
   for (ivstring is=regs.begin();is!=regs.end();is++){
      //cout << *is << endl;
      Sequence seq;
      ifstream fseq;
      fseq.open((*is).c_str());
      string fseqline;
      getline(fseq,fseqline);
      seq.name=fseqline;
      stringstream firstline(fseqline);
      string speciesname;
      firstline >> speciesname;
      string chrom;
      firstline >> chrom;
      if (chrom!=pchrom){
         cout << chrom << endl;
         pchrom=chrom;
      }
      seq.chrom=intfromchrom(chrom.substr(3));
      if (seq.chrom==-1) continue;
      //cout << "chromosome \"" << seq.chrom << "\"\n";
      firstline >> seq.start;
      //cout << "start \"" << seq.start << "\"\n";
      firstline >> seq.stop;
      seq.finame=*is;
      string dum;
      seq.species=vint(nbspecies,0);
      vint dumi;
      seq.iseqs=vvint(nbspecies,dumi);
      seq.imaps=vvint(nbspecies,dumi);
      seq.imapsinv=vvint(nbspecies,dumi);
      while (!fseq.eof()){
         //                cout << fseqline << endl;
         int seqnum=speciestonum(fseqline.substr(1,6));
         getline(fseq,fseqline);
         seq.imaps[seqnum]=alignedtomap(fseqline);
         seq.imapsinv[seqnum]=alignedtorevmap(fseqline);
         string seqwogap=remgaps(fseqline);
         if (seqwogap.size()>30){
            seq.species[seqnum]=1;
         }
         seq.iseqs[seqnum]=stringtoint(seqwogap);
         getline(fseq,fseqline);
      }
      seq.nbN=compN(seq.iseqs[0]);
      seq.nbtb=seq.iseqs[0].size()-seq.nbN;
      
      for (ivcoord ivc=pcoords.begin();ivc!=pcoords.end();ivc++){
         if (ivc->chrom==seq.chrom){ 
            if ((ivc->start>=seq.start && ivc->start<=seq.stop) ||
                  (ivc->stop>=seq.start && ivc->stop <=seq.stop)){
                  int maskstart=max(ivc->start,seq.start);
                  int maskstop=min(ivc->stop,seq.stop);
      //            cout << ivc->name << endl;
      //            cout << vinttostring(seq.iseqs[0]) << endl;
                  for (unsigned int base=maskstart-seq.start;base<=maskstop-seq.start;base++){
                     seq.iseqs[0][base]=4;
                  }
      //            cout << vinttostring(seq.iseqs[0]) << endl;
            }
         }
      }

      //check that ref species contains at least 30bp
      if (seq.species[0]){
         scanseq(seq,motsdef);
      }

      fseq.close();
   }
   unsigned int i=0;

   cout << "Finding Instances" << endl;
   for (ivmot im=motsdef.begin();im!=motsdef.begin()+nbmots_for_score;im++){
      for (unsigned int j=0;j<nbchrom;j++){
         //           cout << "Instances size : " << (*im).instances.size() << "\n";
         sort((*im).instances[j].begin(),(*im).instances[j].end());
         for (ivinst iins=(*im).instances[j].begin();iins!=(*im).instances[j].end();iins++){
            Instance & curinst=*iins;
            int start=0;
            if ((curinst.coord-scanwidth)/scanstep+2>0) start=(curinst.coord-scanwidth)/scanstep+1;
            for (unsigned int k=start;k<curinst.coord/scanstep+1;k++){
               groupedinst[curinst.chrom][k].nbmots[i]++;
               groupedinst[curinst.chrom][k].totmots++;
               groupedinst[curinst.chrom][k].instances.push_back(curinst);
            }
         }
      }
      i++;
   }

   for (unsigned int j=0;j<nbchrom;j++){
      for (ivginst ivg=groupedinst[j].begin();ivg!=groupedinst[j].end();ivg++){
         if ((*ivg).totmots==0){
            (*ivg).discarded=1;
         }
      }
   }
   for (unsigned int j=0;j<nbchrom;j++){
      for (ivginst ivg=groupedinst[j].begin();ivg!=groupedinst[j].end();ivg++){
         if (!(*ivg).discarded){
            potregs.push_back(*ivg);
         }
      }
   }

   cout << "output results" << endl;
   if (args_info.all_mots_only_given){
         outputresults(nbmots_for_score);
   }
   else{
      for (unsigned int i=1;i<nbmots_for_score+1;i++){
         cout << "Motif " << i << endl;
         outputresults(i);
      }
   }

}

   void
scanseqs(ifstream & list,vcoord & coords)
{
   
   vstring regs;
   back_insert_iterator<vstring> dest(regs);
   copy(iisstring(list),iisstring(),dest);

   string pchrom("");
   for (ivstring is=regs.begin();is!=regs.end();is++){
      Sequence seq;
      ifstream fseq;
      fseq.open((*is).c_str());
      string fseqline;
      getline(fseq,fseqline);
      seq.name=fseqline;
      stringstream firstline(fseqline);
      string speciesname;
      firstline >> speciesname;
      string chrom;
      firstline >> chrom;
      if (chrom!=pchrom){
         cout << chrom << endl;
         pchrom=chrom;
      }
      seq.chrom=intfromchrom(chrom.substr(3));
      if (seq.chrom==-1) continue;
      firstline >> seq.start;
      firstline >> seq.stop;
      seq.finame=*is;
      string dum;
      seq.species=vint(nbspecies,0);
      vint dumi;
      seq.iseqs=vvint(nbspecies,dumi);
      seq.imaps=vvint(nbspecies,dumi);
      seq.imapsinv=vvint(nbspecies,dumi);
      while (!fseq.eof()){
         //                cout << fseqline << endl;
         int seqnum=speciestonum(fseqline.substr(1,6));
         getline(fseq,fseqline);
         seq.imaps[seqnum]=alignedtomap(fseqline);
         seq.imapsinv[seqnum]=alignedtorevmap(fseqline);
         string seqwogap=remgaps(fseqline);
         if (seqwogap.size()>30){
            seq.species[seqnum]=1;
         }
         seq.iseqs[seqnum]=stringtoint(seqwogap);
         getline(fseq,fseqline);
      }
      seq.nbN=compN(seq.iseqs[0]);
      seq.nbtb=seq.iseqs[0].size()-seq.nbN;
      
      scanseq(seq,motsdef);

      fseq.close();
   }
   unsigned int i=0;

   cout << "Finding Instances" << endl;
   for (ivmot im=motsdef.begin();im!=motsdef.begin()+nbmots_for_score;im++){
      for (unsigned int j=0;j<nbchrom;j++){
         //           cout << "Instances size : " << (*im).instances.size() << "\n";
         sort((*im).instances[j].begin(),(*im).instances[j].end());
         for (ivinst iins=(*im).instances[j].begin();iins!=(*im).instances[j].end();iins++){
            Instance & curinst=*iins;
            int start=0;
            if ((curinst.coord-scanwidth)/scanstep+2>0) start=(curinst.coord-scanwidth)/scanstep+1;
            for (unsigned int k=start;k<curinst.coord/scanstep+1;k++){
               groupedinst[curinst.chrom][k].nbmots[i]++;
               groupedinst[curinst.chrom][k].totmots++;
               groupedinst[curinst.chrom][k].instances.push_back(curinst);
            }
         }
      }
      i++;
   }

   cout << "Discarding empty Instances" << endl;
   for (unsigned int j=0;j<nbchrom;j++){
      for (ivginst ivg=groupedinst[j].begin();ivg!=groupedinst[j].end();ivg++){
         if ((*ivg).totmots==0){
            (*ivg).discarded=1;
         }
      }
   }
   for (unsigned int j=0;j<nbchrom;j++){
      for (ivginst ivg=groupedinst[j].begin();ivg!=groupedinst[j].end();ivg++){
         if (!(*ivg).discarded){
            potregs.push_back(*ivg);
         }
      }
   }

   cout << "output results" << endl;
   if (args_info.all_mots_only_given){
         outputresults(nbmots_for_score);
   }
   else{
      for (unsigned int i=1;i<nbmots_for_score+1;i++){
         cout << "Motif " << i << endl;
         outputresults(i);
      }
   }
}
   
   void
scanmots()
{
   ifstream potregs;
   if (species==1) potregs.open("/home/santolin/these/files/droso/align/all/align-files.dat");
   else if (species==2) potregs.open("/home/santolin/these/files/mus/epo/align-files.dat");
   //else if (species==2) potregs.open("/home/santolin/these/files/transfac/matrices/align-test.dat");
   vstring regs;
   back_insert_iterator<vstring> dest(regs);
   copy(iisstring(potregs),iisstring(),dest);
   potregs.close();

   string pchrom("");
   system("if ! test -d scanmots;then mkdir scanmots;fi;");      
   unsigned int totlen(0),totlentb(0);
   //random_shuffle(regs.begin(),regs.end());
   for (ivstring is=regs.begin();is!=regs.end();is++){
      Sequence seq;
      ifstream fseq;
      fseq.open((*is).c_str());
      string fseqline;
      getline(fseq,fseqline);
      seq.name=fseqline;
//cout << seq.name << endl;
      stringstream firstline(fseqline);
      string speciesname;
      firstline >> speciesname;
      string chrom;
      firstline >> chrom;
      if (chrom!=pchrom){
         cout << chrom << endl;
         pchrom=chrom;
      }
      seq.chrom=intfromchrom(chrom.substr(3));
      if (seq.chrom==-1) continue;
      firstline >> seq.start;
      firstline >> seq.stop;
      seq.finame=*is;
      string dum;
      seq.species=vint(nbspecies,0);
      vint dumi;
      seq.iseqs=vvint(nbspecies,dumi);
      seq.imaps=vvint(nbspecies,dumi);
      seq.imapsinv=vvint(nbspecies,dumi);
      while (!fseq.eof()){
//cout << fseqline << endl;
         int seqnum=speciestonum(fseqline.substr(1,6));
         getline(fseq,fseqline);
         seq.imaps[seqnum]=alignedtomap(fseqline);
         seq.imapsinv[seqnum]=alignedtorevmap(fseqline);
         string seqwogap=remgaps(fseqline);
         if (seqwogap.size()>30){
            seq.species[seqnum]=1;
         }
         seq.iseqs[seqnum]=stringtoint(seqwogap);
         getline(fseq,fseqline);
      }
      seq.nbN=compN(seq.iseqs[0]);
      seq.nbtb=seq.iseqs[0].size()-seq.nbN;

      scanseqforconsinstances(seq,motsdef);

      fseq.close();
      
      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
         //sort(ivm->refinstances_short.begin(),ivm->refinstances_short.end());
         ostringstream outfs;
         outfs << "scanmots/" << ivm->name << "_" << ivm->motscorethr2 << ".dat";
         ofstream outf;
         if (totlen==0) outf.open(outfs.str().c_str());
         else outf.open(outfs.str().c_str(),ios::app);

         for (ivinst ivi=ivm->refinstances_short.begin();ivi!=ivm->refinstances_short.end();ivi++){
            outf << ivi->site << "\t";
            outf << chromfromint(ivi->chrom) << "\t";
            outf << ivi->coord << "\t";
            outf << ivi->sens << "\t";
            outf << ivi->score << endl;
         }
         outf.close();
         ivm->refinstances_short.clear();
      }

      totlen+=seq.nbtb+seq.nbN;
      totlentb+=seq.nbtb;
   }

   cout << "Total length of the alignement: " << totlen << " bp, including " << totlentb << " unmasked bp" << endl;

   return;
}

bool operator<(const Sequence & seqscr1,const Sequence & seqscr2)
{
   if (seqscr1.score == seqscr2.score)
   {
      return seqscr1.sign > seqscr2.sign; //for identical scores, we want the negatives first (for ROC consistency)
   }
   else
   {
      return seqscr1.score > seqscr2.score; 
   };
}

bool operator<(const Motif & mot1,const Motif & mot2)
{
   if (mot1.optauc == mot2.optauc)
   {
      return mot1.index < mot2.index;//the lower index first 
   }
   else
   {
      return mot1.optauc > mot2.optauc;//the best auc first
   };
}

   void
calcscore(Motif & im,Sequence & seqscr)//(vmot & mots)
{
   unsigned int ncons=im.nbmatchcons(seqscr);//seqscr.motis[im->index];
   //   unsigned int ncons=im.nbmatchmat(seqscr.seq);//seqscr.motis[im->index];
   //   unsigned int ncons=seqscr.motis[im.index];
   double lseq=seqscr.nbtb;
   double score=0;
   //cout << ncons << " ";
   //score= ncons;//*log(im.lambdatrain/im.lambda) + lseq*(im.lambda-im.lambdatrain);
   double lnj=0;
   for (int i=0;i<ncons;i++){
      lnj+=log(i+1);
   }
   seqscr.motis[im.index]=ncons;
   score=ncons;
   //score=ncons*log(im.lambdatrain/im.lambda);
   //score= lnj-ncons*log(im.lambda*lseq) + lseq*im.lambda;
   //   score= (double)ncons/lseq;//*log(im.lambdatrain/im.lambda) + lseq*(im.lambda-im.lambdatrain);
   seqscr.score=score;
   return;
   //cout << seqscr.score << endl;
}

   void
calcscore(vmot & vm,Sequence & seqscr)//(vmot & mots)
{

   Sequence seqtmp=seqscr;
   double lseq=seqscr.nbtb;
   double score=0;
   for (ivmot im=vm.begin();im!=min(vm.begin()+nbmots_for_score,vm.end());im++){
      width=im->bsinit.size();
      unsigned int ncons;
      if (args_info.maskforscore_given){
         ncons=im->nbmatchconsnmask(seqtmp);//seqscr.motis[im->index];
      } else{
         ncons=im->nbmatchcons(seqtmp);//seqscr.motis[im->index];
      }
      //   score+=ncons;
      //      double lnj=0;
      //      for (int i=0;i<ncons;i++){
      //         lnj+=log(i+1);
      //      }
      //   //   score+=lnj-ncons*log(im->lambda*lseq) + lseq*im->lambda;
      //      if (args_info.byaffinity_given){
      //         score+=im->scorematchcons(seqtmp);
      //      } else {
      //         score+=ncons;
      //      }
      if (args_info.weightmots_given){

         score+= ncons*log(im->lambdatrain/im->lambda);// + lseq*(im.lambda-im.lambdatrain);
      }
      else {
         score+= ncons;//*log(im->lambdatrain/im->lambda);// + lseq*(im.lambda-im.lambdatrain);
      }
      seqscr.motis[im->index]=ncons;
   }
   seqscr.score=score;
   return;
   //cout << seqscr.score << endl;
}

   void
calcscore(vmot & vm,vseq & vs)//(vmot & mots)
{
   for (ivseq ivs=vs.begin();ivs!=vs.end();ivs++){
      calcscore(vm,*ivs);
   }
   return;
   //cout << seqscr.score << endl;
}

// allncons is a vint where each int is a number of conserved motifs
   double
calcscore(vmot & vm,vint & allncons,int lseq)
{

   double score=0;
   for (ivmot im=vm.begin();im!=min(vm.begin()+nbmots_for_score,vm.end());im++){
      unsigned int ncons=allncons[im->index];
      double lnj=0;
      for (int i=0;i<ncons;i++){
         lnj+=log(i+1);
      }
      //score+=lnj-ncons*log(im->lambda*lseq) + lseq*im->lambda;
      //    score+=allncons[im->index];
      score+= allncons[im->index];//*log(im->lambdatrain/im->lambda);// + lseq*(im.lambda-im.lambdatrain);
   }
   return score;
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
         dro = seq.speciesname[i];
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
scanseqsforsvg(vseq & align)
{

   system("if ! test -d display;then mkdir display;fi;");      
   system("if ! test -d display/svg;then mkdir display/svg;fi;");      

   for (ivseq is=align.begin();is!=align.end();is++){

      int xsize(0);
      int ysize(0);
      svg s;
      string filename("display/svg/");
      filename+=(*is).name;
      filename+=".svg";
      ofstream svgfile(filename.c_str());

      scanseq(*is,motsdef);

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
outputresults(unsigned int nbmots_score)
{
   cout << "Shuffling" << endl;
   random_shuffle(potregs.begin(),potregs.end());
   for (ivginst ivg=potregs.begin();ivg!=potregs.end();ivg++){
      if (args_info.weightmots_given){
         (*ivg).compscoreweight(motsdef,nbmots_score);
      }
      else{
         (*ivg).compscore(motsdef,nbmots_score);
      }
   }
   cout << "Sorting" << endl;
   sort(potregs.begin(),potregs.end());
   cout << "Discarding" << endl;
   vstring gnames;
   vginst potregs_def;

   // Keep best scoring enhancer per gene
   // OR
   // Keep all enhancers per gene, removing overlapping low score ones
   if (args_info.discard_on_gene_names_given){
      cout << "discard on gene" << endl;
      int test=0;
      for (ivginst ivg=potregs.begin();ivg!=potregs.end();ivg++){
         test=0;
         for (ivstring ivs=gnames.begin();ivs!=gnames.end();ivs++){
            if ((*ivg).besttss.gene==*ivs){
               (*ivg).discarded=1;
               test=1;
               break;
            }
         }
         if (test==0){
            gnames.push_back((*ivg).besttss.gene);
         }
      }
   } 
   else {
      for (ivginst ivg=potregs.begin();ivg!=potregs.end();ivg++){
         int test=0;
         for (ivginst ivg2=potregs_def.begin();ivg2!=potregs_def.end();ivg2++){
            if ((*ivg).distance(*ivg2)<scanwidth){
               test=1;
               break;
            }
         }
         if (test==0){
            potregs_def.push_back(*ivg);
         }
      }
   }
   cout << "Displaying" << endl;
   stringstream filename;
   filename << "result" << nbmots_score << ".dat";
   ofstream res(filename.str().c_str());
   for (ivginst ivg=potregs_def.begin();ivg!=potregs_def.end();ivg++){
      if (!(*ivg).discarded){
         res << (*ivg).score << " " << chromfromint((*ivg).chrom) << ":" << (*ivg).start << ".." << (*ivg).stop << " ";
         res << (*ivg).besttss.gene << " ";
         for (ivTSS ivt=(*ivg).TSSs.begin();ivt!=(*ivg).TSSs.end();ivt++){
            res << (*ivt).gene << ";";
         }
         res << " ";
         for (ivint ivi=(*ivg).nbmots.begin();ivi!=(*ivg).nbmots.end();ivi++){
            res << *ivi << " ";
         }
         res << "\n";
         //         for (ivinst ivi=(*ivg).instances.begin();ivi!=(*ivg).instances.end();ivi++){
         //            cout << "\t" << chromfromint((*ivi).chrom) << ":" << (*ivi).coord << "\n";
         //         }
      }
   }
   res << "\n";
   res.close();

   //	stringstream filenameseqs;
   //	filenameseqs << "seqs-" << nbmots_score << ".dat";
   //	ofstream seqs(filenameseqs.str().c_str());
   //	unsigned int count=0;
   //	for (ivginst ivg=potregs_def.begin();ivg!=potregs_def.end();ivg++){
   //		if (!(*ivg).discarded){
   //			unsigned int mini=scanwidth;
   //			unsigned int max=0;
   //			seqs << "chrom : " << (*ivg).chrom << "\n";
   //			seqs << "start : " << (*ivg).start << "\n";
   //			seqs << "stop : " << (*ivg).stop << "\n";
   //			for (ivinst ivi=(*ivg).instances.begin();ivi!=(*ivg).instances.end();ivi++){
   //				if ((*ivi).coord<mini) mini=(*ivi).coord;
   //				if ((*ivi).coord>max) max=(*ivi).coord;
   //			}
   //			seqs << "mini : " << mini << "\n";
   //			seqs << "max : " << max << "\n";
   //			unsigned int start=(mini+max)/2-scanwidth/2;
   //			unsigned int stop=(mini+max)/2+scanwidth/2;
   //			for (ivinst ivi=(*ivg).instances.begin();ivi!=(*ivg).instances.end();ivi++){
   //				seqs << "motif " << (*ivi).motindex << " at position " << (*ivi).coord-start << "\n";;
   //			}
   //			seqs << ">seq_" << count << " annotated " << (*ivg).besttss.gene;
   //			seqs << " ";
   //			seqs << "\n";
   //			seqs << chromints[(*ivg).chrom].seq.substr(start,stop-start) << "\n";
   //			count++;
   //			if (count>100) break;
   //		}
   //	}
   //	seqs << "\n";
   //	seqs.close();

   if (args_info.phenotype_given){
      stringstream filehist;
      filehist << "hist" << nbmots_score << ".dat";
      ofstream histo(filehist.str().c_str());
      displayhist(potregs_def,histo);
      histo.close();
      if (args_info.print_histo_sets_given){
         stringstream filehist_back;
         filehist_back << "hist-back" << nbmots_score << ".dat";
         ofstream histo_back(filehist_back.str().c_str());
         displayhist_set(potregs_def,gbacks,histo_back);
         histo_back.close();

         stringstream filehist_interest;
         filehist_interest << "hist-interest" << nbmots_score << ".dat";
         ofstream histo_interest(filehist_interest.str().c_str());
         displayhist_set(potregs_def,phenos,histo_interest);
         histo_interest.close();
      }
   }

}


   void
motiftomat(vint & seq,Motif & mot)
{
   vd dum(4,0);
   mot.matrice=vvd(width,dum);
   vvd::iterator imat=mot.matrice.begin();
   for (vint::iterator iseq=seq.begin();iseq!=seq.end();iseq++){
      int base=*iseq;
      double a,t,c,g;
      a=t=alpha/(1+2*alpha+2*beta);
      c=g=beta/(1+2*alpha+2*beta);
      if (base==0){
         a=(1+alpha)/(1+2*alpha+2*beta);
      }
      else if (base==1){
         t=(1+alpha)/(1+2*alpha+2*beta);
      }
      else if (base==2){
         c=(1+beta)/(1+2*alpha+2*beta);
      }
      else if (base==3){
         g=(1+beta)/(1+2*alpha+2*beta);
      }
      (*imat)[0]=log(a/conca);
      (*imat)[1]=log(t/conct);
      (*imat)[2]=log(c/concc);
      (*imat)[3]=log(g/concg);
      imat++;
   }
}


vseq regtests;
vseq regints;

   void
extracttofastawfullname(string folder)
{
   ofstream outf;
   for (ivseq iv=regints.begin();iv!=regints.end();iv++){
      Sequence seq=*iv;
      if (seq.species[0] && seq.nbtb>0){ 
         stringstream file;
         file << folder;
         file << seq.name << "_";
         file << chromfromint(seq.chrom) << "_";
         file << seq.start << "_";
         file << seq.stop << ".fa";
         outf.open(file.str().c_str());
         Sequence & s=seq;
         for (int i=0;i<nbspecies;i++){
            if (s.species[i]){
               if (i==0) outf << ">" << numtospecies(i) << " " <<
                  "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
               else  outf << ">" << numtospecies(i) << endl;
               outf << s.seqsrealigned[i] << endl;
            }; 
         }
         outf.close();
      }
   }
}

   void
extracttofasta(string folder)
{
   ofstream outf;
   string pname("");
   int pnum(1);
   for (ivseq iv=regints.begin();iv!=regints.end();iv++){
      Sequence seq=*iv;
      if (seq.species[0] && seq.nbtb>0){ 
         stringstream file;
         file << folder;
         file << seq.name;
         if (seq.name==pname){ 
            file << "_";
            file << pnum;
            pnum++;
         }
         else {
            pnum=1;
            pname=seq.name;
         }
         file << ".fa";
         outf.open(file.str().c_str());
         Sequence & s=seq;
         for (int i=0;i<nbspecies;i++){
            if (s.species[i]){
               if (i==0) outf << ">" << numtospecies(i) << " " <<
                  "chr" << chromfromint(seq.chrom) << " " <<  seq.start << " " << seq.stop << endl;
               else  outf << ">" << numtospecies(i) << endl;
               outf << s.seqsrealigned[i] << endl;
            }; 
         }
         outf.close();
      }
   }
}


   double
distcv(vvd& mat1, vvd& mat2)
{
   double max(0);
   int col(0);
   for (ivvd iv=mat1.begin();iv!=mat1.end();iv++){
      int line(0);
      for (ivd il=(*iv).begin();il!=(*iv).end();il++){
         double diff;
         diff=fabs((*il)-mat2[col][line]);
         if (diff>max) max=diff;
         line++;
      }
      col++;
   }
   return max;
}

   void
matflank(Motif & mot,int numflank)
{
   vd dum(4,0);
   vvd flank_left(numflank,dum);
   vvd flank_right(numflank,dum);

   for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){
      vint vdum;
      for (int i=1;i<=numflank;i++){
         if ((*ivm).strand==1){
            ivint iv_left=(*ivm).matchespos[0]-i;
            if (iv_left>(*ivm).seq_start) flank_left[numflank-i][(*iv_left)]+=1;
            ivint iv_right=(*ivm).matchespos[0]+width+i;
            if (iv_right<(*ivm).seq_stop) flank_right[i-1][(*iv_right)]+=1;
         }
         else if ((*ivm).strand==-1){
            ivint iv_left=(*ivm).matchespos[0]+width+i;
            if (iv_left<(*ivm).seq_stop){
               int brc;
               if (*iv_left==0) brc=1;
               if (*iv_left==1) brc=0;
               if (*iv_left==2) brc=3;
               if (*iv_left==3) brc=2;
               flank_left[numflank-i][brc]+=1;
            }  
            ivint iv_right=(*ivm).matchespos[0]-i;
            if (iv_right>(*ivm).seq_start){ 
               int brc;
               if (*iv_right==0) brc=1;
               if (*iv_right==1) brc=0;
               if (*iv_right==2) brc=3;
               if (*iv_right==3) brc=2;
               flank_right[i-1][brc]+=1;
            }  
         }
      }

      for (ivint iv=(*ivm).matchespos[0];iv!=(*ivm).matchespos[0]+width;iv++){
         vdum.push_back(*iv);
      }
      // cout << flank_left << "\t" << vinttostring(vdum) << "\t" << flank_right << endl; 
   }

   alpha=conca;
   beta=concc;
   countfreq(flank_left);
   freqtolog(flank_left);
   countfreq(flank_right);
   freqtolog(flank_right);
   //   cout << flank_left <<  "\t" << flank_right << endl; 

   vvd finmat;
   for (ivvd iv=flank_left.begin();iv!=flank_left.end();iv++){
      finmat.push_back(*iv);
   }
   for (ivvd iv=mot.matprec.begin();iv!=mot.matprec.end();iv++){
      finmat.push_back(*iv);
   }
   for (ivvd iv=flank_right.begin();iv!=flank_right.end();iv++){
      finmat.push_back(*iv);
   }

   width+=2*numflank;
   mot.matprec=finmat;
   mot.matprecrevcomp=reversecomp(finmat);

   return;
}


   void
motaffin(vmot & mots,ofstream & streamfile)
{
   double icbinit=scorethr2/10.;
   for (ivmot ivm=mots.begin();ivm!=mots.end();ivm++){

      scorethr2=ivm->motscorethr2;
      compalpha();

      vvd mattot;
      string bstot;
      Motif currmot=(*ivm);
      width=currmot.bsinit.size();
      vint bs=stringtoint(currmot.bsinit);
      vvd pmat=currmot.matprec;

      //cout << "Motif # " << currmot.index << "..." << endl;

      double icb=icbinit;
      double max=0.01;
      int iter(1);
      while(max!=0){
         if (iter>1){
            scorethr2=icb*width;
            scorethr=(icb-0.1)*width;
            scorethrcons=scorethr;
            compalpha();
            currmot.matinit(scorethr2);
         }
         else{
            scorethr2=(icb-0.1)*width;
            scorethr=scorethr2;
            scorethrcons=scorethr;
            compalpha();
            currmot.matinit(scorethr);
         }
         //         if (currmot.nbmot<3){
         //            scorethr-=1;
         //            scorethrcons=scorethr;  
         //            iter=1;
         //            cout << endl;
         //            continue;
         //           // break;
         //         }
         if (currmot.nbmot<1) break;

         if (args_info.MCMCtest_given){
            cout << "TEST!" << endl;
            currmot.compprec_test();
            exit(1);
         }

         currmot.compprec();
         max=distcv(currmot.matprec,pmat);
         //cout << max << "\n";
         pmat=currmot.matprec;
         iter++;
         if (iter==2) break;
      }

      int numflank=5;
      //     matflank(currmot,numflank);

      currmot.matinit(scorethr2);
      currmot.pvaluecomp();
      currmot.display(streamfile);
      //cout << currmot.bsinit << " CV in " << iter << ". Nb mot: " << currmot.nbmot << " " << endl;
   }
   //cout << endl;
}

   void
seqanalysis(Sequence & currseq,ofstream & streamfile)
{
   unsigned int i=0;
//      for (int j=0;j<nbspecies;j++){
//      cout << currseq.iseqs[j] << endl;
//      }
   for (vint::iterator istr=currseq.iseqs[0].begin();istr!=currseq.iseqs[0].end()-width+1;istr++){
      //cout << "\r" << i+1 << "/" << currseq.iseqs[0].size()-width+1 ; 
      vint bs(istr,istr+width);
     //cout << i << " " << bs << endl;
      if (compN(bs)>0) continue;
      Motif currmot;
      currmot.bsinit=vinttostring(bs);
      currmot.seqinit=currseq.name;
      currmot.pos=i;
      motiftomat(bs,currmot);
      currmot.matricerevcomp=reversecomp(currmot.matrice);
      currmot.matprec=currmot.matrice;
      currmot.matprecrevcomp=currmot.matricerevcomp;
      vvd pmat=currmot.matprec;
      unsigned int nbconv(0);
      for (int nb=1;nb<=nbiter;nb++){
         double max=0.01;
         int iter(0);
         while(max>0){
            if (nb>2) currmot.matinit(scorethr2);
            else currmot.matinit(scorethr);
            if (currmot.nbmot<1) break;
            
            if (args_info.MCMCtest_given){
               cout << "TEST!" << endl;
               currmot.compprec_test();
               exit(1);
            }
            
            currmot.compprec();
            max=distcv(currmot.matprec,pmat);
            pmat=currmot.matprec;
            iter++;
            if (iter==20) break;
            nbconv++;
         }
         if (nb==1){
            currmot.corrprec();
            pmat=currmot.matprec;
            currmot.matprecrevcomp=reversecomp(currmot.matprec);
         }
      }
      //cout << currmot.nbmot << " " ; 
      //cout.flush();

      currmot.matinit(scorethr2);
      if (currmot.nbmot>2){
         currmot.pvaluecomp();
         currmot.display(streamfile);
         //		for (ivma ima=currmot.seqs.begin();ima!=currmot.seqs.end();ima++){
         //  cout << (*ima).alignseq[0] << endl;
         //		}
         //cout << currmot.matprec << endl;
      }
      i++;
   }
   cout << endl;
}
   
   void
refseqanalysis(Sequence & currseq,ofstream & streamfile)
{
   unsigned int i=0;
   double scoresave=scorethr2;
   //   for (int j=0;j<nbspecies;j++){
   //   cout << currseq.iseqs[j] << endl;
   //   }
   for (vint::iterator istr=currseq.iseqs[0].begin();istr!=currseq.iseqs[0].end()-width+1;istr++){
      //cout << "\r" << i+1 << "/" << currseq.iseqs[0].size()-width+1 ; 
      vint bs(istr,istr+width);
      //cout << i << " " << bs << endl;
      if (compN(bs)>0) continue;
      Motif currmot;
      currmot.bsinit=vinttostring(bs);
      currmot.seqinit=currseq.name;
      currmot.pos=i;
      motiftomat(bs,currmot);
      currmot.matricerevcomp=reversecomp(currmot.matrice);
      currmot.matprec=currmot.matrice;
      currmot.matprecrevcomp=currmot.matricerevcomp;
      vvd pmat=currmot.matprec;
      unsigned int nbconv(0);
      for (int nb=1;nb<=nbiter;nb++){
         double max=0.01;
         int iter(0);
         while(max>0){
            if (nb>2) scorethr2=scoresave;
            else scorethr2=scorethr;
            
            currmot.matinit(scorethr2);
            if (currmot.nbmot<1) break;
            
            currmot.comprefinstances(regints,0);
            currmot.comprefmot();
            max=distcv(currmot.matprec,pmat);
            pmat=currmot.matprec;
            //        cout << max << endl;
            iter++;
            if (iter==20) break;
            nbconv++;
         }
         if (nb==1){
            currmot.corrprec();
            pmat=currmot.matprec;
            currmot.matprecrevcomp=reversecomp(currmot.matprec);
         }
      }
      cout << currmot.nbmot << " " ; 
      cout.flush();

      currmot.comprefinstances(regints,0);
      currmot.matinit(scorethr2);
      if (currmot.nbmot>2){
         currmot.pvaluecomp();
         currmot.display(streamfile);
         //		for (ivma ima=currmot.seqs.begin();ima!=currmot.seqs.end();ima++){
         //  cout << (*ima).alignseq[0] << endl;
         //		}
         //cout << currmot.matprec << endl;
      }
      i++;
   }
   cout << endl;
}

vTSS TSSall;

   void
loadannots()
{
   ifstream annots;
   annots.open(args_info.phenotype_arg);

   back_insert_iterator<vstring> dest(phenos);
   copy(iisstring(annots),iisstring(),dest);

   annots.close();

   ifstream glist;
   //!!!! FIRST scangens were done with genelist-wosensory
   //   if (species==1) glist.open("/home/santolin/these/files/droso/plaza/phenos/neg-ovo-pheno.dat");
   //   else if (species==2) glist.open("/home/santolin/these/files/mus/affymetrix/e10/genelist-wo-pos-down-uniq.dat");
   if (args_info.phenoback_given){
      glist.open(args_info.phenoback_arg);
   }
   else {
      if (species==1) glist.open("/home/rouault/these/sequence/genomes/genelist.dat");
      else if (species==2) glist.open("/home/santolin/these/files/mus/biomart/genelist-protein-coding+miRNA.dat");
   }

   back_insert_iterator<vstring> destg(gbacks);
   copy(iisstring(glist),iisstring(),destg);

   glist.close();

   for (ivstring ivs=phenos.begin();ivs!=phenos.end();ivs++){
      for (ivstring ivs2=gbacks.begin();ivs2!=gbacks.end();ivs2++){
         if (*ivs2==*ivs){
            gbacks.erase(ivs2);
            break;
         }
      }
   }

   //   for (ivstring ivs=phenos.begin();ivs!=phenos.end();ivs++){
   //      cout << *ivs << endl;
   //   }
}

   vcoord
loadcoords()
{
   ifstream inf;
   inf.open(args_info.coord_file_arg);
   vcoord coords;
   back_insert_iterator<vcoord> dest(coords);
   copy(iiscoord(inf),iiscoord(),dest);
   return coords;
}


   void
loadchroms()
{
   ifstream chromsf;
   if (species==1){
      chromsf.open("/home/rouault/these/sequence/genomes/melano-only/dmel-all-chromosome-r4.3.fasta");
   } else if (species==2){
      chromsf.open("/home/santolin/these/files/mus/genome/fasta/all/all.chromosome_no_MT.fasta");
      //chromsf.open("/home/santolin/these/files/mus/genome/fasta/all/chr4.fa");
      //chromsf.open("/home/santolin/these/files/mus/genome/fasta/repeat_masked/all.chromosome_no_MT.fa");
   }
   back_insert_iterator<vchrom> dest(chromints);
   copy(iischrom(chromsf),iischrom(),dest);
   sort(chromints.begin(),chromints.end());

   for (int i=0;i<chromints.size();i++){
      cout << "chrom " << chromints[i].name << " ";
      cout << "of length " << chromints[i].seq.size() << " : ";
      cout << chromints[i].seq.substr(0,50) << " ... ";
      cout << chromints[i].seq.substr(chromints[i].seq.size()-50);
      cout << "\n";
   }


   chromsf.close();
}

   vint 
loadlengthchrom()
{

   vint lchr(nbchrom,0);

   if (species==1){

      lchr[0]=22410834;//chr2L
      lchr[1]=20769785;//chr2R
      lchr[2]=23774897;//chr3L
      lchr[3]=27908053;//chr3R
      lchr[4]=1284640;//chr4
      lchr[5]=22227390;//chrX


   } else if (species==2){

      lchr[0]=197195432;
      lchr[1]=181748087;
      lchr[2]=159599783;
      lchr[3]=155630120;
      lchr[4]=152537259;
      lchr[5]=149517037;
      lchr[6]=152524553;
      lchr[7]=131738871;
      lchr[8]=124076172;
      lchr[9]=129993255;
      lchr[10]=121843856;
      lchr[11]=121257530;
      lchr[12]=120284312;
      lchr[13]=125194864;
      lchr[14]=103494974;
      lchr[15]=98319150;
      lchr[16]=95272651;
      lchr[17]=90772031;
      lchr[18]=61342430;
      lchr[19]=166650296;//X
      lchr[20]=15902555;//Y
   }

   return lchr;

}

   void
initgroupedinst()
{
   cout << "Loading length chroms..." << endl;
   lengthchrom=loadlengthchrom();

   vginst dumvginst;
   groupedinst=vvginst(nbchrom,dumvginst);
   for (unsigned int i=0;i<nbchrom;i++){
      for (unsigned int j=0;j<lengthchrom[i]/scanstep+1;j++){
         groupedinst[i].push_back(GroupInstance(j*scanstep,j*scanstep+scanwidth,i));
      }
   }
   ifstream annots;
   if (species==1){
      annots.open("/home/rouault/these/sequence/genomes/regres-wellform-all.dat");
   } else if (species==2){
      annots.open("/home/santolin/these/files/mus/biomart/genes-n-strand-protein-coding+miRNA-no-MT-n-NT-TSS.dat");
   }
   importTSS(TSSall,annots);
   annots.close();
   // assign CRMs to genes in an annotextent region
   cout << "Assign CRMs to genes in annotextend..." << endl;
   for (ivTSS ivt=TSSall.begin();ivt!=TSSall.end();ivt++){
      int start=0;
      if (((*ivt).coord-annotextent-scanwidth/2)/scanstep+1>0) start=((*ivt).coord-annotextent-scanwidth/2)/scanstep;
      int stop=((*ivt).coord+annotextent-scanwidth/2)/scanstep+1;
      if (stop>lengthchrom[(*ivt).chrom]/scanstep+1) stop=lengthchrom[(*ivt).chrom]/scanstep+1;
      for (unsigned int i=start;i<stop;i++){
         groupedinst[(*ivt).chrom][i].TSSs.push_back(*ivt);
         //       cout << (*ivt).gene << endl;
      }
   }
   // assign nearest TSS to genes
   cout << "Assign CRMs to nearest gene..." << endl;
   for (unsigned int i=0;i<nbchrom;i++){
      for (ivginst ivg=groupedinst[i].begin();ivg!=groupedinst[i].end();ivg++){
         (*ivg).compbestannot();
      }
   }

   // attribute phenotype to CRMs
   cout << "Attribute phenotype to CRM..." << endl;
   for (unsigned int i=0;i<nbchrom;i++){
      for (ivginst ivg=groupedinst[i].begin();ivg!=groupedinst[i].end();ivg++){
         (*ivg).isdiscarded();
         if (!(*ivg).discarded){
            for (ivstring ivs=phenos.begin();ivs!=phenos.end();ivs++){
               if ((*ivg).besttss.gene==*ivs){
                  (*ivg).goodpheno=1;
                  break;
               } else {
                  (*ivg).goodpheno=0;
               }
            }
         } else {
            (*ivg).goodpheno=0;
         }
      }
   }

}

   void
getsites(ifstream & file, Sequence & bds)
{
   string dum;
   string dumtest;
   int i=0;
   cout << "seq " << endl;
   getline(file,dum);
   while(!file.eof()){
      bds.sitenames.push_back(dum.substr(1));
      getline(file,dum);
      if (i>0 && dum.length()!=dumtest.length()) {
         cout << dum << endl;
         cout << "Error in site length" << endl;
         exit(2);
      }
      dumtest=dum;
      i++;
      bds.seqs.push_back(dum);
      cout << dum << endl;
      bds.iseqs.push_back(stringtoint(dum));
      getline(file,dum);
   }
}

   void
disptex(vseq & seqs)
{
   system("if ! test -d display;then mkdir display;fi;");      
   //   system("if ! test -d display/tex;then mkdir display/tex;fi;");      

   double scoreinit=scorethr2;
   string filename("display/");
   filename+="results.tex";
   ofstream outf(filename.c_str());
   disptexinit(outf);
   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){ 
      cout << "Scanning " << (*ivs).name << endl;
      dispseqwmots(*ivs,motsdef,outf,scoreinit);
   }
   disptexclose(outf);

   outf.close();
}

   void
disptexwgaps(vseq & align)
{
   system("if ! test -d display;then mkdir display;fi;");      

   cout << "Scanning sequences for instances..." << endl;
   scanseqsforinstancesnmask(align,motsdef);

   cout << "Defining conserved instances..." << endl;
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
      ivs->instances2instancescons();
      //cout << ivs->name << "\n" << ivs->instancescons;
   }

   string folder("display/");
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){ 
      stringstream file;
      file << folder;
      file << ivs->name << "_";
      file << chromfromint(ivs->chrom) << "_";
      file << ivs->start << "_";
      file << ivs->stop;
      file << ".tex";
      ofstream outf(file.str().c_str());
      disptexinit(outf);
      cout << "Scanning " << (*ivs).name << endl;
      dispseqwmotswgaps(*ivs,outf);
      disptexclose(outf);
      outf.close();
   }

}

   Sequence
findbestseq(vmot & vm, Sequence & seq)
{

   Sequence bestseq=seq;

   //Find the best 1kb piece (if possible!) maximizing the sites score function
   bestseq.score=-1.e6;
   vint bestmotis(nbmots_for_score,0);

   //we define the 1kb piece with best score
   for (ivvinstseq iv=seq.instancescons.begin();iv!=seq.instancescons.end();iv++){

      vinstseq vic=*iv;
      vint bestmotistmp(nbmots_for_score,0);
      int start=seq.imaps[0][vic[0].pos]-neighbext;
      if (start<0) start=0;
      int stop=start+scanwidth-1;
      if (stop>seq.stop-seq.start){
         stop=seq.stop-seq.start;
         start=stop-scanwidth+1;
         if (start<0) start=0;
      }
      for (ivvinstseq iv1=seq.instancescons.begin();iv1!=seq.instancescons.end();iv1++){
         vinstseq vic1=*iv1;
         int start1=seq.imaps[0][vic1[0].pos];
         if (start1>=start && start1<stop-width+1){
            bestmotistmp[vic1[0].motindex]++;
         }
      }

      double score;
      int lseq=stop-start+1;
      score=calcscore(vm,bestmotistmp,lseq);

      if (score>bestseq.score){
         bestseq.score=score;
         bestseq.start=seq.start+start;
         bestseq.stop=seq.start+stop;
         bestmotis=bestmotistmp;
      }
   }

   if (seq.instancescons.size()==0) bestseq.score=calcscore(vm,bestmotis,seq.nbtb);

   //   cout << bestseq.name << " " << bestseq.score  <<" " << bestseq.start << " "  <<bestseq.stop << endl;
   //   cout << seq.name << " "  << seq.start << " "  <<seq.stop << endl;

   //cut original sequence into best piece
   int rawstart,rawstop;
   // Coorects minor bug due to extraction (incorrect stop due to gap etc...)
   if (bestseq.stop-seq.start+1>=seq.imapsinv[0].size()) bestseq.stop=seq.start+seq.imapsinv[0].size()-1;
   rawstart=seq.imapsinv[0][bestseq.start-seq.start];
   rawstop=seq.imapsinv[0][bestseq.stop-seq.start];
   for (int i=0;i<nbspecies;i++){
      if (seq.species[i]){
         civint istart=seq.iseqs[i].begin()+seq.imaps[i][rawstart];
         civint istop=seq.iseqs[i].begin()+seq.imaps[i][rawstop]+1;
         if (istop>=seq.iseqs[i].end()) istop=seq.iseqs[i].end();
         vint itmp(istart,istop);
         bestseq.iseqs[i]=itmp;
         bestseq.seqs[i]=vinttostring(itmp);
         bestseq.seqsrealigned[i]=seq.seqsrealigned[i].substr(rawstart,rawstop-rawstart+1);
         bestseq.imaps[i]=alignedtomap(bestseq.seqsrealigned[i]);
         bestseq.imapsinv[i]=alignedtorevmap(bestseq.seqsrealigned[i]);
         //         cout << bestseq.iseqs[i] << endl;
         //         cout << bestseq.seqs[i] << endl;
         //         cout << bestseq.seqsrealigned[i] << endl;
         int nbtb=bestseq.iseqs[i].size()-compN(bestseq.iseqs[i]);
         if (nbtb>width+neighbext){//>5){}
            bestseq.species[i]=1;
         } else {
            bestseq.species[i]=0;
         }
      }
   }
   bestseq.nbN=compN(bestseq.iseqs[0]);
   bestseq.nbtb=bestseq.iseqs[0].size()-bestseq.nbN;
   //   cout << bestseq.nbN << endl;
   //   cout << bestseq.nbtb << endl;

   return bestseq;

}

   vseq 
findbestseqs(vmot & vm, vseq & seqs1)
{
   // INIT
   vseq seqs=seqs1;

   // We need different enhancers of the same gene to have same name

   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
      if ((*ivs).name.find('-')!=string::npos)
      {
         (*ivs).name.erase((*ivs).name.begin()+(*ivs).name.find('-'),(*ivs).name.end());
      }
   }

   // SCANNING OLD SEQS
   //
   //   cout << "Scanning sequences for instances..." << endl;
   scanseqsforinstancesnmask(seqs,vm);

   // KEEPING BEST SEQS
   //
   //   cout << "Defining conserved instances..." << endl;
   //For each sequence find best part and cut it into vtmp
   vseq vtmp;
   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
      ivs->instances2instancescons();
      vtmp.push_back(findbestseq(vm,*ivs));
   }

   // UNICITY: ONE SEQUENCE PER GENE
   //
   vseq vbest;
   // only keep the best seq among all surrounding / intronic regions of a given gene
   if (vtmp.size()>1){
      Sequence bestseq=*(vtmp.begin());
      for (ivseq ivs=vtmp.begin()+1;ivs!=vtmp.end();ivs++){ 
         if (ivs->name==bestseq.name){ 
            if (ivs->score>bestseq.score){
               bestseq=*ivs;
            }
         } else{
            vbest.push_back(bestseq);
            bestseq=*ivs;
            if (ivs==vtmp.end()-1) vbest.push_back(*ivs);
         }
      }
   } else {
      vbest=vtmp;
   }

   // RE-SCANNING BEST SEQUENCES
   //
   for (ivseq iv=vbest.begin();iv!=vbest.end();iv++){
      iv->instances.clear();
      iv->instancescons.clear();
      vint dum(nbmots_for_score,0);
      iv->motis=dum;
   }
   //   cout << "Scanning best 1kbs for instances..." << endl;
   //scanseqsforinstancesnmask(vbest,vm);
   scanseqsforinstances(vbest,vm);

   //   cout << "Defining best 1kbs conserved instances..." << endl;
   for (ivseq ivs=vbest.begin();ivs!=vbest.end();ivs++){
      ivs->instances2instancescons();
      //   cout << ivs->name << endl;
      //  cout << ivs->instancescons << endl;
   }


   // RE-DEFINE SIZEPOS AND SIZENEG
   vseq vscorepos;
   sizepos=0;
   for (ivseq ivs=vbest.begin();ivs!=vbest.end();ivs++){
      if (ivs->sign==1){
         sizepos++;
      }
   }
   sizeneg=vbest.size()-sizepos;

   return vbest;
}

   void
dispbestseqs(vseq & seqs, ofstream & outf)
{
   for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
      outf << chromfromint(ivs->chrom) << "\t";
      outf << ivs->start << "\t";
      outf << ivs->stop << "\t";
      outf << ivs->name << "\n";
   }
}

   double
calcareaforscore(vseq & vscore)
{

   double area(0.);
   double pFP(0.),pTP(0.);
   double TP(0.),FP(0.);
   double pscore=vscore[0].score;
   for (ivseq ivs=vscore.begin()+1;ivs!=vscore.end();ivs++){
      if ((*ivs).sign==1){
         TP+=1.;
      }
      else if ((*ivs).sign==-1){
         FP+=1.;
      }
      if (ivs->score!=pscore){
         area+=(FP-pFP)*(TP+pTP)/2;
         pFP=FP;
         pTP=TP;
         pscore=ivs->score;
      }
      if (ivs->score==pscore && ivs==vscore.end()-1){
         area+=(FP-pFP)*(TP+pTP)/2;
      }
   }
   area/=(sizepos*sizeneg);
   return area;
}


   double
calcarea(vvd & rates)
{

   double area(0.);
   double FPtemp((*rates.begin())[0]);
   double TPtemp((*rates.begin())[1]);

   for (ivvd iv=rates.begin()+1;iv!=rates.end();iv++){
      double FP((*iv)[0]);
      double TP((*iv)[1]);
      if (TPtemp!=0 || TP!=0){
         area+=(FP-FPtemp)*(TPtemp+TP)/2;
      }
      FPtemp=FP;
      TPtemp=TP;
   }

   return area;
}

   vd 
calcinflex(vvd & rates)
{
   vd bestpoint(3,0.); //(FP,TP,thr)

   double dist(0.),distmax(0.),bestTP(0.),bestFP(0.);

   //rates for all thresholds, with values (FP,TP)
   for (ivvd iv=rates.begin();iv!=rates.end();iv++){
      double FP,TP,thr;
      FP=(*iv)[0];
      TP=(*iv)[1];
      thr=(*iv)[2];
      dist=(TP-FP)/sqrt(2);

      if (dist>distmax && thr>0.2*width){
         distmax=dist;
         bestpoint[0]=thr;
         bestpoint[1]=FP;
         bestpoint[2]=TP;
      }
   }

   return bestpoint;
}

//!! test sequences need to be >= width !!
   void    
bestpwm(vseq & vscore,vmot & vm)
{
   ofstream streamfile;
   streamfile.open("bestuniq-auc.dat");
   ofstream auc;
   auc.open("mot-auc.dat");

   //compute AUCs
   for (ivmot iv=vm.begin();iv!=vm.end();iv++){
      for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
         vint seqint=ivs->iseqs[0];
         ivs->score=scorefshift(seqint,iv->matprec);
         if (scorefshift(seqint,iv->matprecrevcomp)>ivs->score){
            ivs->score=scorefshift(seqint,iv->matprecrevcomp);
         } 
         //      cout << seqint << " " << ivs->score << endl;
      }
      sort(vscore.begin(),vscore.end());
      iv->optauc=calcareaforscore(vscore);
   }
   sort(vm.begin(),vm.end());

   for (ivmot iv=vm.begin();iv!=vm.end();iv++){
      auc << iv->index+1 << "\t" << iv->optauc << endl;
      cout << iv->index+1 << "\t" << iv->optauc << endl;
      iv->display(streamfile);
   }

   streamfile.close();
   auc.close();

}

   void    
rocinfo(vseq & vscore,vmot & vm)
{


   ofstream streamfile;
   streamfile.open("score/roc-infos.dat");
   streamfile << "Rank" << "\t";
   streamfile << "Nbmot" << "\t";
   streamfile << "AUC" << "\t";
   //   streamfile << "Score_cutoff" << "\t";
   streamfile << "NumFalsePos" << "\t";
   streamfile << "NumFalseNeg" << "\n";

   //compute AUCs
   for (ivmot iv=vm.begin();iv!=min(vm.begin()+nbmots_for_score,vm.end());iv++){

      for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
         calcscore(*iv,*ivs);
      }
      vseq vscorebest;      
      vmot vm;
      vm.push_back(*iv);
      vscorebest=findbestseqs(vm,vscore);
      //vscorebest=vscore;
      sort(vscorebest.begin(),vscorebest.end());
      //cout << "sizepos=" <<sizepos << " sizeneg=" << sizeneg << endl;

      double dist,area(0.),distmax(0.),FP(0.),TP(0.),pFP(0.),pTP(0.),pscore(vscorebest.begin()->score),scoremax(0.);

      for (ivseq ivs=vscorebest.begin();ivs!=vscorebest.end();ivs++){
         if ((*ivs).sign==1) TP+=1./sizepos;
         else if ((*ivs).sign==-1) FP+=1./sizeneg;
         dist=TP-FP;
         //         //cout << "dist " << dist << endl;
         //         if (dist<0){
         //            isneg=1;
         //         }
         //         //cout << ivs->score << endl;
         //         //// we do not want uninstructive, late peaks
         //         //// we also allow one wrong element in the beginning
         //         if (FP>1. && dist>distmax){
         //            // we avoid the situation where we select 0 motif
         //           // if (!(scoremax>0 && ivs->score==0)){
         //           if (!isneg){
         //               distmax=dist;
         //               //cout << "distmax " << dist << endl;
         //               scoremax=ivs->score;
         //               //cout << "scoremax=" << scoremax << endl;
         //            }
         //         //   }
         //         }
         if (ivs==vscorebest.begin() || ivs==vscorebest.end()-1){
            area+=(FP-pFP)*(TP+pTP)/2;
         }
         else {
            if (ivs->score!=pscore){
               area+=(FP-pFP)*(TP+pTP)/2;
               pFP=FP;
               pTP=TP;
            }
         }
         pscore=ivs->score;
      }

      for (double n=0;n<3;n++){
         int numFalsePos(0),numFalseNeg(0);
         //      cout << "motif " << iv->index << " scoremax=" << scoremax << " distmax=" << distmax << endl;
         for (ivseq ivs=vscorebest.begin();ivs!=vscorebest.end();ivs++){
            //         if (ivs->score>=scoremax && ivs->sign==-1) numFalsePos++;
            //         if (ivs->score<scoremax && ivs->sign==1) numFalseNeg++;
            if (ivs->motis[iv->index]>n && ivs->sign==-1) numFalsePos++;
            if (ivs->motis[iv->index]<=n && ivs->sign==1) numFalseNeg++;
         }

         streamfile << iv->index+1 << "\t";
         streamfile << n+1 << "\t";
         streamfile << area << "\t";
         //      streamfile << scoremax << "\t";
         streamfile << numFalsePos << "\t";
         streamfile << numFalseNeg << "\n";
      }
   }

   streamfile.close();

}

   void    
roconsites(vseq & vscore,vmot & vm)
{
   for (ivmot iv=vm.begin();iv!=vm.end();iv++){
      ofstream rocf;
      ostringstream rocfs;
      rocfs << "score/roc-on-sites/roc-" << iv->index << ".dat";
      rocf.open(rocfs.str().c_str());
      for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
         double sum=0;
         civint civ=ivs->iseqs[0].begin();
         ivs->score=scoref(civ,iv->matprec);
         if (scoref(civ,iv->matprecrevcomp)>ivs->score){
            ivs->score=scoref(civ,iv->matprecrevcomp);
         }
      }
      sort(vscore.begin(),vscore.end());
      double FP(0.),TP(0.),scoretmp;
      for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
         if ((*ivs).sign==1){
            TP+=1./sizepos;
         }
         else if ((*ivs).sign==-1){
            FP+=1./sizeneg;
         }
         if (ivs==vscore.begin()){
            rocf << ivs->score+1 << "\t" << 0 << "\t" << 0 << endl;
            rocf << ivs->score << "\t" << FP << "\t" << TP << endl;
         }
         if (ivs==vscore.end()-1){
            rocf << ivs->score << "\t" << FP << "\t" << TP << endl;
         }
         else {
            if (ivs->score!=scoretmp){
               rocf << ivs->score << "\t" << FP << "\t" << TP << endl;
            }
         }
         scoretmp=ivs->score;
      }
      rocf.close();
   }

}

void correlations(vseq & vscore,vmot & vm)
{
   system("if ! test -d score/correlations;then mkdir score/correlations;fi;");      
   vseq vscoretmp;
   //   keeping only same size sequences
   int sizerejpos(0),sizerejneg(0);
   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
      if (ivs->nbtb==2000){
         vscoretmp.push_back(*ivs);
      }
      else{
         if (ivs->sign==1) sizerejpos++;
         if (ivs->sign==-1) sizerejneg++;
      }
   }
   cout << "Rejected " << sizerejpos << " positive and " << sizerejneg << " negative promoters too short" << endl;
   //   exit(9);

   // stopping at 50 pos 50 neg
   for (ivmot iv=vm.begin();iv!=min(vm.end(),vm.begin()+nbmots_for_score);iv++){ 
      vscore=vscoretmp;
      cout << "Motif " << iv->index << endl;
      for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
         double sum=0;
         calcscore(*iv,*ivs);
      }
      sort(vscore.begin(),vscore.end());
      double FP(0.),TP(0.);
      for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
         if ((*ivs).sign==1) TP+=1.;
         else if ((*ivs).sign==-1) FP+=1.;

         if (TP<=50){
            for (ivseq ivs1=vscoretmp.begin(); ivs1!=vscoretmp.end();ivs1++){
               if (ivs1->name==ivs->name){
                  ivs1->signpermot[iv->index]=1;
               }
            }
         } else{
            for (ivseq ivs1=vscoretmp.begin(); ivs1!=vscoretmp.end();ivs1++){
               if (ivs1->name==ivs->name){
                  ivs1->signpermot[iv->index]=-1;
               }
            }
         }
      }
   }

   ofstream auc;
   auc.open("score/correlations/correlations.dat");
   for (ivseq ivs=vscoretmp.begin();ivs!=vscoretmp.end();ivs++){
      if (ivs->sign==1){
         cout << ivs->name << "\t";
         auc << ivs->name << "\t";
         for (int i=0;i<nbmots_for_score;i++){
            cout << ivs->signpermot[i] << "\t";
            auc << ivs->signpermot[i] << "\t";
         }
         cout << endl;
         auc << endl;
      }
   }
   auc.close();

   double corrmoy(0),corrvar(0);
   double counter(0);
   for (int i=0;i<1;i++){//nbmots_for_score-1;i++){}
      for (int j=i+1;j<nbmots_for_score;j++){
         double moy1(0);
         double moy2(0);
         double prod12(0);
         for (ivseq ivs=vscoretmp.begin();ivs!=vscoretmp.end();ivs++){
            if (ivs->sign==1){
               moy1+=ivs->signpermot[i];
               moy2+=ivs->signpermot[j];
               prod12+=ivs->signpermot[i]*ivs->signpermot[j];
            }
         }
         moy1/=sizepos;
         moy2/=sizepos;
         prod12/=sizepos;
         double correlation;
         correlation=prod12-moy1*moy2;
         cout << "Motif "<< i << "VS" << j <<
            " Moy" << i << " " << moy1 <<
            " Moy" << j << " " << moy2 <<
            " prod " << i << j << " " << prod12 <<
            " corr " << correlation << endl;
         corrmoy+=correlation;
         corrvar+=correlation*correlation;
         counter++;
      }
   }
   cout << counter << endl;
   corrmoy/=counter;
   corrvar/=counter;
   corrvar-=corrmoy*corrmoy;
   corrvar=sqrt(corrvar);
   cout << "Correlation mean: " << corrmoy << endl;
   cout << "Correlation variance: " << corrvar << endl;

}

void rocpermot(vseq & vscore,vmot & vm)
{
   system("if ! test -d score/roc-per-mot;then mkdir score/roc-per-mot;fi;");      
   //   vseq vscoretmp;
   //   //   keeping only same size sequences
   ////   int sizerejpos(0),sizerejneg(0);
   //   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
   //  //    if (ivs->nbtb<=2000){
   //  //       vscoretmp.push_back(*ivs);
   //    //  }
   ////      else{
   ////         if (ivs->sign==1) sizerejpos++;
   ////         if (ivs->sign==-1) sizerejneg++;
   ////      }
   //   }
   // //  cout << "Rejected " << sizerejpos << " positive and " << sizerejneg << " negative promoters too short" << endl;
   //   //   exit(9);
   //   //vscore=vscoretmp;

   // computing roc for each motif
   ofstream auc;
   auc.open("score/roc-per-mot/roc_auc.dat");
   for (ivmot iv=vm.begin();iv!=vm.begin()+nbmots_for_score;iv++){ 
      ofstream outf;
      ofstream outf1;
      ostringstream outfs;
      ostringstream outfs1;
      cout << "Motif " << iv->index+1 << endl;
      width=iv->bsinit.size();
      scorethr2=iv->motscorethr2;
      scorethr=iv->motscorethr;
      scorethrcons=iv->motscorethrcons;
      cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;
      double area(0);
      ofstream rocf;
      ostringstream rocfs;
      rocfs << "score/roc-per-mot/roc_" << iv->index << ".dat";
      outfs << "score/roc-per-mot/bestseqs-coords-" << iv->index <<".dat";
      outfs1 << "score/roc-per-mot/bestseqs-" << iv->index << ".dat";
      outf.open(outfs.str().c_str());
      outf1.open(outfs1.str().c_str());
      rocf.open(rocfs.str().c_str());
      //      ofstream rocrf;
      //      ostringstream rocrfs;
      //      rocrfs << "score/roc-per-mot/rocr_" << iv->index << ".dat";
      //      rocrf.open(rocrfs.str().c_str());

      for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
         if (args_info.wocons_given) ivs->score=iv->nbmatchmat(*ivs);
         else ivs->score=iv->nbmatchcons(*ivs);
      }
      vseq vscorebest;      
      vmot vm;
      vm.push_back(*iv);
      //vscorebest=findbestseqs(vm,vscore);
      vscorebest=vscore;
      sort(vscorebest.begin(),vscorebest.end());
      vseq vscorepos;
      for (ivseq ivs=vscorebest.begin();ivs!=vscorebest.end();ivs++){
         if (ivs->sign==1){
            vscorepos.push_back(*ivs);
         }
      }
      dispbestseqs(vscorepos,outf);
      outf1 << vscorepos << endl;
      //cout << vscorepos << endl;
      //cout << vscorebest << endl;

      sizepos=vscorepos.size();
      sizeneg=vscorebest.size()-sizepos;
      //cout << sizepos << " " << sizeneg << endl;
      double FP(0.),TP(0.),pFP(0.),pTP(0.),scoretmp;
      rocf << vscorebest.begin()->score+1 << "\t" << FP << "\t" << TP << endl;
      for (ivseq ivs=vscorebest.begin();ivs!=vscorebest.end();ivs++){
         //         rocrf << ivs->score << "\t" << ivs->sign << endl;
         //      cout << ivs->nbtb << " ";
         if ((*ivs).sign==1) TP+=1./sizepos;
         else if ((*ivs).sign==-1) FP+=1./sizeneg;

         if (ivs==vscorebest.begin() || ivs==vscorebest.end()-1){
            rocf << ivs->score << "\t" << FP << "\t" << TP << endl;
            area+=(FP-pFP)*(TP+pTP)/2;
         }
         else {
            if (ivs->score!=scoretmp){
               rocf << ivs->score << "\t" << FP << "\t" << TP << endl;
               area+=(FP-pFP)*(TP+pTP)/2;
               pFP=FP;
               pTP=TP;
            }
         }
         scoretmp=ivs->score;
      }
      auc << iv->index << "\t" << area << endl;

      rocf.close();
      outf.close();
      outf1.close();
      //      rocrf.close();
   }
   auc.close();
}

void rocpermotvarthr(vseq & vscore,vmot & vm)
{
   system("if ! test -d score/roc-per-mot-var-thr;then mkdir score/roc-per-mot-var-thr;fi;");      

   // computing roc for each motif
   ofstream auc;
   auc.open("score/roc-per-mot-var-thr/roc_auc.dat");

   for (ivmot iv=vm.begin();iv!=vm.begin()+nbmots_for_score;iv++){ 

      cout << "Motif " << iv->index << endl;
      width=iv->bsinit.size();
      ofstream rocf;
      ostringstream rocfs;
      rocfs << "score/roc-per-mot-var-thr/roc_" << iv->index << ".dat";
      rocf.open(rocfs.str().c_str());



      double FP(0.),TP(0.),pFP(0.),pTP(0.),area(0);
      rocf << 12. << "\t" << 0. << "\t" << 0. << endl;

      for (double thr=11.5;thr>0.;thr--){

         FP=0;
         TP=0;   
         iv->motscorethr2=iv->motwidth*thr/10;
         iv->motscorethrcons=iv->motscorethr2;

         for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
            if (args_info.wocons_given){
               if (iv->nbmatchmat(*ivs)==0) continue;
            } else {
               if (iv->nbmatchcons(*ivs)==0) continue;
            }
            if ((*ivs).sign==1) TP+=1./sizepos;
            else if ((*ivs).sign==-1) FP+=1./sizeneg;
         }
         rocf << thr << "\t";
         rocf<< FP << "\t";
         rocf << TP << "\n";
         area+=(FP-pFP)*(TP+pTP)/2;
         pFP=FP;
         pTP=TP;
      }
      rocf << -10 << "\t";
      rocf<< 1. << "\t";
      rocf << 1. << "\n";
      area+=(1.-pFP)*(1.+pTP)/2;
      auc << iv->index << "\t" << area << endl;

      rocf.close();
   }
   auc.close();
}

void roccumulmotvarthr(vseq & vscore,vmot & vm)
{
   system("if ! test -d score/roc-cumul-mot-var-thr;then mkdir score/roc-cumul-mot-var-thr;fi;");      

   // computing roc for each motif
   ofstream auc("score/roc-cumul-mot-var-thr/roc_auc.dat");
   ofstream rocf("score/roc-cumul-mot-var-thr/roc_cumul.dat");

   double FP(0.),TP(0.),pFP(0.),pTP(0.),area(0);
   rocf << 12. << "\t" << 0. << "\t" << 0. << endl;

   for (double thr=11.5;thr>0.;thr--){

      FP=0;
      TP=0;   

      for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){

         unsigned int nummatch(0);
         for (ivmot iv=vm.begin();iv!=vm.begin()+nbmots_for_score;iv++){ 
            width=iv->bsinit.size();
            iv->motscorethr2=iv->motwidth*thr/10;
            iv->motscorethrcons=iv->motscorethr2;

            if (args_info.wocons_given){
               nummatch+=iv->nbmatchmat(*ivs);
            } else {
               nummatch+=iv->nbmatchcons(*ivs);
            }
         }

         //if (nummatch==0) continue; // OR
         if (nummatch<nbmots_for_score) continue; // AND

         if ((*ivs).sign==1) TP+=1./sizepos;
         else if ((*ivs).sign==-1) FP+=1./sizeneg;
      }
      rocf << thr << "\t";
      rocf<< FP << "\t";
      rocf << TP << "\n";
      area+=(FP-pFP)*(TP+pTP)/2;
      pFP=FP;
      pTP=TP;
   }
   rocf << -10 << "\t";
   rocf<< 1. << "\t";
   rocf << 1. << "\n";
   area+=(1.-pFP)*(1.+pTP)/2;
   auc << "cumul" << "\t" << area << endl;

   rocf.close();
   auc.close();
}

void roccumulmot(vseq & vscore,vmot & vm)
{
   system("if ! test -d score/roc-cumul-mot;then mkdir score/roc-cumul-mot;fi;");      

   // computing roc for cumulated motif
   ofstream auc;
   auc.open("score/roc-cumul-mot/roc_auc.dat");
   for (unsigned int i=1;i<=nbmots_for_score;i++){ 
      cout << "Cumul " << i << endl;

      double area(0);
      ofstream outf;
      ofstream outf1;
      ofstream rocf;
      ostringstream outfs;
      ostringstream outfs1;
      ostringstream rocfs;
      rocfs << "score/roc-cumul-mot/roc_" << i << ".dat";
      outfs << "score/roc-cumul-mot/bestseqs-coords-" << i <<".dat";
      outfs1 << "score/roc-cumul-mot/bestseqs-" << i << ".dat";
      outf.open(outfs.str().c_str());
      outf1.open(outfs1.str().c_str());
      rocf.open(rocfs.str().c_str());

      vseq vscorebest;      
      vmot vmtmp;
      for (int j=0;j<i;j++){
         vmtmp.push_back(vm[j]);
      }
      vscorebest=findbestseqs(vmtmp,vscore);
      sort(vscorebest.begin(),vscorebest.end());
      vseq vscorepos;
      for (ivseq ivs=vscorebest.begin();ivs!=vscorebest.end();ivs++){
         if (ivs->sign==1){
            vscorepos.push_back(*ivs);
         }
      }
      dispbestseqs(vscorepos,outf);
      outf1 << vscorepos << endl;
      //cout << vscorepos << endl;
      cout << vscorebest << endl;

      sizepos=vscorepos.size();
      sizeneg=vscorebest.size()-sizepos;
      double FP(0.),TP(0.),pFP(0.),pTP(0.),scoretmp;
      rocf << vscorebest.begin()->score+1 << "\t" << FP << "\t" << TP << endl;
      for (ivseq ivs=vscorebest.begin();ivs!=vscorebest.end();ivs++){
         if ((*ivs).sign==1) TP+=1./sizepos;
         else if ((*ivs).sign==-1) FP+=1./sizeneg;

         if (ivs==vscorebest.begin() || ivs==vscorebest.end()-1){
            rocf << ivs->score << "\t" << FP << "\t" << TP << endl;
            area+=(FP-pFP)*(TP+pTP)/2;
         }
         else {
            if (ivs->score!=scoretmp){
               rocf << ivs->score << "\t" << FP << "\t" << TP << endl;
               area+=(FP-pFP)*(TP+pTP)/2;
               pFP=FP;
               pTP=TP;
            }
         }
         scoretmp=ivs->score;
      }
      auc << i << "\t" << area << endl;

      rocf.close();
      outf.close();
      outf1.close();
      //      rocrf.close();
   }
   auc.close();
}

//Computes best threshold detection, ROC area, new lambdas
//and scores at optimal thr
   vseq
motoptibestseqs(vseq & vscore, Motif & mot)
{
   ofstream auc;
   ostringstream aucs;
   system("if ! test -d score/auc;then mkdir score/auc;fi;");      
   aucs << "score/auc/auc_" << mot.index << ".dat";
   auc.open(aucs.str().c_str());
   double thrmax(10),thrmin(2);
   double bestauc(0),bestthr(thrmin),bestlt,bestlb;
   vseq vscorebest;
   vmot vtmp;
   vtmp.push_back(mot);
   for (double thr=thrmax;thr>thrmin;thr-=.5){
      //vtmp for findbest1kb's use
      vtmp[0].motscorethr2=mot.bsinit.size()*thr/10;
      vtmp[0].motscorethr=mot.bsinit.size()*(thr-1.)/10;
      vtmp[0].motscorethrcons=mot.bsinit.size()*(thr-1.)/10;
      //      mot.calclambda();
      if (mot.lambda==0 || mot.lambdatrain==0){
         //        cout << "Mot #" << mot.index << " Thr= "<< thr << "\tBeware: lambda goes to zero ! " << "lambda_t = " << mot.lambdatrain << " lambda_b = " << mot.lambda << endl;
         continue;
      }
      vseq vscoretmp;
      vscoretmp=findbestseqs(vtmp,vscore);
      sort(vscoretmp.begin(),vscoretmp.end());
      double area;
      area=calcareaforscore(vscoretmp);
      auc << thr << "\t" << area << endl;
      if (area > bestauc){
         bestauc=area;
         bestthr=thr;
         bestlt=mot.lambdatrain;
         bestlb=mot.lambda;
         vscorebest=vscoretmp;
      }
   }
   auc.close();

   mot.motscorethr2=bestthr;
   mot.motscorethr=mot.bsinit.size()*(bestthr-1.)/10;
   mot.motscorethrcons=mot.bsinit.size()*(bestthr-1.)/10;
   //   mot.lambda=bestlb;
   //   mot.lambdatrain=bestlt;
   mot.optthr=bestthr;
   mot.optauc=bestauc;

   return vscorebest;
}


//Motif independent ROCs at optimnal thr
   void
rocopti(vseq & vscore, Motif & mot)
{

   ofstream rocf;
   ostringstream rocfs;
   rocfs << "score/rocopti/roc_" << mot.index << ".dat";
   rocf.open(rocfs.str().c_str());
   sort(vscore.begin(),vscore.end());

   double FP(0.),TP(0.),pFP(0.),pTP(0.),scoretmp;

   rocf << vscore.begin()->score+1 << "\t" << FP << "\t" << TP << endl;
   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
      if ((*ivs).sign==1) TP+=1./sizepos;
      else if ((*ivs).sign==-1) FP+=1./sizeneg;
      if (ivs==vscore.begin() || ivs==vscore.end()-1){
         rocf << ivs->score << "\t" << FP << "\t" << TP << endl;
      }
      else {
         if (ivs->score!=scoretmp){
            rocf << ivs->score << "\t" << FP << "\t" << TP << endl;
            pFP=FP;
            pTP=TP;
         }
      }
      scoretmp=ivs->score;
   }
   rocf.close();
}

//Motif independent optimization
   void
motopti(vseq & vscore)
{
   //   system("if ! test -d score/rocopti;then mkdir score/rocopti;fi;");      
   //   system("if ! test -d score/dispscore;then mkdir score/dispscore;fi;");      

   ofstream outf("score/motopti.dat");
   outf << "Motif\tAUC\tThr\tlambda_train\tlambda_back\n";
   for (ivmot im=motsdef.begin();im!=motsdef.begin()+nbmots_for_score;im++){

      vseq vscorebest;
      vscorebest=motoptibestseqs(vscore,*im);

      // display
      outf << im->index <<  "\t" << im->optthr << "\t" << im->optauc << "\t" << im->lambdatrain << "\t" << im->lambda << endl;
      cout << "Motif " << im->index << "\t" << "AUC " << im->optauc << "\t" << "thr " << im->optthr << "\tlambda_train " << im->lambdatrain << "\tlambda_back " << im->lambda << endl;

      rocopti(vscorebest,*im);

      //Display scores at optimal point for each motif 
      //      dispscore(vscorebest,*im);
   }
   outf.close();
}

//Motif independent ROCs at optimnal thr
   void
roccombine(vseq & vscore)
{

   ofstream rocf;
   ostringstream rocfs;
   rocfs << "score/roc-combine/roc_";
   for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      rocfs << ivm->index;
      if (ivm!=motsdef.end()-1){
         rocfs << "+";
      }
   }
   rocfs << ".dat";
   rocf.open(rocfs.str().c_str());


   double FP(0.),TP(0.),scoretmp;

   rocf << vscore.begin()->score << "\t" << FP << "\t" << TP << endl;
   scoretmp=vscore.begin()->score;
   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
      if ((*ivs).sign==1) TP+=1./sizepos;
      else if ((*ivs).sign==-1) FP+=1./sizeneg;
      if (ivs!=vscore.begin()){
         if (ivs->score!=scoretmp){
            rocf << ivs->score << "\t" << FP << "\t" << TP << endl;
         }
         else {
            if (ivs==vscore.end()-1){
               rocf << ivs->score << "\t" << FP << "\t" << TP << endl;
            }
         }
      }
      scoretmp=ivs->score;
   }
   rocf.close();

}

   void
loadseqsforscore(vseq & vscore)
{

   if (args_info.posalign_given){ 
      ifstream inf;
      inf.open(args_info.posalign_arg);
      vseq seqs;
      seqs=loadsequencesconserv(inf);
      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
         if (ivs->species[0] && ivs->nbtb>0){
            ivs->sign=1;
            vscore.push_back(*ivs);
         }
      }
      sizepos=vscore.size();
      cout << "Positive set: " << sizepos << "\n";
      inf.close();
   }

   if (args_info.poscoords_given){
      //Working on coordinate files
      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align-masked.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align-masked.dat");

      //POSITIVES
      ifstream inf;
      inf.open(args_info.poscoords_arg);
      alignscoord=loadcoordconserv(align);
      align.close();
      vcoord pcoords;
      back_insert_iterator<vcoord> pdest(pcoords);
      copy(iiscoord(inf),iiscoord(),pdest);
      for (ivcoord ivc=pcoords.begin();ivc!=pcoords.end();ivc++){
         Sequence seqtoimport=coordtoseq(*ivc);
         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
            seqtoimport.sign=1;
            vscore.push_back(seqtoimport);
         }
      }
      inf.close();
      sizepos=vscore.size();
      cout << "Positive set: " << sizepos << "\n";
   }

   if (args_info.negalign_given){
      ifstream inf;
      inf.open(args_info.negalign_arg);
      vseq seqs;
      seqs=loadsequencesconserv(inf);
      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
         if (ivs->species[0] && ivs->nbtb>0){
            ivs->sign=-1;
            vscore.push_back(*ivs);
         }
      }
      sizeneg=vscore.size()-sizepos;
      cout << "Negative set: " << sizeneg << "\n";
      inf.close();
   }

   if (args_info.negcoords_given){
      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align-masked.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align-masked.dat");

      //NEGATIVES
      ifstream inf;
      inf.open(args_info.negcoords_arg);
      alignscoord=loadcoordconserv(align);
      align.close();
      vcoord ncoords;
      back_insert_iterator<vcoord> ndest(ncoords);
      copy(iiscoord(inf),iiscoord(),ndest);
      for (ivcoord ivc=ncoords.begin();ivc!=ncoords.end();ivc++){
         Sequence seqtoimport=coordtoseq(*ivc);
         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
            seqtoimport.sign=-1;
            vscore.push_back(seqtoimport);
         }
      }
      inf.close();
      sizeneg=vscore.size()-sizepos;
      cout << "Negative set: " << sizeneg << "\n";
   }

   if (!(args_info.negcoords_given||args_info.negalign_given)){
      cout << "Using background as negative" << endl;
      ifstream inf;
      if (species==1) inf.open("/home/santolin/these/files/droso/backreg/regs2000-align-1000.dat"); 
      else if (species==2) inf.open("/home/santolin/these/files/mus/backreg/regs2000-align-1000.dat");//regs2000.fa");
      vseq seqs;
      seqs=loadsequencesconserv(inf);
      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
         if (ivs->species[0] && ivs->nbtb>0){
            ivs->sign=-1;
            vscore.push_back(*ivs);
         }
      }
      sizeneg=vscore.size()-sizepos;
      cout << "Negative set: " << sizeneg << "\n";
      inf.close();
   }

   if (!(args_info.poscoords_given || args_info.posalign_given)) {
      cout << "Please give positive file" << endl;
      exit(1);
   }

   //   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
   //      if (!(ivs->sign==-1 && ivs->name.find("reg")!=string::npos)){
   ////         if ((*ivs).name.find('_')!=string::npos)
   ////         {
   ////            if (args_info.display_given || args_info.byscore_given){
   ////               (*ivs).name.insert((*ivs).name.find('_'),"\\");
   ////            } else{
   ////               (*ivs).name.erase((*ivs).name.begin()+(*ivs).name.find('_'),(*ivs).name.end());
   ////            }
   ////         }
   //         //      if ((*ivs).name.find('-')!=string::npos)
   //         //      {
   //         //         (*ivs).name.erase((*ivs).name.begin()+(*ivs).name.find('-'),(*ivs).name.end());
   //         //      }
   //      }
   //   }

   return;
}

   void
loadseqsforscorenotmasked(vseq & vscore)
{

   if (args_info.posalign_given){ 
      ifstream inf;
      inf.open(args_info.posalign_arg);
      vseq seqs;
      seqs=loadsequencesconserv(inf);
      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
         if (ivs->species[0] && ivs->nbtb>0){
            ivs->sign=1;
            vscore.push_back(*ivs);
         }
      }
      sizepos=vscore.size();
      cout << "Positive set: " << sizepos << "\n";
      inf.close();
   }

   if (args_info.poscoords_given){
      //Working on coordinate files
      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align.dat");

      //POSITIVES
      ifstream inf;
      inf.open(args_info.poscoords_arg);
      alignscoord=loadcoordconserv(align);
      align.close();
      vcoord pcoords;
      back_insert_iterator<vcoord> pdest(pcoords);
      copy(iiscoord(inf),iiscoord(),pdest);
      for (ivcoord ivc=pcoords.begin();ivc!=pcoords.end();ivc++){
         Sequence seqtoimport=coordtoseq(*ivc);
         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
            seqtoimport.sign=1;
            vscore.push_back(seqtoimport);
         }
      }
      inf.close();
      sizepos=vscore.size();
      cout << "Positive set: " << sizepos << "\n";
   }

   if (args_info.negalign_given){
      ifstream inf;
      inf.open(args_info.negalign_arg);
      vseq seqs;
      seqs=loadsequencesconserv(inf);
      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
         if (ivs->species[0] && ivs->nbtb>0){
            ivs->sign=-1;
            vscore.push_back(*ivs);
         }
      }
      sizeneg=vscore.size()-sizepos;
      cout << "Negative set: " << sizeneg << "\n";
      inf.close();
   }

   if (args_info.negcoords_given){
      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align.dat");

      //NEGATIVES
      ifstream inf;
      inf.open(args_info.negcoords_arg);
      alignscoord=loadcoordconserv(align);
      align.close();
      vcoord ncoords;
      back_insert_iterator<vcoord> ndest(ncoords);
      copy(iiscoord(inf),iiscoord(),ndest);
      for (ivcoord ivc=ncoords.begin();ivc!=ncoords.end();ivc++){
         Sequence seqtoimport=coordtoseq(*ivc);
         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
            seqtoimport.sign=-1;
            vscore.push_back(seqtoimport);
         }
      }
      inf.close();
      sizeneg=vscore.size()-sizepos;
      cout << "Negative set: " << sizeneg << "\n";
   }

   if (!(args_info.negcoords_given||args_info.negalign_given)){
      cout << "Using background as negative" << endl;
      ifstream inf;
      if (species==1) inf.open("/home/santolin/these/files/droso/backreg/regs2000-align-1000.dat"); 
      else if (species==2) inf.open("/home/santolin/these/files/mus/backreg/regs2000-align-1000.dat");//regs2000.fa");
      vseq seqs;
      seqs=loadsequencesconserv(inf);
      for (ivseq ivs=seqs.begin();ivs!=seqs.end();ivs++){
         if (ivs->species[0] && ivs->nbtb>0){
            ivs->sign=-1;
            vscore.push_back(*ivs);
         }
      }
      sizeneg=vscore.size()-sizepos;
      cout << "Negative set: " << sizeneg << "\n";
      inf.close();
   }

   if (!(args_info.poscoords_given || args_info.posalign_given)) {
      cout << "Please give positive file" << endl;
      exit(1);
   }

   //   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
   //      if (!(ivs->sign==-1 && ivs->name.find("reg")!=string::npos)){
   ////         if ((*ivs).name.find('_')!=string::npos)
   ////         {
   ////            if (args_info.display_given || args_info.byscore_given){
   ////               (*ivs).name.insert((*ivs).name.find('_'),"\\");
   ////            } else{
   ////               (*ivs).name.erase((*ivs).name.begin()+(*ivs).name.find('_'),(*ivs).name.end());
   ////            }
   ////         }
   //         //      if ((*ivs).name.find('-')!=string::npos)
   //         //      {
   //         //         (*ivs).name.erase((*ivs).name.begin()+(*ivs).name.find('-'),(*ivs).name.end());
   //         //      }
   //      }
   //   }

   return;
}

//Loads align-file (fasta) or coord-file (name/chrom/start/stop)
   vseq
loadseqs()
{
   vseq seqs;

   if (args_info.align_file_given){

      ifstream inf;
      inf.open(args_info.align_file_arg);
      seqs=loadsequencesconserv(inf);
      inf.close();
   }
   else if (args_info.coord_file_given){

      ifstream inf;
      inf.open(args_info.coord_file_arg);

      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align-masked.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align-masked.dat");

      alignscoord=loadcoordconserv(align);

      align.close();

      vcoord coords;
      back_insert_iterator<vcoord> dest(coords);
      copy(iiscoord(inf),iiscoord(),dest);
      for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
         Sequence seqtoimport=coordtoseq(*ivc);
         //         cout << (*ivc).name << " " << (*ivc).start << endl;
         //         cout << seqtoimport.name << " " << seqtoimport.start << endl;
         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
            seqs.push_back(seqtoimport);
         }
         //         for (int i=0;i<nbspecies;i++){
         //         cout << seqtoimport.iseqs[i] << endl; 
         //         }
         //         exit(9);
      }
      inf.close();
   }
   else {
      cout << "Error in loadseqs: please give a coord/alignment file. Exiting..." << endl;
      exit(1);
   }
   return seqs;
}

//Loads align-file (fasta) or coord-file (name/chrom/start/stop)
   vseq
loadseqscons()
{
   vseq seqs;

   if (args_info.align_file_given){

      ifstream inf;
      inf.open(args_info.align_file_arg);
      seqs=loadsequencesconservonly(inf);
      inf.close();
   }
   else if (args_info.coord_file_given){

      ifstream inf;
      inf.open(args_info.coord_file_arg);

      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align-masked.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align-masked.dat");

      alignscoord=loadcoordconserv(align);

      align.close();

      vcoord coords;
      back_insert_iterator<vcoord> dest(coords);
      copy(iiscoord(inf),iiscoord(),dest);
      for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
         Sequence seqtoimport=coordtoseq(*ivc);
         //         cout << (*ivc).name << " " << (*ivc).start << endl;
         //         cout << seqtoimport.name << " " << seqtoimport.start << endl;
         vint spe;
         int i=0;
         for (ivint iv=seqtoimport.species.begin();iv!=seqtoimport.species.end();iv++){
            if (*iv==1) spe.push_back(i);
            i++;
         }
         if (iscons(spe) && seqtoimport.species[0] && seqtoimport.nbtb>0){
            seqs.push_back(seqtoimport);
         }
         if (!iscons(spe)) cout << ivc->name << ": sequence not conserved\n";
         //         for (int i=0;i<nbspecies;i++){
         //         cout << seqtoimport.iseqs[i] << endl; 
         //         }
         //         exit(9);
      }
      inf.close();
   }
   else {
      cout << "Error in loadseqs: please give a coord/alignment file. Exiting..." << endl;
      exit(1);
   }
   return seqs;
}

//Loads align-file (fasta) or coord-file (name/chrom/start/stop)
   vseq
loadseqsnotmasked()
{
   vseq seqs;

   if (args_info.align_file_given){

      ifstream inf;
      inf.open(args_info.align_file_arg);
      seqs=loadsequencesconserv(inf);
      inf.close();
   }
   else if (args_info.coord_file_given){

      ifstream inf;
      inf.open(args_info.coord_file_arg);

      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align.dat");

      alignscoord=loadcoordconserv(align);

      align.close();

      vcoord coords;
      back_insert_iterator<vcoord> dest(coords);
      copy(iiscoord(inf),iiscoord(),dest);
      for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
         //         cout << (*ivc).name << " " << (*ivc).start << endl;
         Sequence seqtoimport=coordtoseq(*ivc);
         //         cout << seqtoimport.name << " " << seqtoimport.start << endl;
         //         for (int i=0;i<nbspecies;i++){
         //         cout << seqtoimport.iseqs[i] << endl; 
         //         }
         //         exit(9);
         vint spe;
         int i=0;
         for (ivint iv=seqtoimport.species.begin();iv!=seqtoimport.species.end();iv++){
            if (*iv==1) spe.push_back(i);
            i++;
         }
         if (seqtoimport.species[0] && seqtoimport.nbtb>0){
            if (iscons(spe)){
               seqtoimport.cons=1;
            }
            else {
               //cout << ivc->name << ": sequence not conserved\n";
            }
            //      cout << "\r" <<  ivc->name;
            //      cout.flush();
            seqs.push_back(seqtoimport);
         }
      }
      inf.close();
   }
   else {
      cout << "Error in loadseqs: please give a coord/alignment file. Exiting..." << endl;
      exit(1);
   }
   return seqs;
}

//Loads align-file (fasta) or coord-file (name/chrom/start/stop) 
// keep only cons alignments!!
   vseq
loadseqsnotmaskedcons()
{
   vseq seqs;

   if (args_info.align_file_given){

      ifstream inf;
      inf.open(args_info.align_file_arg);
      seqs=loadsequencesconservonly(inf);
      inf.close();
   }
   else if (args_info.coord_file_given){

      ifstream inf;
      inf.open(args_info.coord_file_arg);

      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align.dat");

      alignscoord=loadcoordconserv(align);

      align.close();

      vcoord coords;
      back_insert_iterator<vcoord> dest(coords);
      copy(iiscoord(inf),iiscoord(),dest);
      for (ivcoord ivc=coords.begin();ivc!=coords.end();ivc++){
         Sequence seqtoimport=coordtoseq(*ivc);
         //         cout << (*ivc).name << " " << (*ivc).start << endl;
         //         cout << seqtoimport.name << " " << seqtoimport.start << endl;
         vint spe;
         int i=0;
         for (ivint iv=seqtoimport.species.begin();iv!=seqtoimport.species.end();iv++){
            if (*iv==1) spe.push_back(i);
            i++;
         }
         if (iscons(spe) && seqtoimport.species[0] && seqtoimport.nbtb>0){
            seqs.push_back(seqtoimport);
         }
         if (!iscons(spe)) cout << ivc->name << ": sequence not conserved\n";
      }
      inf.close();
   }
   else {
      cout << "Error in loadseqs: please give a coord/alignment file. Exiting..." << endl;
      exit(1);
   }
   return seqs;
}

   void 
statcounts(vseq & align, vmot & mots)
{

   vseq vscore;
   for (ivseq iv=align.begin();iv!=align.end();iv++){
      if ((*iv).name.find('_')!=string::npos)
      {
         (*iv).name.insert((*iv).name.find('_'),"\\");
      }
      if ((*iv).name.find('-')!=string::npos)
      {
         (*iv).name.erase((*iv).name.begin()+(*iv).name.find('-'),(*iv).name.end());
      }
      vscore.push_back(*iv);
   }
   vseq vbest;
   vbest=findbestseqs(mots,vscore);

   ofstream outf("scores.dat");
   for (ivseq ivs=vbest.begin();ivs!=vbest.end();ivs++){
      outf << ivs->score << endl;
   }
   outf.close();
}

   void 
statshifts(vseq & align, vmot & mots)
{
   vinstseq vshifts;
   for (ivseq iv=align.begin();iv!=align.end();iv++){
      vinstseq vinit;
      vinstseq vcomp;
      //We initialize shift at 0 for ref and scanwidth for other species
      for (ivinstseq ivs=iv->instances.begin();ivs!=iv->instances.end();ivs++){
         ivs->shift=scanwidth;
         if (ivs->species==0){
            ivs->shift=0;
            vinit.push_back(*ivs);
         } else{
            ivs->shift=scanwidth;
            vcomp.push_back(*ivs);
         }
      }
      //We compute the minimal shift for each instance
      for (ivinstseq ivi=vinit.begin();ivi!=vinit.end();ivi++){
         int mindist(scanwidth);
         int i=0;
         int ibest=-1;
         for (ivinstseq ivc=vcomp.begin();ivc!=vcomp.end();ivc++){
            if (ivc->motindex==ivi->motindex){
               int dist=ivc->pos-ivi->pos;
               if (fabs(dist)<fabs(mindist)){
                  mindist=dist;
                  ibest=i;
               }
            }
            i++;
         }
         if (ibest>-1){
            if (fabs(mindist)<fabs(vcomp[ibest].shift)){
               vcomp[ibest].shift=mindist;
            }
         }
      }

      for (ivinstseq ivc=vcomp.begin();ivc!=vcomp.end();ivc++){
         if (fabs(ivc->shift)<scanwidth){
            vshifts.push_back(*ivc);
         }
      }
   }

   for (ivmot ivm=mots.begin();ivm!=mots.begin()+nbmots_for_score;ivm++){
      for (int i=1;i<nbspecies;i++){

         ofstream outf;
         ostringstream outs;
         outs << "data/shift_" << numtospecies(i)<<"_" << ivm->index << ".dat";
         outf.open(outs.str().c_str());

         for (ivinstseq ivf=vshifts.begin();ivf!=vshifts.end();ivf++){
            if (ivf->species==i && ivf->motindex==ivm->index){
               outf << ivf->shift << endl;
               //cout << ivf->species << " " << ivf->motindex << " " << ivf->shift << endl;
            }
         }
         outf.close();
      }
   }
}

   void 
surrounding_score_stat(vseq & align, vmot & mots)
{
   //int length=700;
   int spe=0; // allows to do this statistics on any species (0=ref)
   int length=args_info.scanwidth_arg;
   vd vscoremean(2*length+1,0.);
   vd percent_gc(2*length+1,0.);
   vint totcounts(2*length+1,0);
   Motif mot;
   mot=mots[0];

   if (mots.size()>1) mot=mots[1];

   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
      if (ivs->instancescons.size()>0){
         for (ivvinstseq ivv=ivs->instancescons.begin();ivv!=ivs->instancescons.end();ivv++){
            int numspe=0;
            int isok=0;
            for (ivinstseq ivin=ivv->begin();ivin!=ivv->end();ivin++){
               if (ivin->species==spe){
                  isok=1;
                  break;
               }
               numspe++;
            }
            if (isok==0) continue;
            Instanceseq inst=ivv->at(numspe);
            int start;
            int stop;

            start=ivs->imaps[spe][inst.pos]-length;
            stop=ivs->imaps[spe][inst.pos]+length;
            int truestart=max(start,0);
            int truestop=min(stop,ivs->iseqs[spe].size()-width+1);   
            int i=truestart-start;
            //we want a total of length scores (include position 0 so as to count it twice like the others)
            //for (civint civ=ivs->iseqs[0].begin()+start;civ!=ivs->iseqs[0].begin()+start+length;civ++){}
            for (civint civ=ivs->iseqs[spe].begin()+truestart;civ!=ivs->iseqs[spe].begin()+truestop;civ++){
               //the Ns lead to artefacts (very low scores)
               int ismasked=0;
               for (int pos=0;pos<mot.motwidth;pos++){
                  if (*(civ+pos)==4) ismasked=1;
               }
               if (ismasked==1){
                  i++;
                  continue;
               }
               double score;
               score=scoref(civ,mot.matprec);
               if (scoref(civ,mot.matprecrevcomp)>score){
                  score=scoref(civ,mot.matprecrevcomp);
               }
               if (*civ==2||*civ==3) percent_gc[i]+=1.;
               //if (score>-10){
               vscoremean[i]+=score;
               totcounts[i]++;
               //}
               i++;
            }


            // FOR A SYMMETRIC CASE
            //            start=ivs->imaps[spe][inst.pos]-length;
            //           // if (start>0){}
            //           int truestart=max(start,0);
            //               int i=truestart-start;
            //               //we want a total of length scores (include position 0 so as to count it twice like the others)
            //               //for (civint civ=ivs->iseqs[0].begin()+start;civ!=ivs->iseqs[0].begin()+start+length;civ++){}
            //               for (civint civ=ivs->iseqs[spe].begin()+truestart;civ!=ivs->iseqs[spe].begin()+ivs->imaps[spe][inst.pos];civ++){
            //                  //the Ns lead to artefacts (very low scores)
            //                  if (*civ==4) continue;
            //                  double score;
            //                  score=scoref(civ,mot.matprec);
            //                  if (scoref(civ,mot.matprecrevcomp)>score){
            //                     score=scoref(civ,mot.matprecrevcomp);
            //                  }
            //                  vscoremean[i]+=score;
            //                  if (*civ==2||*civ==3) percent_gc[i]+=1.;
            //                  totcounts[i]++;
            //                  i++;
            //               }
            //
            //            stop=ivs->imaps[spe][inst.pos]+length;
            //            //if (stop<ivs->stop-ivs->start-width+1){}
            //            int truestop=min(stop,ivs->iseqs[spe].size()-width+1);
            //               //int i=length;
            //               //for (civint civ=ivs->iseqs[0].begin()+ivs->imaps[0][inst.pos];civ!=ivs->iseqs[0].begin()+stop+1;civ++){}
            //               for (civint civ=ivs->iseqs[spe].begin()+ivs->imaps[spe][inst.pos];civ!=ivs->iseqs[spe].begin()+truestop;civ++){
            //                  if (*civ==4) continue;
            //                  double score;
            //                  score=scoref(civ,mot.matprec);
            //                  if (scoref(civ,mot.matprecrevcomp)>score){
            //                     score=scoref(civ,mot.matprecrevcomp);
            //                  }
            //                  vscoremean[i]+=score;
            //                  if (*civ==2||*civ==3) percent_gc[i]+=1.;
            //                  totcounts[i]++;
            //                  i++;
            //            }
         }
      }
   }

   ofstream result;
   system("if ! test -d surround;then mkdir surround;fi");
   result.open("surround/results.dat");
   result << "Pos" << "\t";
   result << "Score" << "\t";
   result << "\%GC" << "\t";
   result << "number_of_instances" << "\n";
   for (int i=0;i<2*length+1;i++){
      vscoremean[i]/=totcounts[i];
      percent_gc[i]/=totcounts[i];
   }
   int window=args_info.window_arg;
   for (int i=0;i<2*length+1-window;i++){
      double mean_score(0);
      double mean_gc(0);
      for (int j=0;j<window;j++){
         mean_score+=vscoremean[i+j];
         mean_gc+=percent_gc[i+j];
      }
      mean_score/=window;
      mean_gc/=window;
      result << i-length << " " << mean_score << " " << mean_gc << " " << totcounts[i] << endl;
   }
   result.close();


}

   void
motcorr(vseq & align)
{
   //FOR TWO MOTIFS
   int distmax=5e2;
   if (species==2){ 
      distmax=2000;
   }
   int indextocomp=1;
   if (nbmots_for_score==1) indextocomp=0;
   vint distint;
   vseq seqs;

   if(args_info.first_interval_given){
      for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
         if (ivs->instancescons.size()>1){
            for (ivvinstseq ivv=ivs->instancescons.begin();ivv!=ivs->instancescons.end()-1;ivv++){
               if ((*ivv)[0].motindex==0){
                  int mindist(1e6);
                  for (ivvinstseq ivv1=ivv+1;ivv1!=ivs->instancescons.end();ivv1++){
                     if ((*ivv1)[0].motindex==indextocomp){
                        int pos1=ivs->imaps[0][(*ivv1)[0].pos];
                        int pos=ivs->imaps[0][(*ivv)[0].pos];
                        int dist=fabs(pos-pos1);
                        if (dist<mindist && dist>width) mindist=dist;
                     }
                  }
                  if (mindist<distmax) distint.push_back(mindist);

               }
            }
         }
      }
   }
   else {
      int numint(0);
      for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
         if (ivs->instancescons.size()>1){
            for (ivvinstseq ivv=ivs->instancescons.begin();ivv!=ivs->instancescons.end();ivv++){
               if ((*ivv)[0].motindex==0){
                  int pos=ivs->imaps[0][(*ivv)[0].pos];
                  int stop=ivs->stop-ivs->start-distmax;
                  // !! Here we need sequences to be at least 4xdistmax (2 distmax of sequence, 2x distmax of flanking to avoid side effects) large!
                  if (pos>distmax && pos <stop){
                     numint++;
                     for (ivvinstseq ivv1=ivs->instancescons.begin();ivv1!=ivs->instancescons.end();ivv1++){
                        if ((*ivv1)[0].motindex==indextocomp && ivv1!=ivv){
                           int pos1=ivs->imaps[0][(*ivv1)[0].pos];
                           int dist=fabs(pos-pos1);
                           if (dist<distmax) {
                              distint.push_back(dist);
                              seqs.push_back(*ivs);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   system("if ! test -d motcorr;then mkdir motcorr;fi");
   ofstream corrf;
   ostringstream corrfs;
   corrfs << "motcorr/corr.dat";
   corrf.open(corrfs.str().c_str());

   int i(0);
   for (ivint iv=distint.begin();iv!=distint.end();iv++){
      corrf << *iv << "\t";   
      corrf << seqs[i].name << "\t";   
      corrf << chromfromint(seqs[i].chrom) << "\t";   
      corrf << seqs[i].start << "\t";   
      corrf << seqs[i].stop << "\n";   
      i++;
   }

   corrf.close();
}


   void
peakcorr(vcoord & vpeaks1,vcoord & vpeaks2)
{
   cout << "Loading chrom lengths..." << endl;
   lengthchrom=loadlengthchrom();

   system("if ! test -d data;then mkdir data;fi");
   ofstream outfraw("data/peakcorr_raw.dat");
   ofstream outfhisto("data/peakcorr_histo.dat");

   int distmax=1e4;
   int binsize=distmax/1e2;
   //   double logbinsize=log10(distmax)/10;
   int numbins=distmax/binsize;
   //   int numlogbins=floor(log10(distmax)/logbinsize);
   vd distbin(numbins,0.);//proportion of peaks per bin
   //   vd distlogbin(numlogbins,0.);//proportion of peaks per bin

   //   vint distpeaks;//all the peaks distances
   cout << "Computing correlation..." << endl;

   for (ivcoord ivc1=vpeaks1.begin();ivc1!=vpeaks1.end()-1;ivc1++){
      int pos1=ivc1->start;
      //      if (pos1<distmax || pos1>lengthchrom[ivc1->chrom]-distmax){
      //         continue;
      //      }
      for (ivcoord ivc2=vpeaks2.begin();ivc2!=vpeaks2.end();ivc2++){
         if (ivc1->start!=ivc2->start && ivc1->chrom==ivc2->chrom){
            int pos2=ivc2->start;
            int mindist=fabs(pos2-pos1);
            //for periodic BC
            //            int dist1=fabs(pos2-lengthchrom[ivc1->chrom]-pos1);
            //            int dist2=fabs(pos2+lengthchrom[ivc1->chrom]-pos1);
            //            int mindist=min(dist,min(dist1,dist2));
            if (mindist<distmax){
               //               distpeaks.push_back(dist);
               int index=mindist/binsize;
               distbin[index]+=1./vpeaks1.size();
               //               int logindex=floor(log10(dist)/logbinsize);
               //               distlogbin[logindex]+=1./vpeaks.size();
               outfraw << mindist << endl;
            }
         }
      }
   }

   for (int i=0;i<numbins;i++){
      outfhisto << binsize*(i+0.5) << "\t" << distbin[i] << endl;   
   }

}

   vseq
addscorepergene(vseq & vscore)
{
   vseq vuniq;
   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
      if (ivs->name.find('_')!=string::npos)
      {
         ivs->name.erase(ivs->name.begin()+ivs->name.find('_')-2,ivs->name.end());
      }
   }
   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
      int iscounted=0;
      if (vuniq.size()>0){
         for (ivseq ivs1=vuniq.begin();ivs1!=vuniq.end();ivs1++){
            if (ivs1->name==ivs->name){
               ivs1->score+=ivs->score;
               for (int i=0;i<ivs1->motis.size();i++){
                  ivs1->motis[i]+=ivs->motis[i];
               }
               iscounted=1;
            }
         }
      }
      if (iscounted==0){
         vuniq.push_back(*ivs);
      }
   }

   return vuniq;
}

   void
jaspar2pwm()
{

   ifstream fmotifs;
   fmotifs.open(args_info.jaspar2pwm_arg);
   ofstream out;
   out.open("jaspar2pwm-output-wname.dat");
   string dum;
   getline(fmotifs,dum);
   while (!fmotifs.eof()){
      Motif mot;
      mot.name=dum.substr(1);
      mot.id=mot.name.substr(0,mot.name.find(" "));
      while (mot.name.find(" ")!=string::npos){
         mot.name.replace(mot.name.find(" "),1,"_");
      }
      //cout << mot.id << " " << mot.name << endl;

      vd rowdum;
      vvd mat(4,rowdum);
      // read the matrix
      for (int j=0;j<4;j++){
         getline(fmotifs,dum);
         stringstream row(dum);
         row >> dum; // base
         row >> dum; // [
         if (dum.size()>1){ // in case first number is stuck to the [
            dum=dum.substr(1);
         } else {
            row >> dum; // freq
         }
         while (dum.find(']')==string::npos){ // until ]
            stringstream dumtonum(dum);
            double num;
            dumtonum >> num;
            mat[j].push_back(num);
            row >> dum;
         }
         //cout << mat[j] ;
      }

      //cout << "Computing pseudo-count..." << endl;
      mot.motwidth=mat[0].size();
      width=mot.motwidth;
      double sc=scorethr2;
      scorethr2=sc/10*width;
      compalpha();
      scorethr2=sc;

      vd rowzeros(4,0.);
      mot.matprec=vvd(mot.motwidth,rowzeros);
      for (int i=0;i<mot.motwidth;i++){
         double sum=0;
         sum+=mat[0][i];
         sum+=mat[1][i];
         sum+=mat[2][i];
         sum+=mat[3][i];
         sum+=2*alpha+2*beta;
         mot.matprec[i][0]=log((mat[0][i]+alpha)/sum/conca);
         mot.matprec[i][1]=log((mat[3][i]+alpha)/sum/conca);
         mot.matprec[i][2]=log((mat[1][i]+beta)/sum/concc);
         mot.matprec[i][3]=log((mat[2][i]+beta)/sum/concc);

         string letter="A";
         double max(mot.matprec[i][0]);
         if (mot.matprec[i][1]>max) letter="T";
         else if (mot.matprec[i][2]>max) letter="C";
         else if (mot.matprec[i][3]>max) letter="G";
         mot.bsinit.append(letter);
      }
      //cout << mot.matprec << endl;
      mot.displaywname(out);
      getline(fmotifs,dum);
   }
   out.close();
   fmotifs.close();

}

   void
furlong2pwm()
{


   ofstream out;
   out.open("pwms_final.dat");
   ifstream fmotifs;
   fmotifs.open(args_info.furlong2pwm_arg);
   string dum;
   unsigned int icount(1);
   while (!fmotifs.eof()){

      string name;
      getline(fmotifs,name);
      getline(fmotifs,dum);
      istringstream fwmot(dum);
      fwmot >> width;

      Motif mot;
      mot.motwidth=width;
      mot.name=name;

      getline(fmotifs,dum);
      vd dumd(4,0.0);
      vvd mat(mot.motwidth,dumd);
      mot.matprec=mat;
      for (int row=0;row<4;row++)
      {
         istringstream fmot(dum);
         fmot >> dum;
         for (int col=0;col<width;col++)
         {
            fmot >> mat[col][row];
         }
         getline(fmotifs,dum);
      }

      for (unsigned int i=0;i<mot.motwidth;i++){
         vd & col=mat[i];
         double sum=0;
         sum+=col[0];
         sum+=col[1];
         sum+=col[2];
         sum+=col[3];
         sum+=2*alpha+2*beta;
         double col0=log((col[0]+alpha)/(conca*sum));//A
         double col1=log((col[1]+beta)/(concc*sum));//C
         double col2=log((col[2]+beta)/(concc*sum));//G
         double col3=log((col[3]+alpha)/(conca*sum));//T
         //convert in A/T/C/G format.
         col[0]=col0;//floor(100*col0);//A
         col[1]=col3;//floor(100*col2);//T
         col[2]=col1;//floor(100*col3);//C
         col[3]=col2;//floor(100*col1);//G
         //Consensus sequence
         string letter="A";
         double max(col[0]);
         if (col[1]>max) letter="T";
         else if (col[2]>max) letter="C";
         else if (col[3]>max) letter="G";
         mot.bsinit.append(letter);
         mot.matprec[i]=col;
      }
      mot.displaywname(out);
      icount++;
   }
   fmotifs.close();
   out.close();
}
   
   void
transfac2pwm()
{
   ofstream out;
   out.open("pwms.dat");
   ifstream fmotifs;
   fmotifs.open(args_info.transfac2pwm_arg);
   string dum;
   string name;
   getline(fmotifs,name);
   while (!fmotifs.eof()){

      Motif mot;
      mot.name=name;
      getline(fmotifs,dum);
      vd col(4,0.0);
      vvd mat;
      mot.motwidth=0;
      while(dum.substr(0,8)!="numsites"){
         istringstream fmot(dum);
         fmot >> col[0];
         fmot >> col[1];
         fmot >> col[2];
         fmot >> col[3];
         double sum=col[1]+col[2]+col[3]+col[0];
         for (int i=0;i<4;i++) col[i]/=sum;
         mot.matfreq.push_back(col);
         mot.motwidth++;
         getline(fmotifs,dum);
      }
      width=mot.motwidth;
         
      istringstream fns(dum);
      fns >> dum;
      fns >> mot.nbmot;
      // for genome-wide data:
      if (mot.nbmot==0) mot.nbmot=1000;

      cout << mot.name << endl;

      if (mot.name.substr(0,1)=="V"){
         conca=0.268;
         concc=0.232;
         compalpha();
      }
      else if (mot.name.substr(0,1)=="I"){
         conca=0.3;
         concc=0.2;
         compalpha();
      }
      else {
         alpha=1;
         beta=1;
      }

      if (alpha<=0 || beta<=0){
         alpha=1;
         beta=1;
      }

      for (unsigned int i=0;i<mot.motwidth;i++){
            col=mot.matfreq[i];
            for (unsigned int j=0;j<4;j++){
               col[j]*=mot.nbmot;
            }
            double sum=mot.nbmot;
            sum+=2*alpha+2*beta;
            col[0]=log((col[0]+alpha)/(conca*sum));//A
            col[1]=log((col[1]+alpha)/(conca*sum));//T
            col[2]=log((col[2]+beta)/(concc*sum));//C
            col[3]=log((col[3]+beta)/(concc*sum));//G
            //Consensus sequence
            string letter="A";
            double max(col[0]);
            if (col[1]>max) letter="T";
            else if (col[2]>max) letter="C";
            else if (col[3]>max) letter="G";
            mot.bsinit.append(letter);
            mot.matprec.push_back(col);
      }

      mot.displaywname(out);
      getline(fmotifs,name);
         
   }
   fmotifs.close();
   out.close();
}

   void
motistoscore(Combination & combi,vseq & vscore)
{
   vint motis_curr=combi.motis;

   unsigned int numpos(0),numneg(0);
   unsigned int numtotpos(0),numtotneg(0);
   int state; //-1 neg, 0 nothing, +1 pos
   for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){

      if (ivs->sign==1){
         state=1;
         numtotpos++;
      } else {
         state=-1;
         numtotneg++;
      }

      for (int i=0;i<motis_curr.size();i++){
         if (motis_curr[i]>0){
            if (ivs->motis[i]<motis_curr[i]){ // get sequences that have AT LEAST the muber of motifs
               state=0;
               break;
            }
         }

      }

      if (state==1){
         numpos++;
         combi.posnames.push_back(ivs->name);
         combi.posseqs.push_back(*ivs);
      }
      else if (state==-1) numneg++;
   }

   combi.TP=numpos;
   combi.FP=numneg;
   combi.TPR=(double)numpos/numtotpos;
   combi.FPR=(double)numneg/numtotneg;
   if (combi.FP!=0){
      combi.score=combi.TPR/combi.FPR;
   } else {
      if (combi.TP!=0) combi.score=1e6;  
      else combi.score=0;
   }

   return;
}

//recursive function to enumarate all possible vint of integers
//here we have motis.size() motifs with each having a maximal number of motis_max[i] of sites on a sequence
//n is for the nth motif. It stops at motis.size().
   void
recursive_combination(int n,vint motis_max,vint& motis_curr,vseq & vscore,vcombi & vcomb)
{
   if (n==motis_max.size())
   {  // final combination
      Combination combitmp;
      combitmp.motis=motis_curr;
      motistoscore(combitmp,vscore);
      if (combitmp.TP>50 || combitmp.TPR>0.01) vcomb.push_back(combitmp);
   }
   else
   {     // adding terms to combination
      //for(int j=-1;j<=motis_max[n];j++){}
      for(int j=0;j<=motis_max[n];j++)
      {  // -1 for any number of motifs, 
         vint motitemp(motis_curr);
         if (n==0)
         {
            motitemp.clear();
            motitemp.push_back(j);
         }
         else
         {
            motitemp.push_back(j);
         }
         recursive_combination(n+1,motis_max,motitemp,vscore,vcomb);
      }
   }
   return;
}

   void
countnmers(ofstream & outf)
{
   
   unsigned int numseq=1;
   for (ivseq iseq=regints.begin();iseq!=regints.end();iseq++){
      Sequence & seq=*iseq;
      for (vint::const_iterator istr=seq.iseqs[0].begin();istr!=seq.iseqs[0].end()-width+1;istr++){
         vint word=vint(istr,istr+width);
         outf << numseq << "\t" << vinttostring(word) << endl;
         vint revword=reversecomp(word);
         outf << numseq << "\t" << vinttostring(revword) << endl;
      }
      numseq++;
   }
   return;
}


   void
countwords(Motif & mot)
{
   Motif motinit=mot;

   cout << "Finding instances on interest sequences..." << endl;
   mot.comprefinstances(regints,0);
   displaymat(mot.matfreq);
   
//   ofstream motmeldb("ref-PWM.dat");
//   mot.display(motmeldb);
//   motmeldb.close();

//   ofstream outf("sites.dat");
//   for (ivstring iv=mot.refinstances.seqs.begin();iv!=mot.refinstances.seqs.end();iv++){
//      outf << *iv << endl;
//   }
//   outf.close();

   cout << "Enumerating all possible sites above threshold and rank them..." << endl;
   vint site;
   vtfbs sites;
   // Here probability is with correlations (initial matrix)
   mot.enumeratewthres(mot.motwidth,site,sites);
   double totprob(0); 
   for (ivtfbs ivt=sites.begin();ivt!=sites.end();ivt++){
      totprob+=ivt->prob;
   }
   for (ivtfbs ivt=sites.begin();ivt!=sites.end();ivt++){
      ivt->prob/=totprob;
   }
   sort(sites.begin(),sites.end());
   cout << "There are " << sites.size() << " sites above threshold" << endl;


   cout << "Counting words in dataset..." << endl;
   for (ivvint ivv=mot.refinstances.iseqs.begin();ivv!=mot.refinstances.iseqs.end();ivv++){
      for (ivtfbs ivt=sites.begin();ivt!=sites.end();ivt++){
         if (*ivv==ivt->site) ivt->num++;
      }
   }

   unsigned int numobs=mot.refinstances.iseqs.size();

   // NOW WE DEFINE A REFINED MATRIX, THAT GET RIDS OF BASES CORRELATIONS
   mot.comprefmot();
   double prob=0,sigma=0,nsigma;

   cout << "binding_site" << "\t";
   cout << "score" << "\t";
   cout << "Ini_PWM" << "\t"; //expected from PWM
   cout << "Thr_PWM" << "\t"; // expected from correlation due to threshold
   //   cout << "Mix_PWM" << "\t"; // expected from mixture of PWMs
   cout << "obs_num"  << "\t"; // observed
   cout << "n-sigma"  << "\n"; // n-sigma
   
   ofstream outf("sites-statistics-obsVSpwm.dat");
   // obs VS pwm
   //
   vd true_prob;
   vd thr_prob;
   for (ivtfbs ivt=sites.begin();ivt!=sites.end();ivt++){
      cout << *ivt << "\t";

      prob=ivt->prob;
      sigma=sqrt(prob*(1.-prob)*numobs);
      cout << floor(numobs*prob+0.5) << ""; //expected from INTIAL PWM
      cout << "+-" << floor(sigma+0.5) << "\t"; // stdev

      prob=compprob(ivt->site,mot.matfreq);
      thr_prob.push_back(prob);
      sigma=sqrt(prob*(1.-prob)*numobs);
      cout << floor(numobs*prob+0.5) << ""; //expected from PWM W CORRELATION
      cout << "+-" << floor(sigma+0.5) << "\t"; // stdev
         
      double probpwm=prob;

      //sigma=sqrt(prob*(1.-prob)*numobs);
      //cout << floor(numobs*prob+0.5) << ""; //expected from PWM mixture
      //cout << "+-" << floor(sigma+0.5) << "\t"; // stdev

      cout << floor(ivt->num+0.5) << "\t"; // observed
      true_prob.push_back((double)ivt->num/numobs);
      nsigma=abs((ivt->num-numobs*prob)/sigma);
      if (nsigma<5){
         cout << "<" << ceil(abs((ivt->num-numobs*prob)/sigma)) << " sigma" << "\n";
      } else {
         cout << ">5 sigma" << "\n";
      } 
      
      if (ivt->num>0){
         outf << vinttostring(ivt->site) << "\t"; // SITE
         outf << (double)ivt->num/numobs << "\t"; // OBS
         outf << probpwm << "\n";// PWM
      }
   }
   outf.close();
   //cout << sites << endl;
   //
   double distkl=0;
   for (unsigned int count=0;count<true_prob.size();count++){
      distkl+=true_prob[count]*log(true_prob[count]/thr_prob[count]);
   }
   cout << "Dkl_thr: " << distkl << endl;

}
   
   void
countwords(vmot & mots,vd kprobs)
{
   Motif refmot=mots[0];

   cout << "Finding instances on interest sequences..." << endl;
   refmot.comprefinstances(regints,0);
   
   ofstream outfs("sites.dat");
   for (ivstring iv=refmot.refinstances.seqs.begin();iv!=refmot.refinstances.seqs.end();iv++){
      outfs << *iv << endl;
   }
   outfs.close();

   cout << "Enumerating all possible sites above threshold and rank them..." << endl;
   vint site;
   vtfbs sites;
   // Here probability is with correlations (initial matrix)
   refmot.enumeratewthres(refmot.motwidth,site,sites);
   cout << "There are " << sites.size() << " sites above threshold" << endl;


   cout << "Counting words in dataset..." << endl;
   for (ivvint ivv=refmot.refinstances.iseqs.begin();ivv!=refmot.refinstances.iseqs.end();ivv++){
      for (ivtfbs ivt=sites.begin();ivt!=sites.end();ivt++){
         if (*ivv==ivt->site) ivt->num++;
      }
   }
   vtfbs tempsites;
   for (ivtfbs ivt=sites.begin();ivt!=sites.end();ivt++){
      if (ivt->num>0) tempsites.push_back(*ivt);
   }
   sites=tempsites;
   double totprob(0); 
   // INITIAL PWM PROBABILITIES
   for (ivtfbs ivt=sites.begin();ivt!=sites.end();ivt++){
      totprob+=ivt->prob;
   }
   for (ivtfbs ivt=sites.begin();ivt!=sites.end();ivt++){
      ivt->prob/=totprob;
   }
   sort(sites.begin(),sites.end());

   unsigned int numobs=refmot.refinstances.iseqs.size();

   // NOW WE DEFINE A REFINED MATRIX, THAT GET RIDS OF BASES CORRELATIONS
   refmot.comprefmot();
   double prob=0,sigma=0,nsigma;

//   cout << "binding_site" << "\t";
//   cout << "score" << "\t";
//   cout << "Ini_PWM" << "\t"; //expected from PWM
//   cout << "Thr_PWM" << "\t"; // expected from correlation due to threshold
//   cout << "Mix_PWM" << "\t"; // expected from mixture of PWMs
//   cout << "obs_num"  << "\t"; // observed
//   cout << "n-sigma"  << "\n"; // n-sigma

   ofstream outf("sites-statistics.dat");
   outf << "ini" << "\t";
   outf << "thr" << "\t";
   outf << "mix" << "\t";
   outf << "obs" << "\n";
   
   vd ini_prob,true_prob,thr_prob,mix_prob;

   double totmixprob=0;
   double totthrprob=0;
   for (ivtfbs ivt=sites.begin();ivt!=sites.end();ivt++){
         prob=0;
         for (unsigned int j=0;j<kprobs.size();j++){
            double score=1;
            unsigned int pos=0;
            for (ivint iv=ivt->site.begin();iv!=ivt->site.end();iv++){
               const int base=*iv;
               if (base == 4){
                  score*=exp(-100);
               }
               else score*=motsdef[j+1].matfreq[pos][base];
               pos++;
            }
            prob+=kprobs[j]*score;
            //if (compprob(ivt->site,motsdef[j+1].matfreq)>prob) prob=compprob(ivt->site,motsdef[j+1].matfreq);
         }
         totmixprob+=prob;
         totthrprob+=compprob(ivt->site,refmot.matfreq);
   //   }
   }
   for (ivtfbs ivt=sites.begin();ivt!=sites.end();ivt++){
     // if (ivt->num>0){
         //cout << *ivt << "\t";

      // INITIAL PWM PROBABILITY
         prob=ivt->prob;
         ini_prob.push_back(prob);
         sigma=sqrt(prob*(1.-prob)*numobs);
         //cout << floor(numobs*prob+0.5) << ""; //expected from INTIAL PWM
         outf << floor(numobs*prob+0.5) << "\t"; //expected from INTIAL PWM
         //cout << "+-" << floor(sigma+0.5) << "\t"; // stdev

         // THR PWM PROBABILITY
         prob=compprob(ivt->site,refmot.matfreq)/totthrprob;
         thr_prob.push_back(prob);
         sigma=sqrt(prob*(1.-prob)*numobs);
         //cout << floor(numobs*prob+0.5) << ""; //expected from PWM W CORRELATION
         outf << floor(numobs*prob+0.5) << "\t"; //expected from PWM W CORRELATION
         //cout << "+-" << floor(sigma+0.5) << "\t"; // stdev


         // MIXTURE MODEL WITH K-MEANS PWMs 
         double pj;
         prob=0;
         for (unsigned int j=0;j<kprobs.size();j++){
            double score=1;
            unsigned int pos=0;
            for (ivint iv=ivt->site.begin();iv!=ivt->site.end();iv++){
               const int base=*iv;
               if (base == 4){
                  score*=exp(-100);
               }
               else score*=motsdef[j+1].matfreq[pos][base];
               pos++;
            }
            prob+=kprobs[j]*score;
            //if (compprob(ivt->site,motsdef[j+1].matfreq)>prob) prob=compprob(ivt->site,motsdef[j+1].matfreq);
         }
         prob/=totmixprob;

         // MIXTURE PWM PROBABILITY
         mix_prob.push_back(prob);
         sigma=sqrt(prob*(1.-prob)*numobs);
         //cout << floor(numobs*prob+0.5) << ""; //expected from PWM mixture
         outf << floor(numobs*prob+0.5) << "\t"; //expected from PWM mixture
         //cout << "+-" << floor(sigma+0.5) << "\t"; // stdev

         // TRUE COUNT PROBABILITY
         //cout << floor(ivt->num+0.5) << "\t"; // observed
         outf << floor(ivt->num+0.5) << "\n"; // observed
         true_prob.push_back((double)ivt->num/numobs);
         nsigma=abs((ivt->num-numobs*prob)/sigma);
         if (nsigma<5){
            //cout << "<" << ceil(abs((ivt->num-numobs*prob)/sigma)) << " sigma" << "\n";
         } else {
            //cout << ">5 sigma" << "\n";
         }
   }
   //cout << accumulate(true_prob.begin(),true_prob.end(),0.) << endl;
   //cout << accumulate(ini_prob.begin(),ini_prob.end(),0.) << endl;
   //cout << accumulate(thr_prob.begin(),thr_prob.end(),0.) << endl;
   //cout << accumulate(mix_prob.begin(),mix_prob.end(),0.) << endl;
   //cout << sites << endl;
   double distkl=0;
   for (unsigned int count=0;count<true_prob.size();count++){
      distkl+=ini_prob[count]*log(ini_prob[count]/true_prob[count]);
   }
   cout << "Dkl_ini: " << distkl/log(2) << endl;
   distkl=0;
   for (unsigned int count=0;count<true_prob.size();count++){
      distkl+=thr_prob[count]*log(thr_prob[count]/true_prob[count]);
   }
   cout << "Dkl_thr: " << distkl/log(2) << endl;
   distkl=0;
   for (unsigned int count=0;count<true_prob.size();count++){
      distkl+=mix_prob[count]*log(mix_prob[count]/true_prob[count]);
   }
   cout << "Dkl_mix: " << distkl/log(2) << endl;

}
   
   void
countwords_mixture(vmot & mots,vd kprobs)
{

   cout << "Finding instances on interest sequences..." << endl;
   Motif refmot=mots[0];
   refmot.comprefinstances(regints,0);
   refmot.comprefmot();
   Motif motini=mots[0];
   vmot motmix=mots;
   motmix.erase(motmix.begin());
   
   // NOW WE DEFINE A REFINED MATRIX, THAT GET RIDS OF BASES CORRELATIONS
   cout << "Enumerating all possible sites and computing probs..." << endl;
   vint site;
   vtfbs sites;
   ofstream outf("sites-statistics.dat");
   outf << "site" << "\t";
   outf << "obs" << "\t";
   outf << "ini" << "\t";
   outf << "thr" << "\t";
   outf << "mix" << "\n";
   enumerateallsitesncomp(refmot.motwidth,site,sites,motini,refmot,motmix,kprobs,outf);

   return;
}
   
   void
countwords_numk(Motif & motini,unsigned int numk)
{
   vvd vkprobs;
   vvmot vmotmix;
   for (unsigned int k=2;k<=numk;k++){
      vmot motmix;
      stringstream ssm;
      ssm << "kmeans/k_" << k << "/matrices.dat";
      loadmotswnames(ssm.str().c_str(),motmix);
      motmix.erase(motmix.begin());
      vmotmix.push_back(motmix);

      stringstream ssk;
      ssk << "kmeans/k_" << k << "/kprobs.dat";
      vd kprobs;
      ifstream inf;
      inf.open(ssk.str().c_str());
      back_insert_iterator<vd> dest(kprobs);
      copy(iidouble(inf),iidouble(),dest);
      inf.close();
      vkprobs.push_back(kprobs);
   }

   cout << "Finding instances on interest sequences..." << endl;
   Motif refmot=motini;
   refmot.comprefinstances(regints,0);
   refmot.comprefmot();
   
   // NOW WE DEFINE A REFINED MATRIX, THAT GET RIDS OF BASES CORRELATIONS
   cout << "Enumerating all possible sites and computing probs..." << endl;
   vint site;
   vtfbs sites;
   ofstream outf("sites-statistics.dat");
   outf << "site" << "\t";
   outf << "obs" << "\t";
   outf << "ini" << "\t";
   outf << "thr" << "\t";
   outf << "mix_k=2.." << numk << "\n";
   enumeratesitesforkmeans(refmot.motwidth,site,sites,motini,refmot,vmotmix,vkprobs,outf);

   return;
}
   
   void
countwords_compare_pseudocount(Motif & motini)
{
   int ks[6]={2,3,4,5,10,30};
   //int ks[]={2,3,4,5};
   vvmot vmotmix_nops;
   vvd vkprobs_nops;
   vvmot vmotmix_t10;
   vvd vkprobs_t10;
   vvmot vmotmix_t6;
   vvd vkprobs_t6;

   for (unsigned int i=0;i<6;i++){
      int k=ks[i];

      cout << k << endl;
      vmot motmix_nops;
      stringstream ssm;
      ssm << "kmeans/k_" << k << "/nops/matrices.dat";
      loadmotswnames(ssm.str().c_str(),motmix_nops);
      motmix_nops.erase(motmix_nops.begin());
      vmotmix_nops.push_back(motmix_nops);

      vd kprobs_nops;
      stringstream ssk;
      ssk << "kmeans/k_" << k << "/nops/kprobs.dat";
      ifstream inf;
      inf.open(ssk.str().c_str());
      back_insert_iterator<vd> dest(kprobs_nops);
      copy(iidouble(inf),iidouble(),dest);
      inf.close();
      vkprobs_nops.push_back(kprobs_nops);
      
      vmot motmix_t10;
      stringstream ssm1;
      ssm1 << "kmeans/k_" << k << "/t10/matrices.dat";
      loadmotswnames(ssm1.str().c_str(),motmix_t10);
      motmix_t10.erase(motmix_t10.begin());
      vmotmix_t10.push_back(motmix_t10);

      vd kprobs_t10;
      stringstream ssk1;
      ssk1 << "kmeans/k_" << k << "/t10/kprobs.dat";
      inf.open(ssk1.str().c_str());
      back_insert_iterator<vd> dest1(kprobs_t10);
      copy(iidouble(inf),iidouble(),dest1);
      inf.close();
      vkprobs_t10.push_back(kprobs_t10);
      
      vmot motmix_t6;
      stringstream ssm2;
      ssm2 << "kmeans/k_" << k << "/t6/matrices.dat";
      loadmotswnames(ssm2.str().c_str(),motmix_t6);
      motmix_t6.erase(motmix_t6.begin());
      vmotmix_t6.push_back(motmix_t6);

      vd kprobs_t6;
      stringstream ssk2;
      ssk2 << "kmeans/k_" << k << "/t6/kprobs.dat";
      inf.open(ssk2.str().c_str());
      back_insert_iterator<vd> dest2(kprobs_t6);
      copy(iidouble(inf),iidouble(),dest2);
      inf.close();
      vkprobs_t6.push_back(kprobs_t6);
   }

   cout << "Finding instances on interest sequences..." << endl;
   Motif refmot=motini;
   refmot.comprefinstances(regints,0);
   refmot.comprefmot();

   // NOW WE DEFINE A REFINED MATRIX, THAT GET RIDS OF BASES CORRELATIONS
   cout << "Enumerating all possible sites and computing probs..." << endl;
   vint site;
   vtfbs sites;
   ofstream outf("sites-statistics-t6.dat");
   outf << "site" << "\t";
   outf << "obs" << "\t";
   outf << "ini" << "\t";
   outf << "thr" << "\t";
   outf << "mix_k" << "\n";
   enumeratesitesforkmeans(refmot.motwidth,site,sites,motini,refmot,vmotmix_t6,vkprobs_t6,outf);
   outf.close();
   outf.open("sites-statistics-t10.dat");
   outf << "site" << "\t";
   outf << "obs" << "\t";
   outf << "ini" << "\t";
   outf << "thr" << "\t";
   outf << "mix_k" << "\n";
   enumeratesitesforkmeans(refmot.motwidth,site,sites,motini,refmot,vmotmix_t10,vkprobs_t10,outf);
   outf.close();
   outf.open("sites-statistics-nops.dat");
   outf << "site" << "\t";
   outf << "obs" << "\t";
   outf << "ini" << "\t";
   outf << "thr" << "\t";
   outf << "mix_k" << "\n";
   enumeratesitesforkmeans(refmot.motwidth,site,sites,motini,refmot,vmotmix_nops,vkprobs_nops,outf);
   outf.close();

   return;
}
         
void
calcoptthr(Motif & mot,double fdrthr,ofstream & outf)
{
   double maxinfo(0);
   int j(0);
   for (ivvd ivv=mot.matprec.begin();ivv!=mot.matprec.end();ivv++){
      int i(0);
      double maxcol(-10);
      for (ivd iv=ivv->begin();iv!=ivv->end();iv++){
         if (*iv>maxcol) maxcol=*iv;
         i++;
      }
      j++;
      maxinfo+=maxcol;
   }

   mot.motscorethr2=maxinfo;
   mot.motscorethrcons=mot.motscorethr2;

   while (1){

      mot.matinit(mot.motscorethr2);
      double fdr=0;
      for (ivseq ivs=regints.begin();ivs!=regints.end();ivs++){
         if (ivs->sign==-1 && ivs->nbmot>0) fdr++;
      }
      
      if (fdr/sizeneg<fdrthr){
         mot.motscorethr2-=0.1;
         mot.motscorethrcons=mot.motscorethr2;
         continue;
      }
      
      double tpr=0;
      for (ivseq ivs=regints.begin();ivs!=regints.end();ivs++){
         if (ivs->sign==1 && ivs->nbmot>0) tpr++;
      }

      outf << mot.name << " " << mot.motscorethr2 << " " << tpr/sizepos << " " << fdr/sizeneg << endl;
      
      break;
   }

   return;

}

   void
extractcoordstofasta(ifstream & inf)
{
   //   cout << "Loading Chromosomes..." << endl;
   //   loadchroms();
   //   
   //   cout << "Loading coords..." << endl;
   //   vcoord vcds;
   //   vcds=loadcoordconserv(inf);
   vseq vs;
   vs=loadseqsnotmasked();

   cout << "Writing Sequences..." << endl;
   ofstream outf("seqs.dat");
   int length;
   //for (ivcoord ivc=vcds.begin();ivc!=vcds.end();ivc++){}
   for (ivseq ivs=vs.begin();ivs!=vs.end();ivs++){

      //length=ivc->stop-ivc->start;
      //string seq=chromints[ivc->chrom].seq.substr(ivc->start+1,length);

      outf << ">";
      //      outf << ivc->name << "\t";
      //      outf << chromfromint(ivc->chrom) << "\t";
      //      outf << ivc->start << "\t";
      //      outf << ivc->stop << "\n";
      //      outf << seq << "\n";
      outf << ivs->name << "\t";
      outf << chromfromint(ivs->chrom) << "\t";
      outf << ivs->start << "\t";
      outf << ivs->stop << "\n";
      outf << ivs->seqs[0] << "\n";
      exit(9);
   }

   outf.close();
}

void
scangen2cons(){

   ifstream inf;
   inf.open(args_info.coord_file_arg);

   ifstream align;
   if (species==1) align.open("/home/santolin/these/files/droso/align/all/align.dat");
   else if (species==2) align.open("/home/santolin/these/files/mus/epo/align.dat");
   alignscoord=loadcoordconserv(align);
   align.close();

   string rawline,sdum;
   int idum;
   getline(inf,rawline);
   Coordinate coord;
   Sequence seq;

   ofstream outf("result-cons.dat");

   while (!inf.eof()){

      stringstream line(rawline);
      line >> sdum;
      line >> sdum;
      int pos1,pos2;
      pos1=sdum.find(":");
      pos2=sdum.find("..");
      coord.chrom=intfromchrom(sdum.substr(0,pos1));
      coord.start=atoi(sdum.substr(pos1+1,pos2-pos1-1).c_str());
      coord.stop=atoi(sdum.substr(pos2+2).c_str());
      line >> coord.name;

      seq=coordtoseq(coord);
      if (!seq.species[0] || seq.nbtb==0){
         getline(inf,rawline);
         continue;
      }

      regints.clear();
      regints.push_back(seq);

      // HERE FOR DROSO, with MOT1=ovo stringent and MOT2=ovo rouge
      Motif mot1(motsdef[0]),mot2(motsdef[1]);
      mot1.matinit(mot1.motscorethr2);
      mot2.matinit(mot2.motscorethr2);
      if (mot1.seqs.size()>0 || mot2.seqs.size()>0){
         cout << rawline << endl;
         outf << rawline << endl;
      }
      getline(inf,rawline);
   }

   outf.close();

   return;
}


   void
weightmots()
{
   vseq seqs;
   cout << "Loading pos/neg sequences..." << endl;
   loadseqsforscore(seqs);

   for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      cout << ivm->name << endl;
      ivm->nbmot=0;
      width=ivm->motwidth;
      scorethr2=ivm->motscorethr2;
      ivm->calclambdaposneg(seqs);
      double scoresave=ivm->motscorethr2;
      while (ivm->lambda==0 || ivm->lambdatrain==0){
         ivm->motscorethr2=ivm->motscorethr2-0.5;
         ivm->calclambdaposneg(seqs);
         if (ivm->motscorethr2<0){
            ivm->motscorethr2=scoresave;
            ivm->calclambdaposneg(seqs);
            if (ivm->lambda==0 && ivm->lambdatrain>0){ivm->lambda=exp(1.)*ivm->lambdatrain;}
            else if (ivm->lambda>0 && ivm->lambdatrain==0){ivm->lambdatrain=exp(-1.)*ivm->lambda;}
            else { ivm->lambdatrain=1.e-4;ivm->lambda=ivm->lambdatrain/exp(1.);}
            break;
         }
      }
      cout << "Motif " << ivm->index+1 << " score " << ivm->motscorethr2 <<  " lb " << ivm->lambda << " lt " << ivm->lambdatrain << " log(lt/lb) " << log(ivm->lambdatrain/ivm->lambda) << endl;
      ivm->motscorethr2=scoresave;
   }
   return;
}

   void
args_init()
{
   width=args_info.width_arg;
   scorethr2=args_info.threshold_arg;
   species=args_info.species_arg;
   if (species==1){
      conca=0.3; 
      concc=0.2;
      nbspecies=12;
      cout << "Species: droso" << endl;
   }
   else if (species==2){
      conca=0.263; 
      concc=0.5-conca;
      nbspecies=10;
      cout << "Species: mus" << endl;
   }
   conct=conca;
   concg=concc;
   evolutionary_model=args_info.evolutionary_model_arg;

   scorethrcons=scorethr2-1.0;
   scorethr=scorethr2-1.0;

   neighbext=args_info.neighbext_arg;
   nbmots_for_score=args_info.nbmots_arg;
   scanwidth=args_info.scanwidth_arg;

   if (species==1) nbchrom=6;
   else if (species==2) nbchrom=21;

}

   void
args_init_scangen()
{
   scanwidth=args_info.scanwidth_arg;
   scanstep=args_info.scanstep_arg;

   nbmots_for_score=args_info.nbmots_arg;
}

   int
main(int argc, char** argv)
{
   if ( cmdline_parser(argc,argv, & args_info)!=0) 
      exit(1);

   args_init();

   rnginit();

   //                        compalpha();
   //                        exit(9);
   //   loadmots(args_info.motifs_arg,motsdef);
   //   inittreedist();

   if (args_info.scangen_given){

      compalpha();
      args_init_scangen();
      cout << "Loading Motifs" << endl;
      if (args_info.mots_w_names_given){
         loadmotswnames(args_info.mots_w_names_arg,motsdef);
      } 
      else{
         loadmots(args_info.motifs_arg,motsdef); 
      }
      cout << "Loaded " << motsdef.size() << " motifs." << endl;
      cout << "Nb mots for score: " << nbmots_for_score  << endl;

      if (args_info.opt_thr_file_given){
         ifstream opt;
         opt.open(args_info.opt_thr_file_arg);
         string dum;
         getline(opt,dum);
         int index=0;
         while (getline(opt,dum)){
            istringstream line(dum);
            line >> dum;
            double thr;
            line >> thr;
            motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
            motsdef[index].motscorethr=motsdef[index].motwidth*(thr-1)/10;
            motsdef[index].motscorethrcons=motsdef[index].motwidth*(thr-1)/10;
            index++;
         }
         opt.close();
      }
      else if (args_info.mean_info_given){
         cout << "Setting threshold to mean info..." << endl;
         for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
            ivm->setscorethr2meaninfo();
            cout << ivm->name << " mean thr = " << ivm->motscorethr2 << endl;
         }
      }
      else {
         cout << "Setting threshold to given value..." << endl;
         cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;
      }
      if (args_info.weightmots_given){
         cout << "Weighting Motifs" << endl;
         weightmots();
      }
      //      cout << "Setting threshold to mean info..." << endl;
      //      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      //         ivm->setscorethr2meaninfo();
      //         ivm->id.append("_");
      //         ivm->id.append(ivm->name);
      //         ivm->name=ivm->id;
      //      }

      if (args_info.phenotype_given){

         //cout << "Loading chroms" << "\n";
         //loadchroms();

         cout << "Loading phenotypes" << endl;
         loadannots();

         cout << "Loading instances" << endl;
         initgroupedinst();

         ifstream potregs;
         //if (species==1) potregs.open("/home/rouault/these/sequence/genomes/seqs.dat");
         //if (species==1) potregs.open("/home/santolin/these/files/droso/align/all/align-masked-files.dat");
         if (species==1) potregs.open("/home/santolin/these/files/droso/align/all/align-files.dat");
         else if (species==2) potregs.open("/home/santolin/these/files/mus/epo/align-files.dat");
         vstring regs;
         back_insert_iterator<vstring> dest(regs);
         copy(iisstring(potregs),iisstring(),dest);
         potregs.close();

         cout << "Scanning seqs" << endl;
         scanseqs(regs);

      } 
      else if (args_info.coord_file_given){

         cout << "Loading chromosomes" << "\n";
         loadchroms();

         cout << "Loading coordinates" << "\n";
         vcoord coords;
         coords=loadcoords();

         cout << "Loading instances" << endl;
         initgroupedinst();

         cout << "Loading aligments" << endl;
         ifstream align;
         if (species==1) align.open("/home/santolin/these/files/droso/align/all/align.dat");
         else if (species==2) align.open("/home/santolin/these/files/mus/epo/align.dat");
         alignscoord=loadcoordconserv(align);
         align.close();

         vstring regs;
         for (ivcoord ivc=alignscoord.begin();ivc!=alignscoord.end();ivc++){
            for (ivcoord ivc1=coords.begin();ivc1!=coords.end();ivc1++){
               if (ivc->chrom != ivc1->chrom || ivc->stop <= ivc1->start || ivc->start >= ivc1->stop) continue;
               else regs.push_back(ivc->name);
            }
         }

         for (ivstring iv=regs.begin();iv!=regs.end();iv++){
            cout << *iv << endl;
         }


         cout << "Scanning seqs" << endl;
         scanseqs(regs);
      }
   } 
   else if (args_info.scangen2cons_given){
      cout << "Loading Motifs" << endl;
      if (args_info.mots_w_names_given) loadmotswnames(args_info.mots_w_names_arg,motsdef); 
      else loadmots(args_info.motifs_arg,motsdef); 

      // !! HERE SCORETHRCONS IS SET TO SCORETHR2
      if (args_info.opt_thr_file_given){
         ifstream opt;
         opt.open(args_info.opt_thr_file_arg);
         string dum;
         getline(opt,dum);
         int index=0;
         while (getline(opt,dum)){
            istringstream line(dum);
            line >> dum;
            double thr;
            line >> thr;
            motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
            motsdef[index].motscorethr=motsdef[index].motwidth*(thr-1)/10;
            motsdef[index].motscorethrcons=motsdef[index].motwidth*(thr-1)/10;
            motsdef[index].motscorethrcons=motsdef[index].motwidth*thr/10;
            index++;
         }
         opt.close();
      }
      else if (args_info.mean_info_given){
         cout << "Setting threshold to mean info..." << endl;
         for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
            ivm->setscorethr2meaninfo();
            cout << ivm->name << " mean thr = " << ivm->motscorethr2 << endl;
            ivm->motscorethrcons=ivm->motscorethr2;
         }
      }
      else {
         cout << "Setting threshold to given value..." << endl;
         scorethrcons=scorethr2;
         cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;
      }

      cout << "Keeping only conserved CRMs..." << endl;
      scangen2cons();

   }
   else if (args_info.scanmots_given){
      cout << "Loading Motifs" << endl;
      if (args_info.mots_w_names_given) loadmotswnames(args_info.mots_w_names_arg,motsdef); 
      else loadmots(args_info.motifs_arg,motsdef); 
      cout << "Loaded " << motsdef.size() << " motifs." << endl;

      // !! HERE SCORETHRCONS IS SET TO SCORETHR2
      if (args_info.opt_thr_file_given){
         ifstream opt;
         opt.open(args_info.opt_thr_file_arg);
         string dum;
         getline(opt,dum);
         int index=0;
         while (getline(opt,dum)){
            istringstream line(dum);
            line >> dum;
            double thr;
            line >> thr;
            motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
            motsdef[index].motscorethr=motsdef[index].motwidth*(thr-1)/10;
            motsdef[index].motscorethrcons=motsdef[index].motwidth*(thr-1)/10;
            motsdef[index].motscorethrcons=motsdef[index].motwidth*thr/10;
            index++;
         }
         opt.close();
      }
      else if (args_info.mean_info_given){
         cout << "Setting threshold to mean info..." << endl;
         for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
            ivm->setscorethr2meaninfo();
            cout << ivm->name << " mean thr = " << ivm->motscorethr2 << endl;
            ivm->motscorethrcons=ivm->motscorethr2;
         }
      }
      else {
         cout << "Setting threshold to given value..." << endl;
         scorethrcons=scorethr2;
         cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;
      }

      cout << "Finding conserved instances..." << endl;
      scanmots();

   }
   else if (args_info.display_given){


      cout << "Loading (unmasked) alignments " << endl;
      vseq align;
      align=loadseqsnotmasked();
      cout << "Nb sequences to scan: " << align.size() << endl;

      if (args_info.jaspardb_given){

         cout << "Loading Jaspar Motifs..." << endl;
         loadjaspardb(motsdef);
         cout << "Setting threshold to mean info..." << endl;
         for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
            ivm->setscorethr2meaninfo();
            ivm->id.append("_");
            ivm->id.append(ivm->name);
            ivm->name=ivm->id;
         }
         nbmots_for_score=motsdef.size();
         cout << "Loaded " << motsdef.size() << " motifs." << endl;

         cout << "Scanning sequences for instances..." << endl;
         scanseqsforinstances(align,motsdef);

         cout << "Defining conserved instances..." << endl;

         system("if ! test -d TFBS-pos;then mkdir TFBS-pos;fi;");      
         string folder="TFBS-pos/";
         ofstream outf;
         cout << "Writing TFBS-pos files..." << endl;
         for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
            ivs->instances2instancescons();
            if (ivs->instancescons.size()){
               stringstream file;
               file << folder;
               file << ivs->name << "_";
               file << chromfromint(ivs->chrom) << "_";
               file << ivs->start << "_";
               file << ivs->stop << ".dat";
               outf.open(file.str().c_str());
               outf << ivs->instancescons;
               outf.close();
            }
         }
         stringstream file;
         file << folder;
         file << "all-motifs.dat";
         outf.open(file.str().c_str());
         for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
            bool show(0);
            for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
               if (ivs->instancescons.size()){
                  for (ivvinstseq ivv=ivs->instancescons.begin();ivv!=ivs->instancescons.end();ivv++){
                     Instanceseq ist=ivv->at(0);
                     if (ist.motindex!=ivm->index) continue;
                     if (!show){
                        outf << ivm->name << "\t";
                        show=1;
                     }
                     outf << "chr" << ivs->chrom << ":";
                     int pos=ivs->start+ist.truepos;
                     if (ist.sens==1) outf << pos << "-" << pos+ivm->motwidth-1;
                     else if (ist.sens==-1) outf << pos-ivm->motwidth+1 << "-" << pos;
                     outf << " (" <<ivs->name << ",";
                     outf <<  ist.seq << ",";
                     outf << ist.score << ")";
                     //if (iv!=ivs->instancescons[0].end()-1) 
                     outf <<",\t";
                  }
               }
            }
            if (show==1) outf << endl;
         }
         outf.close();
      }
      else{

         cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;
         cout << "Loading Motifs" << endl;
         if (args_info.mots_w_names_given){
            loadmotswnames(args_info.mots_w_names_arg,motsdef);
         } else{
            loadmots(args_info.motifs_arg,motsdef); 
         }

         if (nbmots_for_score<motsdef.size()) motsdef.erase(motsdef.begin()+nbmots_for_score,motsdef.end());
         cout << "Loaded " << motsdef.size() << " motifs." << endl;
         for (ivmot iv=motsdef.begin();iv!=motsdef.end();iv++){
            iv->motscorethrcons=iv->motscorethr2;
         }
         if (args_info.opt_thr_file_given){
            ifstream opt;
            opt.open(args_info.opt_thr_file_arg);
            string dum;
            getline(opt,dum);
            int index=0;
            while (getline(opt,dum)){
               istringstream line(dum);
               line >> dum;
               double thr;
               line >> thr;
               motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
               motsdef[index].motscorethr=motsdef[index].motwidth*thr/10;
               motsdef[index].motscorethrcons=motsdef[index].motwidth*thr/10;
               index++;
            }
            opt.close();
         }


         if (args_info.disp_tex_given){

            ofstream outf;
            cout << "Creating fasta/tex files... " << endl;
            if (args_info.wgaps_given){
               disptexwgaps(align);
            } else{
               disptex(align);
            }
         }
         if (args_info.disp_svg_given){
            scanseqsforsvg(align);
         }
      }
   }
   else if (args_info.disp_mean_info_given){
      cout << "Loading Motif(s)..." << endl;
      if (args_info.mots_w_names_given){
         loadmotswnames(args_info.mots_w_names_arg,motsdef);
      }
      else {
         loadmots(args_info.motifs_arg,motsdef);
      }
      cout << "Setting scorethr2 to mean info..." << endl;
      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
         ivm->setscorethr2meaninfo();
         cout << ivm->name << " mean info = " << ivm->motscorethr2 << endl;
      }
   }
   else if (args_info.dispscore_given){
      cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;

      system("if ! test -d score;then mkdir score;fi;");      

      cout << "Loading Motifs..." << endl;
      if (args_info.mots_w_names_given) loadmotswnames(args_info.mots_w_names_arg,motsdef);
      else loadmots(args_info.motifs_arg,motsdef);
      cout << "Loaded " << motsdef.size() << " motifs." << endl;
      cout << "Nb mots for score: " << nbmots_for_score  << endl;

      if (args_info.opt_thr_file_given){
         ifstream opt;
         opt.open(args_info.opt_thr_file_arg);
         string dum;
         getline(opt,dum);
         int index=0;
         while (getline(opt,dum)){
            istringstream line(dum);
            line >> dum;
            double thr;
            line >> thr;
            motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
            motsdef[index].motscorethr=motsdef[index].motwidth*thr/10;
            motsdef[index].motscorethrcons=motsdef[index].motwidth*thr/10;
            index++;
         }
         opt.close();
      }
      if (args_info.weightmots_given){
         cout << "Weighting Motifs" << endl;
         weightmots();
      }

      vseq vscore;
      cout << "Loading pos/neg sequences..." << endl;
      loadseqsforscore(vscore);
      for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
         if (ivs->name.find('_')!=string::npos)
         {
            (*ivs).name.insert((*ivs).name.find('_'),"\\");
         }
      }

      system("if ! test -d score;then mkdir score;fi;");      
      system("if ! test -d score/dispscore;then mkdir score/dispscore;fi;");      
      calcscore(motsdef,vscore);
      sort(vscore.begin(),vscore.end());
      dispscore(vscore,motsdef);
      //for (ivmot iv=motsdef.begin();iv!=motsdef.end();iv++){
      //   dispscore(vscore,*iv);
      //}


   }
   else if (args_info.motgen_given){

      scorethr=scorethr2-2.*(double)width/10;//(double)width/10*6.5;
      scorethrcons=scorethr2-1.*(double)width/10;//(double)width/10*6.5;
      compalpha();
      cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;
      cout << "alpha: " << alpha << endl;
      ofstream motmeldb;
      motmeldb.open("motmeldb.txt");

      cout << "Loading background set..." << endl;
      ifstream backreg;
      //      if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000-align.dat"); //regs2000.fa");
      //      else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs2000-align.dat");//regs2000.fa");
      //      regtests=loadsequencesconserv(backreg);//loadsequences(backreg);
      if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000.fa");
      else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs2000-masked.fa");
      regtests=loadsequences(backreg);
      backreg.close();
      cout << "Background set size : " << regtests.size() << endl;

      cout << "Loading training set..." << endl;
      regints=loadseqs();
      cout << "Training set size : " << regints.size() << endl;
      inittreedist();

      for (vseq::iterator iseq=regints.begin();iseq!=regints.end();iseq++){
         cout << (*iseq).name << endl;
         if (args_info.ref_given){
            refseqanalysis(*iseq,motmeldb);
         }
         else{
            seqanalysis(*iseq,motmeldb);
         }
         cout << endl;
      }
      motmeldb.close();
   }
   else if (args_info.motaffin_given){

      ofstream motmeldb;
      motmeldb.open("motmeldb.txt");

      cout << "Loading background set..." << endl;
      ifstream backreg;
      if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000.fa");
      else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs2000.fa");
      regtests=loadsequences(backreg);
      backreg.close();
      cout << "Background set size : " << regtests.size() << endl;

      cout << "Loading training set..." << endl;
      regints=loadseqs();
      cout << "Training set size : " << regints.size() << endl;

      inittreedist();

      vmot motinit;
      cout << "Loading Motifs..." << endl;
      loadmots(args_info.motifs_arg,motinit);

      cout << "Raffinement..." << endl;
      motaffin(motinit,motmeldb);
      motmeldb.close();

   }
   else if (args_info.motcorr_given){
      cout << "Loading Motifs" << endl;
      if (args_info.mots_w_names_given){
         loadmotswnames(args_info.mots_w_names_arg,motsdef);
      } 
      else{
         loadmots(args_info.motifs_arg,motsdef); 
      }
      cout << "Loaded " << motsdef.size() << " motifs." << endl;
      cout << "Nb mots for corr: " << nbmots_for_score  << endl;
      if (nbmots_for_score==1){
         motsdef.erase(motsdef.begin()+1,motsdef.end());
      }

      if (args_info.opt_thr_file_given){
         ifstream opt;
         opt.open(args_info.opt_thr_file_arg);
         string dum;
         getline(opt,dum);
         int index=0;
         while (getline(opt,dum)){
            istringstream line(dum);
            line >> dum;
            double thr;
            line >> thr;
            motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
            motsdef[index].motscorethr=motsdef[index].motwidth*(thr-1)/10;
            motsdef[index].motscorethrcons=motsdef[index].motwidth*(thr-1)/10;
            index++;
         }
         opt.close();
      }

      cout << "Loading alignments " << endl;
      vseq align;
      if (args_info.repeat_masker_given){
         cout << "Using repeat masker..." << endl;
         align=loadseqs();
      }
      else {
         cout << "No repeat masker..." << endl;
         align=loadseqsnotmasked();
      }
      cout << "Nb sequences to scan: " << align.size() << endl;

      //OLD VERSION FOR ENRICHMENT      
      //      cout << "Loading background set..." << endl;
      //      ifstream backreg;
      //      if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000-align-1000.dat"); 
      //      else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs100000-align-100.dat");//regs2000.fa");
      //      regtests=loadsequencesconserv(backreg);//loadsequences(backreg);
      //      cout << "Loaded " << regtests.size() << " sequences." << endl;
      //      backreg.close();

      cout << "Scanning sequences for instances..." << endl;
      //scanseqsforinstancesnmask(align,motsdef);
      // NO MASKING FOR THE CORRELATION
      scanseqsforinstances(align,motsdef);
      //      scanseqsforinstancesnmask(regtests,motsdef);

      cout << "Defining conserved instances..." << endl;
      for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
         ivs->instances2instancescons();
         //cout << ivs->name << "\n" << ivs->instancescons;
      }
      //      for (ivseq ivs=regtests.begin();ivs!=regtests.end();ivs++){
      //         ivs->instances2instancescons();
      //cout << ivs->name << "\n" << ivs->instancescons;
      //      }

      cout << "Computing correlations..." << endl;
      motcorr(align);
      //motcorr(align,regtests);

   }
   else if (args_info.pwmcorr_given){

      cout << "Loading Motif" << endl;
      if (args_info.mots_w_names_given){
         loadmotswnames(args_info.mots_w_names_arg,motsdef);
      } 
      else{
         loadmots(args_info.motifs_arg,motsdef); 
      }
      Motif mot=motsdef[0];
      if (args_info.mean_info_given){
         cout << "Setting scorethr2 to mean info..." << endl;
         mot.setscorethr2meaninfo();
      }
      mot.motscorethrcons=mot.motscorethr2;
      mot.motscorethr=mot.motscorethr2;
      scorethr2=mot.motscorethr2;
      cout << "scorethr2=" << scorethr2 << endl;
      scorethr=scorethr2;
      scorethrcons=scorethr2;
      width=mot.motwidth;
      compalpha();
      cout << "Loading alignments " << endl;
      regints=loadseqs();
      cout << "Nb sequences to scan: " << regints.size() << endl;

      cout << "Computing base correlations..." << endl;
      mot.correlations();

   }
   else if (args_info.peakcorr_given){

      vcoord vpeaks1,vpeaks2;

      cout << "Loading peaks 1 coords..." << endl;
      ifstream ifc;
      ifc.open(args_info.peaks1_arg);
      vpeaks1=loadcoordpeaks(ifc);
      ifc.close();
      cout << "Nb peaks 1 for corr: " << vpeaks1.size() << endl;

      cout << "Loading peaks 2 coords..." << endl;
      ifc.open(args_info.peaks2_arg);
      vpeaks2=loadcoordpeaks(ifc);
      ifc.close();
      cout << "Nb peaks 2 for corr: " << vpeaks2.size() << endl;

      cout << "Computing correlations..." << endl;
      peakcorr(vpeaks1,vpeaks2);

   }
   else if (args_info.roc_given){
      //      compalpha();
      cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;

      system("if ! test -d score;then mkdir score;fi;");      

      cout << "Loading Motifs..." << endl;
      if (args_info.mots_w_names_given) loadmotswnames(args_info.mots_w_names_arg,motsdef);
      else loadmots(args_info.motifs_arg,motsdef);
      cout << "Loaded " << motsdef.size() << " motifs." << endl;
      cout << "Nb mots for score: " << nbmots_for_score  << endl;

      if (args_info.opt_thr_file_given){
         ifstream opt;
         opt.open(args_info.opt_thr_file_arg);
         string dum;
         getline(opt,dum);
         int index=0;
         while (getline(opt,dum)){
            istringstream line(dum);
            line >> dum;
            double thr;
            line >> thr;
            motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
            motsdef[index].motscorethr=motsdef[index].motwidth*(thr-1)/10;
            motsdef[index].motscorethrcons=motsdef[index].motwidth*(thr-1)/10;
            index++;
         }
         opt.close();
      }
      else if (args_info.mean_info_given){
         cout << "Setting threshold to mean info..." << endl;
         for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
            ivm->setscorethr2meaninfo();
            ivm->motscorethrcons=ivm->motscorethr2;
            cout << ivm->name << " mean thr = " << ivm->motscorethr2 << endl;
         }
      }
      else {
         for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
            ivm->motscorethrcons=ivm->motscorethr2;
         }
      }

      cout << "Loading background set..." << endl;
      ifstream backreg;
      if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000-align-100.dat"); 
      else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs2000-align-100.dat");//regs2000.fa");
      regtests=loadsequencesconserv(backreg);//loadsequences(backreg);
      cout << "Loaded " << regtests.size() << " sequences." << endl;
      backreg.close();

      cout << "Computing lambda background..." << endl;
      for (ivmot im=motsdef.begin();im!=min(motsdef.end(),motsdef.begin()+nbmots_for_score);im++){
         width=im->bsinit.size();
         scorethr2=im->motscorethr2;
         scorethr=im->motscorethr;
         scorethrcons=im->motscorethrcons;
         im->calclambdaback();
      }

      vseq vscore;
      cout << "Loading pos/neg sequences..." << endl;
      loadseqsforscore(vscore);


      if (args_info.rocopti_given){
         //
         //         cout << "Loading training set (for lambda)..." << endl;
         //         regints=loadseqs();
         //         cout << "Training set size : " << regints.size() << endl;

         cout << "Optimizing motifs..." << endl;
         motopti(vscore);

      }
      if (args_info.combine_mot_given){

         system("if ! test -d score/roc-combine;then mkdir score/roc-combine;fi;");      
         cout << "Combining motifs..." << endl;
         vmot motsdeftmp(motsdef);
         ofstream auc;
         auc.open("score/roc-combine/auc.dat");
         for (ivmot ivm=motsdeftmp.begin();ivm!=motsdeftmp.begin()+nbmots_for_score;ivm++){
            for (ivmot ivm1=ivm;ivm1!=motsdeftmp.begin()+nbmots_for_score;ivm1++){
               ofstream outf;
               ofstream outf1;
               ostringstream outfs;
               ostringstream outfs1;
               outfs << "score/roc-combine/bestseqs-coords-";
               outfs1 << "score/roc-combine/bestseqs-";
               motsdef.clear();
               motsdef.push_back(*ivm);
               auc << ivm->index;
               outfs << ivm->index;
               outfs1 << ivm->index;
               cout << ivm->index;
               if (ivm1!=ivm) {
                  motsdef.push_back(*ivm1);
                  auc << "+" << ivm1->index;
                  outfs << "+" << ivm1->index;
                  outfs1 << "+" << ivm1->index;
                  cout << "+" << ivm1->index;
               }
               cout << endl;
               outfs << ".dat" ;
               outfs1 << ".dat" ;
               vseq vscorebest;
               vscorebest=findbestseqs(motsdef,vscore);
               sort(vscorebest.begin(),vscorebest.end());
               roccombine(vscorebest);
               outf.open(outfs.str().c_str());
               outf1.open(outfs1.str().c_str());
               vseq vscorepos;
               for (ivseq ivs=vscorebest.begin();ivs!=vscorebest.end();ivs++){
                  if (ivs->sign==1){
                     vscorepos.push_back(*ivs);
                  }
               }
               dispbestseqs(vscorepos,outf);
               outf1 << vscorepos << endl;
               cout << vscorepos << endl;
               outf.close();
               outf1.close();
               auc << "\t" << calcareaforscore(vscorebest) << endl;
            }
            cout << endl;
         }
         motsdef=motsdeftmp;
         auc.close();
      }
      else if(args_info.roc_on_sites_given){
         system("if ! test -d score/roc-on-sites;then mkdir score/roc-on-sites;fi;");      

         cout << "Computing ROCs..." << endl;
         roconsites(vscore,motsdef);

      }
      else if(args_info.best_pwm_given){

         cout << "Computing AUCs..." << endl;
         bestpwm(vscore,motsdef);

      }
      else if(args_info.info_given){

         cout << "Computing ROC infos..." << endl;
         rocinfo(vscore,motsdef);

      }
      else if (args_info.per_mot_given){

         cout << "Computing ROCs..." << endl;
         rocpermot(vscore,motsdef);
      }
      else if (args_info.per_mot_var_thr_given){

         cout << "Computing ROCs (we show nb of conserved motifs)..." << endl;
         rocpermotvarthr(vscore,motsdef);
      }
      else if (args_info.cumul_mot_var_thr_given){

         cout << "Computing ROCs (we show nb of conserved motifs)..." << endl;
         roccumulmotvarthr(vscore,motsdef);
      }
      else if (args_info.correlations_given){
         cout << "Computing correlations..." << endl;
         correlations(vscore,motsdef);
      }
      else if (args_info.cumul_mot_given){

         cout << "Computing ROCs (we show nb of conserved motifs)..." << endl;
         roccumulmot(vscore,motsdef);
      }
      else if (args_info.bestcombi_given){
         //system("if test -d score;then rm -r score;fi;");      

         cout << "Computing number of conserved binding sites for each motif..." << endl;
         int i=0;
         vint motis_max(nbmots_for_score,0);

         for (ivmot ivm=motsdef.begin();ivm!=min(motsdef.end(),motsdef.begin()+nbmots_for_score);ivm++){

            double maxinfo(0);
            double meaninfo(0);
            int j(0);

            if (!(args_info.threshold_given||args_info.opt_thr_file_given)){
               for (ivvd ivv=ivm->matprec.begin();ivv!=ivm->matprec.end();ivv++){
                  int i(0);
                  double maxcol(-10);
                  for (ivd iv=ivv->begin();iv!=ivv->end();iv++){
                     if (*iv>maxcol) maxcol=*iv;
                     meaninfo+=ivm->matfreq[j][i]*(*iv);
                     i++;
                  }
                  j++;
                  maxinfo+=maxcol;
               }

               ivm->motscorethr2=min(0.95*maxinfo,meaninfo);
               ivm->motscorethr=ivm->motwidth*(ivm->motscorethr2-1)/10;
               ivm->motscorethrcons=ivm->motwidth*(ivm->motscorethr2-1)/10;
            }
            scorethr2=ivm->motscorethr2;
            scorethr=ivm->motscorethr;
            scorethrcons=ivm->motscorethrcons;
            cout << ivm->name << "\t" << "at threshold " << scorethr2 << endl;


            for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
               ivs->motis[i]=ivm->nbmatchcons(*ivs);
               if (ivs->motis[i]>motis_max[i]){
                  motis_max[i]=min(ivs->motis[i],cutoff_for_combination);
               }
            }
            i++;


         }
         //      cout << vscore << endl;
         
         // UNCOMMENT AND COMMENT NBMATCHCONS UPSTREAM FOR BEST 1kb SEARCH
//         cout << "Keeping best1kb..." << endl;
//         vseq vscorebest;      
//         vscorebest=findbestseqs(motsdef,vscore);
//         sort(vscorebest.begin(),vscorebest.end());
//         vscore=vscorebest;

         for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
            i=0;
            for (ivmot ivm=motsdef.begin();ivm!=min(motsdef.end(),motsdef.begin()+nbmots_for_score);ivm++){
               if (ivs->motis[i]>motis_max[i]){
                  motis_max[i]=min(ivs->motis[i],cutoff_for_combination);
               }
               i++;
            }
         }
         //cout << vscorebest << endl;
         //cout << "sizepos=" <<sizepos << " sizeneg=" << sizeneg << endl;

         //cout << "Writing file for SVM..." << endl;
         //printforsvm(vscore);


         cout << "Computing combinations scores..." << endl;
         vcombi vcomb;
         vint motis_curr;
         recursive_combination(0,motis_max,motis_curr,vscore,vcomb);

         cout << "Sorting combinations..." << endl;
         sort(vcomb.begin(),vcomb.end());

         cout << "Writing results file..." << endl;
         ofstream outf;
         outf.open("score/combis.dat");
         for (ivmot ivm=motsdef.begin();ivm!=min(motsdef.end(),motsdef.begin()+nbmots_for_score);ivm++){
            outf  << ivm->name << "\t";
         }
         outf << "TPR" << "\t";
         outf << "TP" << "\t";
         outf << "FPR" << "\t";
         outf << "FP" << "\t";
         outf << "score" << "\n";
         outf << vcomb << endl;
         outf.close();

         //         for (ivstring ivs=vcomb[0].posnames.begin();ivs!=vcomb[0].posnames.end();ivs++){
         //            cout << *ivs << endl;
         //         }

         outf.open("score/combis-enhancers.dat");
         for (ivmot ivm=motsdef.begin();ivm!=min(motsdef.end(),motsdef.begin()+nbmots_for_score);ivm++){
            outf  << ivm->name << "\t";
         }
         outf << "TPR" << "\t";
         outf << "TP" << "\t";
         outf << "FPR" << "\t";
         outf << "FP" << "\t";
         outf << "score" << "\n";
         // not the last one, which gives all the sequences
         for (ivcombi ivc=vcomb.begin();ivc!=vcomb.end()-1;ivc++){
            outf << "========" << endl;
            outf << *ivc << endl;
            outf << "========" << endl;
            for (ivseq ivs=ivc->posseqs.begin();ivs!=ivc->posseqs.end();ivs++){
               outf << ivs->name << "\t";
               outf << chromfromint(ivs->chrom) << "\t";
               outf << ivs->start << "\t";
               outf << ivs->stop << "\n";
            }
         }
         outf.close();
      }
   }
   else if (args_info.pwm_given){

      Sequence bds;

      //extract the sites
      if (args_info.tfbs_given)
      {
         ifstream filesites(args_info.tfbs_arg);

         getsites(filesites,bds);// to get the width

         if (bds.iseqs[0].size()==0 && !args_info.width_given){
            cout << "Error: give a site sequence or precise a site width" << endl;
            exit(1);
         } else {
            width=bds.iseqs[0].size();
         }
      }
      else{
         if (!args_info.width_given){
            cout << " You have to precise a width" << endl;
            exit(1);
         }
         cout << "Width: " << width << endl;
      }

      scorethr2=width*scorethr2/10; //we set the score to 0.9 per base

      Motif mot;
      vd dum(4,0.0);
      mot.matrice=vvd(width,dum);

      //!!!!!!!!!!!!!!!!!!!!!! HERE I HAVE TO CHANGE THE FUNCTION !!!!!!!!!!!!!!!!!!!!

      if (args_info.mat_given){
         ifstream filemat(args_info.mat_arg);
         getmatrices(filemat,mot);
      }

      if (args_info.tfbs_given){
         // we count the bases in the sites
         countbases(mot,bds);
      }

      // we get the bayesian frequency count
      //   vvd matnorm=vvd(width,dum);
      //cout << "\n sites width : " << width << endl;
      if (args_info.nopseudocount_given){
         alpha=0;
         beta=0;
      }
      else{
         compalpha();
      }
      countfreq(mot.matrice);

      //we convert to pwm
      freqtolog(mot.matrice);

      //  displaymat(mot.matrice);
      //ofstream outnorm("matnorm.txt");
      ofstream outbayes("matbayes.dat");
      mot.matprec=mot.matrice;


      //Consensus sequence
      for (int i=0;i<width;i++){
         string letter="A";
         double max(mot.matprec[i][0]);
         if (mot.matprec[i][1]>max) letter="T";
         else if (mot.matprec[i][2]>max) letter="C";
         else if (mot.matprec[i][3]>max) letter="G";
         mot.bsinit.append(letter);
      }

      motsdef.push_back(mot);
      if (args_info.coord_file_given||args_info.align_file_given){
         cout << "Loading background set..." << endl;
         ifstream backreg;
         if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000.fa");
         else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs2000-masked.fa");
         regtests=loadsequences(backreg);
         backreg.close();
         cout << "Background set size : " << regtests.size() << endl;
         cout << "Loading training set..." << endl;
         regints=loadseqs();
         cout << "Training set size : " << regints.size() << endl;
         mot.matinit(scorethr2);
         mot.pvaluecomp();
      }
      else {
         mot.pvalue=0;
         mot.scorepoiss=0;
         mot.nbmot=0;
         mot.ntrain=0;
      }


      displaymat(mot.matrice);
      mot.display(outbayes);

      mot.matprec=reversecomp(mot.matrice);
      mot.bsinit="";
      for (int i=0;i<width;i++){
         string letter="A";
         double max(mot.matprec[i][0]);
         if (mot.matprec[i][1]>max) letter="T";
         else if (mot.matprec[i][2]>max) letter="C";
         else if (mot.matprec[i][3]>max) letter="G";
         mot.bsinit.append(letter);
      }

      mot.display(outbayes);

   }		
   else if (args_info.extract_given){
      //!!!!!!!!!!names have to be in alphabetical order!!!

      // training set
      if (args_info.onlyref_given){ // EXTRACTS ONLY REFERENCE SPECIES SEQUENCE
         //         ifstream inf;
         //         inf.open(args_info.coord_file_arg);
         //         extractcoordstofasta(inf);
         //         inf.close();
         cout  << "Extracting sequences" << endl;
         vseq align;
         if (args_info.repeat_masker_given){
            cout << "Using repeat masker..." << endl;
            align=loadseqs();
         }
         else {
            cout << "No repeat masker..." << endl;
            align=loadseqsnotmasked();
         }

         cout << "Writing Sequences..." << endl;
         ofstream outf("seqs.dat");
         for (ivseq ivs=align.begin();ivs!=align.end();ivs++){

            outf << ">";
            outf << ivs->name << "\t";
            outf << chromfromint(ivs->chrom) << "\t";
            outf << ivs->start << "\t";
            outf << ivs->stop << "\n";
            outf << ivs->seqs[0] << "\n";
         }

         outf.close();
      } 
      else if (args_info.cons_motifs_given){ // RETURNS ONLY COORDINATES (OR ALIGNMENT FILES NAMES) OF SEQUENCES MATCHED BY GIVEN MOTIFS

         cout << "Loading Motifs" << endl;
         if (args_info.mots_w_names_given){
            loadmotswnames(args_info.mots_w_names_arg,motsdef);
         } 
         else{
            loadmots(args_info.motifs_arg,motsdef); 
         }
         cout << "Nb mots for score: " << nbmots_for_score  << endl;
         if (nbmots_for_score<motsdef.size()) motsdef.erase(motsdef.begin()+nbmots_for_score,motsdef.end());
         cout << "Loaded " << motsdef.size() << " motifs." << endl;

         if (args_info.opt_thr_file_given){
            cout << "Setting threshold to optimum..." << endl;
            ifstream opt;
            opt.open(args_info.opt_thr_file_arg);
            string dum;
            getline(opt,dum);
            int index=0;
            while (getline(opt,dum)){
               istringstream line(dum);
               line >> dum;
               double thr;
               line >> thr;
               motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
               motsdef[index].motscorethr=motsdef[index].motwidth*(thr-1)/10;
               motsdef[index].motscorethrcons=motsdef[index].motwidth*(thr-1)/10;
               index++;
            }
            opt.close();
         }
         else if (args_info.mean_info_given){
            cout << "Setting threshold to mean info..." << endl;
            for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
               ivm->setscorethr2meaninfo();
               cout << ivm->name << " mean thr = " << ivm->motscorethr2 << endl;
            }
         }
         else {
            cout << "Setting threshold to given value..." << endl;
            cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;
         }

         cout  << "Extracting input sequences..." << endl;
         vseq align;
         if (args_info.repeat_masker_given){
            cout << "Using repeat masker..." << endl;
            align=loadseqs();
         }
         else {
            cout << "No repeat masker..." << endl;
            align=loadseqsnotmasked();
         }

         cout << "Scanning sequences for conserved instances..." << endl;
         ofstream outfc;
         outfc.open("coords-cons-motifs.dat");
         ofstream outfa;
         if (args_info.align_file_given) outfa.open("align-file-cons-motifs.dat");
         unsigned int numseqs(0);

         for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
            int isbound=1;

            for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
               unsigned int nmot;
               if (args_info.wocons_given) nmot=ivm->nbmatchmat(*ivs);
               else nmot=ivm->nbmatchcons(*ivs);

               if (nmot>0) continue;
               else {
                  isbound=0;
                  break;
               }
            }

            if (isbound){
               if (args_info.align_file_given) outfa << ivs->finame << endl ;
               outfc << ivs->name << "\t";
               outfc << chromfromint(ivs->chrom) << "\t";
               outfc << ivs->start << "\t";
               outfc << ivs->stop << "\n";
               numseqs++;
            }
         }
         outfc.close();
         if (args_info.align_file_given) outfa.close();

         cout << "Yields " << numseqs << " sequences out of " << align.size() << "!" << endl;
      }
      else{ // EXTRACTS ALIGNMENTS, WITH (ONLYCONS OR ALSOCONS) OR WITHOUT TAKING INTO ACCOUNT CONSERVATION
         cout << "Loading sequences..." << endl;
         if (args_info.onlycons_given){
            if (args_info.repeat_masker_given){
               cout << "Using repeat masker..." << endl;
               regints=loadseqscons();
            }
            else {
               cout << "No repeat masker..." << endl;
               regints=loadseqsnotmaskedcons();
            }
         } 
         else {
            if (args_info.repeat_masker_given){
               cout << "Using repeat masker..." << endl;
               regints=loadseqs();
            }
            else {
               cout << "No repeat masker..." << endl;
               regints=loadseqsnotmasked();
            }
         }
         cout << "Number of sequences to extract : " << regints.size() << endl;

         cout << "Extracting..." << endl;
         system("if ! test -d align;then mkdir align;fi;");      
         extracttofastawfullname("./align/");
         if (args_info.alsocons_given){

            system("if ! test -d align-cons;then mkdir align-cons;fi;");      
            vseq tmp;
            for (ivseq ivs=regints.begin();ivs!=regints.end();ivs++){
               if (ivs->cons==1) tmp.push_back(*ivs);
            }
            regints=tmp;
            extracttofastawfullname("./align-cons/");
         }
      }

   }
   else if (args_info.pwm2jaspar_given){
      compalpha();
      /*
         typedef vector<double> vd;
         typedef vd::iterator ivd;
         typedef vector<vd> vvd;
         typedef vvd::iterator ivvd;

         istream&
         operator >>(istream &is,vvd &matrice)
         {
         char dummy;
         for (vvd::iterator imat=matrice.begin();imat!=matrice.end();imat++){
         is >> (*imat)[0] >> dummy;
         }
         for (vvd::iterator imat=matrice.begin();imat!=matrice.end();imat++){
         is >> (*imat)[1] >> dummy;
         }
         for (vvd::iterator imat=matrice.begin();imat!=matrice.end();imat++){
         is >> (*imat)[2] >> dummy;
         }
         for (vvd::iterator imat=matrice.begin();imat!=matrice.end();imat++){
         is >> (*imat)[3] >> dummy;
         }
         return is;
         }

         ostream&
         operator <<(ostream &os,const vvd &matrice)
         {
         for (vvd::const_iterator imat=matrice.begin();imat!=matrice.end();imat++){
         os << (*imat)[0] << " ";
         }
         os << "\n";
         for (vvd::const_iterator imat=matrice.begin();imat!=matrice.end();imat++){
         os << (*imat)[1] << " ";
         }
         os << "\n";
         for (vvd::const_iterator imat=matrice.begin();imat!=matrice.end();imat++){
         os << (*imat)[2] << " ";
         }
         os << "\n";
         for (vvd::const_iterator imat=matrice.begin();imat!=matrice.end();imat++){
         os << (*imat)[3] << " ";
         }
         os << "\n";
         return os;
         }
         */
      int width=10;
      ofstream out;
      out.open("jaspar.dat");
      ifstream fmotifs;
      fmotifs.open(argv[1]);
      string dum;
      getline(fmotifs,dum);
      for (int i=1;i<6;i++)
      {
         istringstream fmot(dum);
         for (int j=1;j<8;j++)
         {
            fmot >> dum;
         } 
         vd dumd(4,0.0);
         vvd mat(width,dumd);
         fmot >> mat;
         for (ivvd iv=mat.begin();iv!=mat.end();iv++){
            vd & col=*iv;
            double col0=0.3*exp(col[0]);//A
            double col1=0.3*exp(col[1]);//T
            double col2=0.2*exp(col[2]);//C
            double col3=0.2*exp(col[3]);//G
            //convert in A/C/G/T format.
            col[0]=col0;//floor(100*col0);//A
            col[1]=col2;//floor(100*col2);//C
            col[2]=col3;//floor(100*col3);//G
            col[3]=col1;//floor(100*col1);//T
         }
         out << setprecision(2) << mat << endl;

         getline(fmotifs,dum);
      }
      fmotifs.close();
      out.close();
   }
   else if (args_info.jaspar2pwm_given){
      cout << "Converting file..." << endl;
      jaspar2pwm();
   }
   else if (args_info.transfac2pwm_given){
      cout << "Converting file..." << endl;
      transfac2pwm();
   }
   else if (args_info.furlong2pwm_given){
      cout << "Converting file..." << endl;
      furlong2pwm();
   }
   else if (args_info.jaspar_rank_given){

      cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;

      //      cout << "Loading background set for chi2..." << endl;
      //      ifstream backreg;
      //      if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000-align-1000.dat"); 
      //      else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs2000-align-1000.dat");
      //      regtests=loadsequencesconserv(backreg);
      //      cout << "Loaded " << regtests.size() << " sequences." << endl;
      //      backreg.close();

      cout << "Loading Jaspar Motifs..." << endl;
      loadjaspardb(motsdef);
      cout << "Loaded " << motsdef.size() << " motifs." << endl;

      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
         ivm->setscorethr2meaninfo();
      }

      vseq vscore;
      cout << "Loading pos/neg sequences..." << endl;
      loadseqsforscore(vscore);

      cout << "Defining neg sequences as background set for chi2..." << endl;
      for (ivseq ivs=vscore.begin();ivs!=vscore.end();ivs++){
         if (ivs->sign==-1) regtests.push_back(*ivs);
      }

      cout << "Scanning sequences..." << endl;
      cout << "Motif ";
      ofstream outf("results.dat");
      ofstream outf2("motmeldb.txt");

      outf << "ID" << "\t";
      outf << "Name" << "\t";
      outf << "Consensus" << "\t";
      outf << "enrich" << "\t";
      outf << "nbcons" << "\t";
      outf << "pvalue" << "\t";
      outf << "chi2" << "\n";

      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
         cout << "\r";
         cout << ivm->index+1 << "/" << motsdef.size() << " ";
         cout.flush();
         if (ivm->motwidth>20) continue;

         //reset nb of conserved motifs
         ivm->nbmot=0;
         width=ivm->motwidth;
         scorethr2=ivm->motscorethr2;
         ivm->calclambdaposneg(vscore);
         //         cout << ivm->name << " " << ivm->pvalue << endl;

         double enrich;
         if (ivm->lambdatrain==0 || ivm->lambda==0) enrich=0;
         else enrich=ivm->lambdatrain/ivm->lambda;

         //computing chi2 on conserved sites
         ivm->lambdacompcons();
         ivm->calcscorepoiss();

         outf << ivm->id << "\t";
         outf << ivm->name << "\t";
         outf << ivm->bsinit << "\t";
         outf << enrich << "\t";
         outf << ivm->nbmot << "\t";
         outf << ivm->pvalue << "\t";
         outf << ivm->scorepoiss << "\t";
         //         for (unsigned int i=0;i<distwidth;i++){
         //            outf << ivm->distmot[i];
         //            if (i!=distwidth-1) outf << ",";
         //         }
         //         for (ivinst iv=ivm->instances[0].begin();iv!=ivm->instances[0].end();iv++){
         //            outf << "chr:" << iv->chrom << ":";
         //            outf << iv->coord << "-";
         //            if (iv->sens==1) outf << iv->coord+ivm->motwidth-1;
         //            else outf << iv->coord-ivm->motwidth+1;
         //            if (iv!=ivm->instances[0].end()-1) outf << ",";
         //         }
         outf << endl;

         ivm->id+="_";
         ivm->id+=ivm->name;
         ivm->name=ivm->id;
         ivm->displaywname(outf2);
      }
      cout << endl;
      outf.close();
      outf2.close();
   }
   else if (args_info.backreg_given){

      cout << "Loading genomic alignments..." << endl;

      ifstream align;
      if (species==1) align.open("/home/santolin/these/files/droso/align/all/align.dat");
      else if (species==2) align.open("/home/santolin/these/files/mus/epo/align.dat");//-masked.dat");
      alignscoord=loadcoordconserv(align);
      align.close();

      cout << "Loading intergenic coordinates..." << endl;
      vcoord vcds;

      ifstream intgen;
      if (species==1) intgen.open("/home/santolin/these/files/droso/genome/noncoding/nonCDS.dat");
      else if (species==2) intgen.open("/home/santolin/these/files/mus/all/files/noncoding.dat");
      vcds=loadcoordconserv(intgen);
      intgen.close();

      cout << "Random shuffling intergenic sequences..." << endl;

      random_shuffle(vcds.begin(),vcds.end());

      if (args_info.coord_file_given){
         //cout << "Loading available and conserved sequences from coordinate file..." << endl;
         //regints=loadseqsnotmaskedcons();
         cout << "Loading available sequences from coordinate file..." << endl;
         regints=loadseqsnotmasked();
         cout << "Loaded " << regints.size() << " sequences" << endl;
      }

      cout << "Writing background files..." << endl;
      system("if ! test -d backreg;then mkdir backreg;fi;");      
      string folder="backreg/";//"/home/santolin/these/files/mus/backreg/align/reg_";
      if (args_info.coord_file_given){
         printbackregwcoords(vcds,folder);

      } else{
         printbackreg(vcds,folder);
      }
   }
   else if (args_info.byscore_given){
      compalpha();
      cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;

      cout << "Loading Motifs..." << endl;
      if (args_info.mots_w_names_given){
         loadmotswnames(args_info.motifs_arg,motsdef);
      }
      else {
         loadmots(args_info.motifs_arg,motsdef);
      }
      if (nbmots_for_score<motsdef.size()) motsdef.erase(motsdef.begin()+nbmots_for_score,motsdef.end());
      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
         if (ivm->name.find('_')!=string::npos)
         {
            (*ivm).name.insert((*ivm).name.find('_'),"\\");
         }
      }
      cout << "Loaded " << motsdef.size() << " motifs." << endl;
      cout << "Nb mots for score: " << nbmots_for_score  << endl;

      cout << "Loading sequences..." << endl;
      regints=loadseqs();
      cout << "Number of sequences: " << regints.size() << endl;

      cout << "Computing score..." << endl;
      vseq vscr;
      for (ivseq ivs=regints.begin();ivs!=regints.end();ivs++){
         if (ivs->name.find('_')!=string::npos)
         {
            (*ivs).name.insert((*ivs).name.find('_'),"\\");
         }
         // for motis filling
         calcscore(motsdef,*ivs);
         if (ivs->motis[0]>0){// || ivs->motis[1]>0){}
            //        // if (ivs->motis[2]>0 || ivs->motis[3]>0){}
            ////         if (ivs->motis[7]>0){}
            //            vmot motstmp;
            //            motstmp.push_back(motsdef[0]);
            //    //        motstmp.push_back(motsdef[1]);
            ////            motstmp.push_back(motsdef[7]);
            //            vint motis;
            //            motis=ivs->motis;
            //            calcscore(motstmp,*ivs);
            //            ivs->motis=motis;
            vscr.push_back(*ivs);
         }
         //      }
   }

   //   cout << "Adding scores for each gene..." << endl;
   //   vseq vuniq;
   //   vuniq=addscorepergene(vscr);
   //   vscr=vuniq;

   cout << "Sorting" << endl;
   sort(vscr.begin(),vscr.end());


   cout << "Writing file ..." << endl;

   ofstream disp;
   ostringstream disps;
   disps << "ranking.tex";
   disp.open(disps.str().c_str());
   dispscoreforrank(vscr,motsdef,disp);
   disp.close();
   ofstream outf("ranking.dat");
   for (ivseq ivs=vscr.begin();ivs!=vscr.end();ivs++){
      if (ivs->name.find('_')!=string::npos)
      {
         (*ivs).name.erase((*ivs).name.find('\\'),2);
      }
      outf << ivs->name << "\t";
      outf << chromfromint(ivs->chrom) << "\t";
      outf << ivs->start << "\t";
      outf << ivs->stop << "\n";
      //      outf << ivs->motis[0] << "\t";
      //      outf << ivs->motis[1] << "\t";
      //      outf << ivs->score << "\n";
   }
   outf.close();

}
   else if (args_info.rank_given){
      
      cout << "Loading Motifs..." << endl;
      if (args_info.mots_w_names_given){
         loadmotswnames(args_info.mots_w_names_arg,motsdef);
      }
      else {
         loadmots(args_info.motifs_arg,motsdef);
      }
      cout << "Loaded " << motsdef.size() << " motifs." << endl;
      nbmots_for_score=motsdef.size();

      cout << "Loading sequences..." << endl;
      loadseqsforscorenotmasked(regints);
      cout << "Number of sequences: " << regints.size() << endl;

      if (args_info.fdr_given){
         cout << "Computing score..." << endl;
         ofstream outf("motifs-tpr-fdr.txt");
         double fdr=args_info.fdr_arg;
         for (ivmot iv=motsdef.begin();iv!=motsdef.end();iv++){
            calcoptthr(*iv,fdr,outf);
         }
         outf.close();
      }
}
else if (args_info.neargene_given){

   //!! Here we load TSS, i.e start=stop
   cout << "Loading genes coordinates..." << endl;
   vcoord vcdg;
   ifstream gen;
   if (species==1){
      gen.open("/home/santolin/these/files/droso/annot/coords/regres-coords.dat");
      vcdg=loadcoordpeaks(gen);
   }
   else if (species==2){
      //gen.open("/home/santolin/these/files/mus/biomart/genes-n-strand-protein-coding+miRNA.dat");
      //               gen.open("/home/santolin/these/files/mus/biomart/genes-n-strand-only-protein-coding.dat");
      //               vcdg=loadcoordconservwstrand(gen);
      gen.open("/home/santolin/these/files/mus/biomart/genes-n-strand-only-protein-coding-TSS.dat");
      vcdg=loadcoordfromTSS(gen);
   }
   gen.close();

   if (args_info.define_peak_given){
      cout << "Loading peaks coordinates..." << endl;
      vcoord vcdp;
      ifstream pea;
      pea.open(args_info.coord_file_arg);
      vcdp=loadcoordpeaks(pea);
      pea.close();

      cout << "Finding nearest gene..." << endl;
      findnearestgene(vcdg,vcdp);

      cout << "Writing file ..." << endl;
      ofstream outf("neargene.dat");
      for (ivcoord ivc=vcdp.begin();ivc!=vcdp.end();ivc++){
         outf << ivc->name << "\t";
         outf << chromfromint(ivc->chrom) << "\t";
         outf << ivc->start-100 << "\t";
         outf << ivc->start+100 << "\n";
      }
      outf.close();
   } 
   else {
      cout << "Loading coordinates..." << endl;
      vcoord vcdp;
      ifstream pea;
      pea.open(args_info.coord_file_arg);
      vcdp=loadcoordconserv(pea);
      cout << "Loaded " << vcdp.size() << " coords." << endl;
      pea.close();

      cout << "Finding nearest gene..." << endl;
      if (args_info.on_each_side_given){
         findnearestgene_sides(vcdg,vcdp);
      }
      else{
         findnearestgene(vcdg,vcdp);
      }

      cout << "Writing file ..." << endl;
      ofstream outf("neargene.dat");
      if (args_info.on_each_side_given){
         for (ivcoord ivc=vcdp.begin();ivc!=vcdp.end();ivc++){
            outf << ivc->name << "\t";
            outf << chromfromint(ivc->chrom) << "\t";
            outf << ivc->start << "\t";
            outf << ivc->stop << "\t";
            for (ivstring ivs=ivc->neargenes.begin();ivs!=ivc->neargenes.end();ivs++){
               outf << *ivs << "\t";
            }
            outf << "\n";
         }
      }
      else{
         for (ivcoord ivc=vcdp.begin();ivc!=vcdp.end();ivc++){
            outf << ivc->name << "\t";
            outf << chromfromint(ivc->chrom) << "\t";
            outf << ivc->start << "\t";
            outf << ivc->stop << "\n";
         }
      }
      outf.close();


   }
}
else if (args_info.comp_given){


   int stat1(0),stat2(0);
   vstring vpos1,vpos2;
   vcoord vcds1,vcds2;

   if (args_info.refnames_given){
      cout << "Loading ref genes names..." << endl;
      string dum;
      ifstream pos(args_info.refnames_arg);
      getline(pos,dum);
      while (!pos.eof()){
         vpos1.push_back(dum);
         getline(pos,dum);
      }
      pos.close();
      stat1=1;
   } 
   else if (args_info.refcoords_given){
      cout << "Loading ref genes coordinates..." << endl;
      ifstream gen;
      gen.open(args_info.refcoords_arg);
      vcds1=loadcoordconserv(gen);
      sort(vcds1.begin(),vcds1.end());
      //               for (ivcoord ivc=vcds1.begin();ivc!=vcds1.end();ivc++){
      //                  cout << *ivc << endl;
      //               }
      gen.close();
      stat1=2;
   } 
   else {
      cout << "Please give a ref file. Exiting." << endl;
      exit(2);
   }


   if (args_info.compnames_given){
      cout << "Loading comp genes names..." << endl;
      string dum;
      ifstream pos(args_info.compnames_arg);
      getline(pos,dum);
      while (!pos.eof()){
         vpos2.push_back(dum);
         getline(pos,dum);
      }
      pos.close();
      stat2=1;
   } 
   else if (args_info.compcoords_given){
      cout << "Loading comp genes coordinates..." << endl;
      ifstream gen;
      gen.open(args_info.compcoords_arg);
      vcds2=loadcoordconserv(gen);
      sort(vcds2.begin(),vcds2.end());
      gen.close();
      stat2=2;
   } 
   else {
      cout << "Please give a comp file. Exiting." << endl;
      exit(2);
   }

   cout << "Comparing and writing file..." << endl;
   ofstream outf("genecomp.dat");
   if (stat1==1 && stat2==2){
      for (ivcoord ivc=vcds2.begin();ivc!=vcds2.end();ivc++){
         for (ivstring ivs=vpos1.begin();ivs!=vpos1.end();ivs++){
            if (*ivs==ivc->name){
               outf << *ivc;
               break;
            }
         }
      }
   } 
   else if (stat1==2 && stat2==2){

      int cutoff;
      cutoff=args_info.cutoff_for_comp_arg;

      //      outf << "Ref_name\tChrom\tStart\tStop\tComp_name\n";
      int crmcount=1;
      int start,stop;
      int pstart,pstop,pchrom;
      string pname;
      int pending=0;
      int numchip(1);
      pstart=vcds2.begin()->start;
      pstop=vcds2.begin()->stop;
      pchrom=vcds2.begin()->chrom;
      pname=vcds2.begin()->name;
      for (ivcoord ivc=vcds2.begin();ivc!=vcds2.end();ivc++){
         for (ivcoord ivs=vcds1.begin();ivs!=vcds1.end();ivs++){
            if (ivs->chrom==ivc->chrom){

               //if ((ivs->start<=ivc->start && ivs->stop>=ivc->stop) || (ivs->start <= ivc->stop+1000 && ivs->start >= ivc->start) || (ivs->stop >= ivc->start-1000 && ivs->stop <= ivc->stop)){}
               // 2 inclus dans 1, 2 à gauche de 1, 2 à droite de 1, 
               if (  (ivs->start<=ivc->start && ivs->stop>=ivc->stop) || 
                     (ivs->start <= ivc->stop && ivs->start >= ivc->start) ||
                     (ivs->stop >= ivc->start && ivs->stop <= ivc->stop) || 
                     ((ivc->start+ivc->stop)/2-ivs->stop<=cutoff && (ivc->start+ivc->stop)/2-ivs->stop>0 ) || // CAS DE CHIP VS GENE : centre CHIP a moins de cutoff de fin gene
                     (ivs->start-(ivc->start+ivc->stop)/2<=cutoff && ivs->start-(ivc->start+ivc->stop)/2>0) ){ // CAS DE CHIP VS GENE : centre CHIP a moins de cutoff de debut gene
                     //(abs((ivs->start+ivs->stop)/2-(ivc->start+ivc->stop)/2)<=cutoff) ){} // CAS DE CHIP VS CHIP : centres a moins de cutoff=1kb

                  start =  min(ivs->start,ivc->start);
                  stop  =  max(ivs->stop,ivc->stop);
                  numchip++;

                  //check we do not cover previous crm
                  if ((start <= pstop && start >= pstart) || (stop >= pstart && stop <=pstop)){
                     if (pchrom==ivc->chrom){
                        pstart=min(start,pstart);
                        pstop=max(stop,pstop);
                     }
                     else {
                        pstart=start;
                        pstop=stop;
                     }
                     pchrom=ivc->chrom;
                     pending=1;
                  }
                  else {

                     outf << chromfromint(pchrom) << "\t";
                     outf << pstart << "\t";
                     outf << pstop << "\t";
                     //outf << "crm_" << crmcount << "\n";
                     outf << pname << "\n";
                     //                        outf << numchip << "\n";
                     pstart=start;
                     pstop=stop;
                     pchrom=ivc->chrom;
                     pname=ivc->name;
                     crmcount++;
                     pending=0;
                     numchip=1;
                     //ivs->check=1;
                     //             break;

                  }

               }
            }
         }
      }
      if (pending){
         outf << chromfromint(pchrom) << "\t";
         outf << pstart << "\t";
         outf << pstop << "\t";
         //outf << "crm_" << crmcount << "\n";
         outf << pname << "\n";
      }

      //
      //                  ofstream outf1("genecomp-ref-wo-comp.dat");
      //                  for (ivcoord ivs=vcds1.begin();ivs!=vcds1.end();ivs++){
      //                     if (ivs->check==0){
      //                        outf1 << ivs->name << "\t";
      //                        outf1 << chromfromint(ivs->chrom) << "\t";
      //                        outf1 << ivs->start << "\t";
      //                        outf1 << ivs->stop << "\n";
      //                     }
      //                  }
      //                  outf1.close();
   } 
   else if (stat1==1 && stat2==1){
      for (ivstring ivc=vpos1.begin();ivc!=vpos1.end();ivc++){
         for (ivstring ivs=vpos2.begin();ivs!=vpos2.end();ivs++){
            if (*ivs==*ivc){
               outf << *ivs << endl;
               break;
            }
         }
      }
   }
   else if (stat1==2 && stat2==1){
      for (ivcoord ivc=vcds1.begin();ivc!=vcds1.end();ivc++){
         for (ivstring ivs=vpos2.begin();ivs!=vpos2.end();ivs++){
            if (*ivs==ivc->name){
               outf << *ivc;
               break;
            }
         }
      }
   }
   outf.close();


   //   cout << "Writing file ..." << endl;
   //   ofstream outf("neargene.dat");
   //   for (ivcoord ivc=vcdp.begin();ivc!=vcdp.end();ivc++){
   //      outf << ivc->name << "\t";
   //      outf << chromfromint(ivc->chrom) << "\t";
   //      outf << ivc->start-100 << "\t";
   //      outf << ivc->start+100 << "\n";
   //   }
   //   outf.close();
}
else if (args_info.stats_given){
   cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;

   cout << "Loading Motifs..." << endl;
   loadmots(args_info.motifs_arg,motsdef);
   cout << "Loaded " << motsdef.size() << " motifs." << endl;
   cout << "Nb mots for stats: " << nbmots_for_score  << endl;

   cout << "Loading alignments " << endl;
   vseq align;
   align=loadseqsnotmasked();
   cout << "Nb sequences to scan: " << align.size() << endl;

   cout << "Loading background set..." << endl;
   ifstream backreg;
   if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000-align-short.dat");
   else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs2000-align-short.dat");
   regtests=loadsequencesconserv(backreg);
   backreg.close();
   cout << "Background set size : " << regtests.size() << endl;

   cout << "Scanning sequences for instances..." << endl;
   scanseqsforinstancesnmask(align,motsdef);

   cout << "Defining conserved instances..." << endl;
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
      ivs->instances2instancescons();
      cout << ivs->instancescons;
   }

   vseq vbest;
   //      cout << "Keeping best 1kbs..." << endl;
   //      vbest=findbestseqs(motsdef,align);
   //      
   //      for (ivseq iv=vbest.begin();iv!=vbest.end();iv++){
   //         iv->instances.clear();
   //         iv->instancescons.clear();
   //      }
   //      
   //      cout << "Scanning best 1kbs for instances..." << endl;
   //      scanseqsforinstancesnmask(vbest,motsdef);
   //
   //      cout << "Defining best 1kbs conserved instances..." << endl;
   //      for (ivseq ivs=vbest.begin();ivs!=vbest.end();ivs++){
   //         ivs->instances2instancescons();
   //      }
   vbest=align;

   if (args_info.counts_given){
      cout << "Counts statistics..." << endl;
      statcounts(align,motsdef);
   } 
   else if(args_info.shifts_given){
      cout << "Shifts statistics..." << endl;
      system("if ! test -d data;then mkdir data;fi;");      
      statshifts(align,motsdef);
   } 
   else {
      cout << "Please give an option" << endl;
   }
}
else if (args_info.surround_given){

   cout << "Loading Motifs" << endl;
   if (args_info.mots_w_names_given){
      loadmotswnames(args_info.mots_w_names_arg,motsdef);
   } 
   else{
      loadmots(args_info.motifs_arg,motsdef); 
   }
   if (args_info.mean_info_given){
      cout << "Setting scorethr2 to mean info..." << endl;
      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
         ivm->setscorethr2meaninfo();
      }
   }
   cout << "Nb mots for score: " << nbmots_for_score  << endl;
   if (nbmots_for_score<motsdef.size()) motsdef.erase(motsdef.begin()+nbmots_for_score,motsdef.end());
   cout << "Loaded " << motsdef.size() << " motifs." << endl;

   if (args_info.opt_thr_file_given){
      ifstream opt;
      opt.open(args_info.opt_thr_file_arg);
      string dum;
      getline(opt,dum);
      int index=0;
      while (getline(opt,dum)){
         istringstream line(dum);
         line >> dum;
         double thr;
         line >> thr;
         motsdef[index].motscorethr2=motsdef[index].motwidth*thr/10;
         motsdef[index].motscorethr=motsdef[index].motwidth*(thr-1)/10;
         motsdef[index].motscorethrcons=motsdef[index].motwidth*(thr-1)/10;
         index++;
      }
      opt.close();
   }
   else {
      scorethr2=motsdef[0].motscorethr2;
      scorethr=motsdef[0].motscorethr2;
      scorethrcons=motsdef[0].motscorethr2;
      cout << "Thresholds: thr2=" << scorethr2 << " thr=" << scorethr << " thrcons=" << scorethrcons << endl;
   }
   if (motsdef.size()==1){
      cout << "One motif auto-correlation" << endl;
   }
   else if (motsdef.size()==2){
      cout << "Two motifs cross-correlation" << endl;
   }
   else if (motsdef.size()>1){
      cout << "Please give a maximum of two motifs. Exiting." << endl;
      exit(1);
   }

   cout << "Loading alignments " << endl;
   vseq align;
   align=loadseqsnotmasked();
   cout << "Nb sequences to scan: " << align.size() << endl;

   //      cout << "Loading background set..." << endl;
   //      ifstream backreg;
   //      if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000-align-short.dat");
   //      else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs2000-align-short.dat");
   //      regtests=loadsequencesconserv(backreg);
   //      backreg.close();
   //      cout << "Background set size : " << regtests.size() << endl;

   cout << "Scanning sequences for instances..." << endl;
   vmot motstmp;
   motstmp.push_back(motsdef[0]);
   scanseqsforinstances(align,motstmp);

   cout << "Defining conserved instances..." << endl;
   for (ivseq ivs=align.begin();ivs!=align.end();ivs++){
      ivs->instances2instancescons();
      //cout << ivs->instancescons;
   }

   cout << "Surrounding score statistics..." << endl;
   surrounding_score_stat(align,motsdef);
}
else if (args_info.testevol_given){
   // tests similarity between ref_species_only PWM and evolved_PWM
   //
   ofstream motmeldb;
   motmeldb.open("motmeldb.txt");
   cout << "Loading Motif(s)..." << endl;
   if (args_info.mots_w_names_given){
      loadmotswnames(args_info.mots_w_names_arg,motsdef);
   }
   else {
      loadmots(args_info.motifs_arg,motsdef);
   }
   if (args_info.mean_info_given){
      cout << "Setting scorethr2 to mean info..." << endl;
      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
         ivm->setscorethr2meaninfo();
      }
   }

   for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      ivm->motscorethrcons=ivm->motscorethr2;
      ivm->motscorethr=ivm->motscorethr2;
      cout << ivm->name << " thr = " << ivm->motscorethr2 << endl;
   }
   scorethr2=motsdef[0].motscorethr2;
   cout << "scorethr2=" << scorethr2 << endl;
   scorethr=scorethr2;
   scorethrcons=scorethr2;
   width=motsdef[0].motwidth;
   compalpha();

   //cout << scorethr2 << " " << alpha << " " << width << endl;

   if (args_info.synthetic_given){


      cout << "Loading alignments " << endl;
      regints=loadseqs();

      cout << "Nb sequences to scan: " << regints.size() << endl;

      //                  for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      //                     ivm->matinit(ivm->motscorethr2);
      //                  }

      motsdef[0].comprefinstances(regints,0);
      motsdef[0].comprefmot();
      motsdef[0].setscorethr2meaninfo();
      cout << motsdef[0].name << " refined on ref, new thr = " << motsdef[0].motscorethr2 << endl;
      motsdef[0].matinit(motsdef[0].motscorethr2);
      cout << motsdef[0].name << " has " << motsdef[0].refinstances.seqs.size() << " instances." << endl;
      cout << motsdef[0].name << " has " << motsdef[0].seqs.size() << " conserved instances." << endl;
      //                  for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      //                     ivm->matinit(ivm->motscorethr2);
      //                     cout << ivm->name << " has " << ivm->seqs.size() << " conserved instances." << endl;
      //                  }

      Motif mot=motsdef[0];
      //                  for (ivma ivm=mot.seqs.begin();ivm!=mot.seqs.end();ivm++){
      //                     ivm->print();
      //                  }
      //                  exit(9);
      if (args_info.fitdist_given){
         //fitdistperbase(mot);
         fitdistpersite(mot);

      }
      else if (args_info.fitscore_given){
         //fitdistperbase(mot);
         fitscorepersite(mot);

      }
      else if (args_info.fitkmeans_given){
         fitkmeanspersite(mot,args_info.fitkmeans_arg);
      }
      else if (args_info.fitdist_kmeans_given){
         fitdistkmeanspersite(mot,args_info.fitdist_kmeans_arg);
      }
      else{
         if (motsdef.size()==1){
            //evolvebase(mot);
            evolvesite(mot);
            //phylotest(mot);
         } 
         else{
            // FOR MIXTURE MODEL
            evolvesite(motsdef);
         }
      }
   }
   else {

      cout << "Loading background set..." << endl;
      ifstream backreg;
      if (species==1) backreg.open("/home/santolin/these/files/droso/backreg/regs2000.fa");
      else if (species==2) backreg.open("/home/santolin/these/files/mus/backreg/regs2000.fa");
      regtests=loadsequences(backreg);
      backreg.close();
      cout << "Background set size : " << regtests.size() << endl;

      cout << "Loading alignments " << endl;
      regints=loadseqs();
      cout << "Nb sequences to scan: " << regints.size() << endl;


      inittreedist();
      scanseqsforinstances(regints,motsdef);
      ofstream outf;
      if (args_info.pseudocount_given){
         outf.open("testcv.dat");
         width=motsdef[0].motwidth;
         for (double al=0.01;al<1;al+=.1){
            alpha=al;
            beta=concc/conca*alpha;
            motsdef[0].testpseudocount(outf);
         }
         outf.close();
      }
      else if (args_info.timer_given){
         outf.open("testcv-timer.dat");
         motsdef[0].testpwmcvtimer(outf);
         outf.close();
      }
      else{
         outf.open("testcv.dat");
         cout << "Testing convergence..." << endl;
         motsdef[0].testpwmcv(outf,motmeldb);
         outf.close();
      }
   }
}
else if (args_info.countwords_given){
   // First motif = total PWM
   // Following motifs = with k-means


   cout << "Loading Motifs..." << endl;
   if (args_info.mots_w_names_given){
      loadmotswnames(args_info.mots_w_names_arg,motsdef);
   }
   else {
      loadmots(args_info.motifs_arg,motsdef);
   }
   if (args_info.mean_info_given){
      cout << "Setting scorethr2 to mean info..." << endl;
      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
         ivm->setscorethr2meaninfo();
         cout << ivm->name << " mean info = " << ivm->motscorethr2 << endl;
      }
   }
   Motif mot=motsdef[0];
   if (args_info.kprobs_given || args_info.numk_given)  cout << "Using k-means generated motifs to fit sites distribution" << endl;
   //displaymat(mot.matfreq);
   scorethr2=mot.motscorethr2;
   width=mot.motwidth;
   compalpha();

   cout << "Loading alignments " << endl;
   regints=loadseqs();
   cout << "Nb sequences to scan: " << regints.size() << endl;

   if (args_info.kprobs_given){
      vd kprobs;
      ifstream inf;
      inf.open(args_info.kprobs_arg);
      back_insert_iterator<vd> dest(kprobs);
      copy(iidouble(inf),iidouble(),dest);
      inf.close();
      countwords_mixture(motsdef,kprobs);
   } 
   else if (args_info.numk_given){
      unsigned int numk=args_info.numk_arg;
      countwords_numk(mot,numk);
   }
   else if (args_info.compare_pseudo_count_given){
      countwords_compare_pseudocount(mot);
   }
   else countwords(mot);

}
else if (args_info.nmers_given){

   cout << "Loading alignments " << endl;
   regints=loadseqs();
   cout << "Nb sequences to scan: " << regints.size() << endl;

   // just display nmers, sort and counting procedures can be done easily in bash
   ofstream outf("nmers.dat");
   countnmers(outf);
   outf.close();

}
else if (args_info.numseq_given){

   cout << "Loading Motif(s)..." << endl;
   if (args_info.mots_w_names_given){
      loadmotswnames(args_info.mots_w_names_arg,motsdef);
   }
   else {
      loadmots(args_info.motifs_arg,motsdef);
   }
   if (args_info.mean_info_given){
      cout << "Setting scorethr2 to mean info..." << endl;
      for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
         ivm->setscorethr2meaninfo();
      }
   }

   for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      ivm->motscorethrcons=ivm->motscorethr2;
      ivm->motscorethr=ivm->motscorethr2;
      cout << ivm->name << " thr = " << ivm->motscorethr2 << endl;
   }

   nbmots_for_score=motsdef.size();

   cout << "Loading alignments " << endl;
   regints=loadseqs();
   cout << "Nb sequences to scan: " << regints.size() << endl;

   ofstream outf("numcons.dat");
   outf << "numseq" << "\t";
   outf << regints.size() << "\n";
   for (ivmot ivm=motsdef.begin();ivm!=motsdef.end();ivm++){
      ivm->matinit(ivm->motscorethr2);
      unsigned int numcons=0;
      unsigned int num=0;
      for (ivseq ivs=regints.begin();ivs!=regints.end();ivs++){
         if (ivs->motis[ivm->index]>0) numcons++;
         if (ivs->nmot>0) num++;
      }
      outf << ivm->name << "\t";
      outf << num << "\t";
      outf << numcons << "\n";
   }
   outf.close();

}
else if (args_info.test_given){
   cout << "Testing rng convergence..." << endl;
   rngtest();
}
else {
   cout << "No mode given" << endl;
}


gsl_rng_free(gslran);
cmdline_parser_free(&args_info);
cout << "exit normally" << endl;
return 0;
}

