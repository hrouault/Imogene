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
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <dirent.h>

#include "vectortypes.hpp"
#include "sequence.hpp"
#include "tree.hpp"

vint lengthchrom;

vcoord alignscoord;

bool 
operator<(const Chromosome & chr1,const Chromosome & chr2)
{
   return intfromchrom(chr1.name) < intfromchrom(chr2.name); 
}

bool 
operator<(const Coordinate & coord1,const Coordinate & coord2)
{
   if (coord1.chrom==coord2.chrom) return (coord1.start<=coord2.start);
   else return coord1.chrom < coord2.chrom;

}
   
   int
compN(vint & bs)
{
   int N=0;
   for (vint::const_iterator ibs=bs.begin();ibs!=bs.end();ibs++){
      if (*ibs==4) N++;
   }
   return N;
}

   string
inttostring(int iseq)
{
   string seq;
   if (iseq==0){
      seq="A";
   }
   else if (iseq==1){
      seq="C";
   }
   else if (iseq==2){
      seq="G";
   }
   else if (iseq==3){
      seq="T";
   }
   else if (iseq==4){
      seq="N";
   }
   else if (iseq==5){
      seq="-";
   }
   else {
      seq="?";
   }

   return seq;
}
   
   string
vinttostring(vint & iseq)
{
   vint::const_iterator istr;

   string seq;
   for (istr=iseq.begin();istr!=iseq.end();istr++){
      seq+=inttostring(*istr);
   }
   return seq;
}

   vint
stringtoint(string & seq)
{
   string::const_iterator istr;

   vint iseq;
   for (istr=seq.begin();istr!=seq.end();istr++){
      if (*istr=='A'||*istr=='a'){
         iseq.push_back(0);
      }
      else if (*istr=='C'||*istr=='c'){
         iseq.push_back(1);
      }
      else if (*istr=='G'||*istr=='g'){
         iseq.push_back(2);
      }
      else if (*istr=='T'||*istr=='t'){
         iseq.push_back(3);
      }
      else if (*istr=='-'){
         iseq.push_back(5);
      }
      else {
         iseq.push_back(4);
      }
   }
   return iseq;
}

   string
remgaps(string & seq)
{
   string::const_iterator istr;

   string oseq;
   for (istr=seq.begin();istr!=seq.end();istr++){
      if (*istr!='-'){
         oseq.push_back(*istr);
      }
   }
   return oseq;
}

TSS::TSS(int position,char dir,string chr,string gname)
{
   if (chr=="2L"){
      chrom=0;
   } else if (chr=="2R"){
      chrom=1;
   } else if (chr=="3L"){
      chrom=2;
   } else if (chr=="3R"){
      chrom=3;
   } else if (chr=="4"){
      chrom=4;
   } else if (chr=="X"){
      chrom=5;
   }
   coord=position;
   if (dir=='+'){
      sens=1;
   } else sens=-1;
   gene=gname;
}

string
chromfromint(int chr)
{
   if (species=="droso"){
      if (chr==0){
         return "2L";
      } else if (chr==1){
         return "2R";
      } else if (chr==2){
         return "3L";
      } else if (chr==3){
         return "3R";
      } else if (chr==4){
         return "4";
      } else if (chr==5){
         return "X";
      } else {
         return "unknown";
      }
   }
   else if (species=="eutherian"){
      if (chr==0){
         return "1";
      } else if (chr==1){
         return "2";
      } else if (chr==2){
         return "3";
      } else if (chr==3){
         return "4";
      } else if (chr==4){
         return "5";
      } else if (chr==5){
         return "6";
      } else if (chr==6){
         return "7";
      } else if (chr==7){
         return "8";
      } else if (chr==8){
         return "9";
      } else if (chr==9){
         return "10";
      } else if (chr==10){
         return "11";
      } else if (chr==11){
         return "12";
      } else if (chr==12){
         return "13";
      } else if (chr==13){
         return "14";
      } else if (chr==14){
         return "15";
      } else if (chr==15){
         return "16";
      } else if (chr==16){
         return "17";
      } else if (chr==17){
         return "18";
      } else if (chr==18){
         return "19";
      } else if (chr==19){
         return "X";
      } else if (chr==20){
         return "Y";
      } else {
         return "unknown";
      }
   }
}

int
intfromchrom(string chrname)
{
   if (species=="droso"){	
      if (chrname=="2L"){
         return 0;
      } else if (chrname=="2R"){
         return 1;
      } else if (chrname=="3L"){
         return 2;
      } else if (chrname=="3R"){
         return 3;
      } else if (chrname=="4"){
         return 4;
      } else if (chrname=="X"){
         return 5;
      } else {
         return -1;
      }
   }
   else if (species=="eutherian"){	
      if (chrname=="1"){
         return 0;
      } else if (chrname=="2"){
         return 1;
      } else if (chrname=="3"){
         return 2;
      } else if (chrname=="4"){
         return 3;
      } else if (chrname=="5"){
         return 4;
      } else if (chrname=="6"){
         return 5;   
      } else if (chrname=="7"){
         return 6;   
      } else if (chrname=="8"){
         return 7;   
      } else if (chrname=="9"){
         return 8;   
      } else if (chrname=="10"){
         return 9;   
      } else if (chrname=="11"){
         return 10;   
      } else if (chrname=="12"){
         return 11;   
      } else if (chrname=="13"){
         return 12;   
      } else if (chrname=="14"){
         return 13;   
      } else if (chrname=="15"){
         return 14;   
      } else if (chrname=="16"){
         return 15;   
      } else if (chrname=="17"){
         return 16;   
      } else if (chrname=="18"){
         return 17;   
      } else if (chrname=="19"){
         return 18;   
      } else if (chrname=="X"){
         return 19;   
      } else if (chrname=="Y"){
         return 20;   
      } else {
         return -1;
      }
   }
}

TSS::TSS()
{
}

istream &
operator >>(istream &is,TSS & tss)
{
   int pos;
   char dum;
   char dir;
   string gname;
   string chrom;
   is >> pos;
   is >> dir;
   is >> gname;
   is >> chrom;

   if (dir=='+'){
      tss.sens=1;
   } else if (dir=='-'){
      tss.sens=-1;
   } else if (!is.eof()){
      cout << "error!! pos : " << pos << "\n";
   }
   tss.chrom=intfromchrom(chrom);
   
   tss.coord=pos;
   tss.gene=gname;
   return is;
}

istream &
operator >>(istream &is,Chromosome & chrom)
{
   getline(is,chrom.name);
   if (chrom.name.find('>')!=string::npos){
      chrom.name=chrom.name.substr(1);
   }
   chrom.name=chrom.name.substr(0,chrom.name.find(" "));
   string dumseq="";
   getline(is,dumseq,'>');

   string::iterator ist=dumseq.begin();

   chrom.seq="";
   while (ist!=dumseq.end())
   {
      if (*ist!='\n'){
         chrom.seq.push_back(*ist);
      }
      ist++;
   }

   return is;
}

istream &
operator >>(istream &is,Sequence & seq)
{ 
   string dummystr;
   vint dummyvint;

   seq.seqsrealigned.insert(seq.seqsrealigned.begin(),nbspecies,dummystr);
   seq.imaps.insert(seq.imaps.begin(),nbspecies,dummyvint);
   seq.imapsinv.insert(seq.imapsinv.begin(),nbspecies,dummyvint);
   seq.species.insert(seq.species.begin(),nbspecies,0);
   seq.seqs.insert(seq.seqs.begin(),nbspecies,dummystr);
   seq.iseqs.insert(seq.iseqs.begin(),nbspecies,dummyvint);

   string fseqline;
   getline(is,fseqline);

   stringstream firstline(fseqline);
   string dros;
   firstline >> dros;
   string chrom;
   firstline >> chrom;
   seq.chrom=intfromchrom(chrom);
//	cout << chrom << " ";
   firstline >> seq.start;
  // cout << seq.start << " ";
   firstline >> seq.stop;
  // cout << seq.stop << "\n";
   while(!is.eof()){
      //cout << fseqline << endl;
	   int numdro=speciestonum(fseqline.substr(1,6)); 
	   //      cout << numdro << "\n";
      getline(is,fseqline);

      seq.seqsrealigned[numdro]=fseqline;
      seq.imaps[numdro]=alignedtomap(fseqline);
      seq.imapsinv[numdro]=alignedtorevmap(fseqline);
      string seqwogap=remgaps(fseqline);
      seq.seqs[numdro]=seqwogap;
      seq.iseqs[numdro]=stringtoint(seqwogap);
      int nbtb=seq.iseqs[numdro].size()-compN(seq.iseqs[numdro]);
      if (nbtb>width+neighbext){//>5){
         seq.species[numdro]=1;
      } else {
         seq.species[numdro]=0;
      }
      getline(is,fseqline);
   }
   seq.nbN=compN(seq.iseqs[0]);
   seq.nbtb=seq.iseqs[0].size()-seq.nbN;
   return is;
}

istream &
operator >>(istream &is,Coordinate &coord)
{
   is >> coord.name;
   string chromname;
   is >> chromname;
   coord.chrom=intfromchrom(chromname);
   is >> coord.start;
   is >> coord.stop;

   return is;
}

ostream&
operator <<(ostream &os,const Sequence & s)
{
      os << s.name << "\t";
      os << "score: " << s.score << "\t";
      os << "chrom: " << chromfromint(s.chrom) << "\t";
      os << "start: " << s.start << "\t";
      os << "stop: " << s.stop << "\t";
      os << "length: " << s.iseqs[0].size() << "\t";
      if (s.motis.size()!=0){
         unsigned int motindex(1);
         for (civint iv=s.motis.begin();iv!=s.motis.end();iv++){   
            os << "#"<< motindex << ": " << *iv << "\t";
            motindex++;
         }
      }
      if (s.sign==1){
         os << "+" ;
      }
      else if (s.sign==-1){
         os << "-" ;
      }
      os << endl;
      return os;
}

   ostream&
operator <<(ostream &os,const vseq &vs)
{
   for (civseq ivs=vs.begin();ivs!=vs.end();ivs++){
      os << (*ivs).name << "\t";
      os << "score: " << (*ivs).score << "\t";
      os << "chrom: " << chromfromint(ivs->chrom) << "\t";
      os << "start: " << ivs->start << "\t";
      os << "stop: " << ivs->stop << "\t";
      os << "length: " << ivs->iseqs[0].size() << "\t";
      if ((*ivs).motis.size()!=0){
         unsigned int motindex(1);
         for (civint iv=(*ivs).motis.begin();iv!=(*ivs).motis.end();iv++){   
            os << "#"<< motindex << ": " << *iv << "\t";
            motindex++;
         }
      }
      if ((*ivs).sign==1){
         os << "+" ;
      }
      else if ((*ivs).sign==-1){
         os << "-" ;
      }
      os << endl;
   }
   return os;
}

ostream &
operator <<(ostream &os,Coordinate &coord)
{
   os << coord.name << "\t";
   os << chromfromint(coord.chrom) << "\t";
   os << coord.start << "\t";
   os << coord.stop << "\n";
//   os << chromfromint(coord.chrom) << "\t";
//   os << coord.start << "\t";
//   os << coord.stop << "\t";
//   os << coord.name << "\n";

   return os;
}

ostream &
operator <<(ostream &os,Instanceseq &ist)
{
  // os << "Motif #" << ist.motindex+1 << "\t";
   os <<  ist.motname << "\t";
   os <<  ist.seq << "\t";
   os << "sens " << ist.sens << "\t";
   os << "pos (w/gaps) " << ist.pos << "\t";
   os << "pos " << ist.truepos << "\t";
   os << "score " << ist.score << "\t";
   os << "species " << ist.species ;
   os << " (" << numtospecies(ist.species) << ")" << "\t";
   os << "\n";

   return os;
}

   ostream &
operator <<(ostream &os,vinstseq &vist)
{

   for (ivinstseq ivs=vist.begin();ivs!=vist.end();ivs++){
      os << *ivs;
   }

   return os;
}

   ostream &
operator <<(ostream &os,vvinstseq &vvist)
{

   for (ivvinstseq ivvs=vvist.begin();ivvs!=vvist.end();ivvs++){
      for (ivinstseq ivs=ivvs->begin();ivs!=ivvs->end();ivs++){
         os << *ivs;
      }
      os << "\n";
   }

   return os;
}
   

   void
importTSS(vTSS & vt,ifstream & file)
{
   back_insert_iterator<vTSS> dest(vt);
   copy(iisTSS(file),iisTSS(),dest);
//   cout << vt.size() << endl;
}
   
   vint
alignedtomap(string seq)
{
   string::const_iterator istr;

   int i=0;
   vint map;
   for (istr=seq.begin();istr!=seq.end();istr++){
      map.push_back(i);
      if (*istr!='-'){
         i++;
      }
   }  
   return map;
}

   vint
alignedtorevmap(string seq)
{
   string::const_iterator istr;

   int i=0;
   vint map;
   for (istr=seq.begin();istr!=seq.end();istr++){
      if (*istr!='-'){
         map.push_back(i);
      }
      i++;
   }  
   return map;
}

   vint
reversecomp(vint & istr)
{
   vint revistr;
   for (vint::reverse_iterator is=istr.rbegin();is!=istr.rend();is++){
      if (*is==0){
         revistr.push_back(1);
      }
      else if (*is==1){
         revistr.push_back(0);
      }
      else if (*is==2){
         revistr.push_back(3);
      }
      else if (*is==3){
         revistr.push_back(2);
      }
      else revistr.push_back(*is);
   }
   return revistr;
}

   vvd
reversecomp(vvd & matrice)
{
   vvd matrev;
   for (vvd::reverse_iterator imat=matrice.rbegin();imat!=matrice.rend();imat++){
      vd line;
      line.push_back((*imat)[1]);
      line.push_back((*imat)[0]);
      line.push_back((*imat)[3]);
      line.push_back((*imat)[2]);
      matrev.push_back(line);
   }
   return matrev;
}

//used for regs2000.fa
   vseq
loadsequences(ifstream & list)
{
   vseq seqs;
   string id,seqstr;
   int i=0;
   string tmpstring;
   getline(list,tmpstring);
   while (!list.eof()){
      Sequence seq;
      string curline;
      seq.name=tmpstring;
      getline(list,curline);
      seq.seqs.push_back(curline);
      vint iseq=stringtoint(curline);
      seq.iseqs.push_back(iseq);
      seq.nbN=compN(iseq);
      seq.nbtb=curline.size()-seq.nbN;
      seqs.push_back(seq);
      i++;
      getline(list,tmpstring);
   }
   return seqs;
}

   vseq
loadsequencesints(ifstream & list)
{
   vseq seqs;
   string id,seqstr;
   int i=0;
   string tmpstring;
   getline(list,tmpstring);
   while (!list.eof()){
      Sequence seq;
      string curline;
      seq.name=tmpstring;
      getline(list,curline);
      seq.seqsrealigned.push_back(curline);
      seq.imaps.push_back(alignedtomap(curline));
      seq.imapsinv.push_back(alignedtorevmap(curline));
      for (unsigned int j=0;j<3;j++){
         getline(list,curline);
         getline(list,curline);
         seq.seqsrealigned.push_back(curline);
         seq.imaps.push_back(alignedtomap(curline));
         seq.imapsinv.push_back(alignedtorevmap(curline));
      }
      for (unsigned int j=0;j<4;j++){
         getline(list,curline);
         getline(list,curline);
         seq.seqs.push_back(curline);
         seq.iseqs.push_back(stringtoint(curline));
      }
      seq.nbN=compN(seq.iseqs[0]);
      seq.nbtb=seq.iseqs[0].size()-seq.nbN;
      seqs.push_back(seq);
//      cout << i << endl;
      i++;
//      cout << seq.name << endl;
      getline(list,tmpstring);
   }
   return seqs;
}

Sequence loadseqconserv(string & filename)
{
	ifstream fseq;
	fseq.open(filename.c_str());
	Sequence seq;
	fseq >> seq;
	seq.name= filename.c_str(); // for display purpose
	fseq.close();
	
	return seq;
}
   
   vseq
loadsequencesconserv(ifstream & list)
{
   vstring seqsfile;
   back_insert_iterator<vstring> dest(seqsfile);
   copy(iisstring(list),iisstring(),dest);

   vseq seqs;

	//cout << "inputs" << endl;

   unsigned int counter=0;
   for (ivstring is=seqsfile.begin();is!=seqsfile.end();is++){
//      cout << *is << endl;
//	   if (counter>300) break;
	  ifstream fseq;
      fseq.open((*is).c_str());
      Sequence seq;
      // get the filename
      seq.finame=*is;
//      cout << seq.finame << endl;
      fseq >> seq;
      //we capture the sequence name, between final / and .fa
      string tname;
      tname=*is;
      tname.erase(tname.begin(),tname.begin()+tname.rfind('/')+1);
      tname.erase(tname.begin()+tname.rfind('.'),tname.end());
      int found=tname.find('_');
      if (found!=string::npos){
         //if there are mpre than two _
         if (tname.find('_',found+1)!=string::npos){
            //stop
            found=tname.find_last_of('_');
            tname=tname.substr(0,found);
            //start
            found=tname.find_last_of('_');
            tname=tname.substr(0,found);
            //chrom
            found=tname.find_last_of('_');
            tname=tname.substr(0,found);
         }
      }
      seq.name=tname;//(*is).c_str(); // for display purpose
      fseq.close();
      seqs.push_back(seq);
	   counter++;
   }
   return seqs;
}

   vseq
loadsequencesconservonly(ifstream & list)
{
   vstring seqsfile;
   back_insert_iterator<vstring> dest(seqsfile);
   copy(iisstring(list),iisstring(),dest);

   vseq seqs;

	//cout << "inputs" << endl;

   for (ivstring is=seqsfile.begin();is!=seqsfile.end();is++){
      //      cout << *is << endl;
      //	   if (counter>300) break;
      ifstream fseq;
      fseq.open((*is).c_str());
      Sequence seq;
      fseq >> seq;
      vint spe;
      int i=0;
      for (ivint iv=seq.species.begin();iv!=seq.species.end();iv++){
         if (*iv==1) spe.push_back(i);
         i++;
      }
      if (iscons(spe)){ 
         //we capture the sequence name, between final / and .fa
         string tname;
         tname=*is;
         tname.erase(tname.begin(),tname.begin()+tname.rfind('/')+1);
         tname.erase(tname.begin()+tname.rfind('.'),tname.end());
         seq.name=tname;//(*is).c_str(); // for display purpose
         fseq.close();
         seqs.push_back(seq);
      }
   }
   return seqs;
}
   
   vcoord
loadcoordconservwstrand(ifstream & list)
{
   vcoord vcds;
   string dum;
   getline(list,dum);
   while(!list.eof()){
      Coordinate coord;
      stringstream line(dum);
      string chrom;
      line>>chrom;
      coord.chrom=intfromchrom(chrom);
      line>>coord.start;
      line>>coord.stop;
      line>>coord.strand;
      line>>coord.name;
      if (coord.chrom!=-1){
         vcds.push_back(coord);
      }
      getline(list,dum);
   }
   return vcds;
}
   
vcoord
loadcoordfromTSS(ifstream & list)
{
   vcoord vcds;
   string dum;
   getline(list,dum);
   while(!list.eof()){
      Coordinate coord;
      stringstream line(dum);
      string chrom;
      line>>coord.start;
      coord.stop=coord.start;
      line>>dum;
      if (dum=="+") coord.strand=1;
      else coord.strand=-1;
      line>>coord.name;
      line>>chrom;
      coord.chrom=intfromchrom(chrom);
      if (coord.chrom!=-1){
         vcds.push_back(coord);
      }
      getline(list,dum);
   }
   return vcds;
}

int
loadcoordconserv(string folder, vcoord output)
{
   DIR *dp;
   struct dirent *ep;

   dp = opendir ("DATA_DIR");
   if (dp != NULL)
   {
      while (ep = readdir (dp))
         puts (ep->d_name);
      (void) closedir (dp);
   }
   else
      cerr << "Couldn't open the directory" << endl;

   return 0;
}

//Loads seqs from a folder containing .fa aligned sequences
   vseq
loadseqs(const char * folder)
{
   vseq seqs;
   
   DIR *dp;
   struct dirent *ep;

   dp = opendir ( folder );
   if (dp != NULL)
   {
      while (ep = readdir (dp)){
         string file = string(folder);
         file += "/";
         file += ep->d_name;
         if ( file.find(".fa") != string::npos ) {
            Sequence seq=loadseqconserv(file);
            // get rid of root path
            seq.name = string( ep->d_name );
            // get rid of .fa
            seq.name = seq.name.substr( 0 , seq.name.size() - 3 );
            seqs.push_back(seq);
         }
      }
      (void) closedir (dp);
   }
   else
      cerr << "Couldn't open the directory" << endl;

   return seqs;
}

//Loads fasta filenames from a folder
   vstring
loadfilenames(const char * folder)
{
   vstring filenames;
   
   DIR *dp;
   struct dirent *ep;

   dp = opendir ( folder );
   if (dp != NULL)
   {
      while (ep = readdir (dp)){
         string file = string(folder);
         file += "/";
         file += ep->d_name;
         if ( file.find(".fa") != string::npos ) {
            filenames.push_back( file );
         }
      }
      (void) closedir (dp);
   }
   else
      cerr << "Couldn't open the directory" << endl;

   return filenames;
}


   vcoord
loadcoordconserv(ifstream & list)
{
   vcoord vcds;
   string dum;
   getline(list,dum);
   while(!list.eof()){
      Coordinate coord;
      stringstream line(dum);
      string chrom;
      line>>chrom;
      coord.chrom=intfromchrom(chrom);
      line>>coord.start;
      line>>coord.stop;
      line>>dum;
      coord.name=dum;
         
      if (coord.chrom!=-1){
         vcds.push_back(coord);
      }
      getline(list,dum);
   }
   return vcds;
}


   double
scorefshift(vint & seqint, vvd &matrice)
{ 
   int deca(1);
   int dum(5);
   vint seqtmp=seqint;
   seqtmp.insert(seqtmp.begin(),deca,dum);
   seqtmp.insert(seqtmp.end(),deca,dum);

   double score(-100);
   for (ivint ibs=seqtmp.begin();ibs!=seqtmp.end()-width+1;ibs++){
      double sc=0;
      ivint ibs1=ibs;
      for (vvd::const_iterator imat=matrice.begin();imat!=matrice.end();imat++){
      int base=*ibs1;
         if (base==5){
            sc += 0;
         }
         else if (base==4){
            sc -= 100;
         }
         else sc+=(*imat)[base];
         ibs1++;
      }
      if (sc>score) score=sc;
   }
   return score;
}
   
// the matrice is a scoring PWM (matprec)
   double
scoref(vint::const_iterator &iseq, vvd &matrice)
{
   double sc=0;
   vint::const_iterator ibs=iseq;
   for (vvd::const_iterator imat=matrice.begin();imat!=matrice.end();imat++){
      const int base=*ibs;
      if (base == 4){
         sc -= 100;
      }
      else sc+=(*imat)[base];
      ibs++;
   }
   return sc;
}
   
   double
scoref(vint site, vvd &matrice)
{
   double sc=0;
   unsigned int pos=0;
   for (ivint iv=site.begin();iv!=site.end();iv++){
      const int base=*iv;
      if (base == 4){
         sc -= 100;
      }
      else sc+=matrice[pos][base];
      pos++;
   }
   return sc;
}
   
// here give a frequency matrix
   double
compprob(vint site, vvd &matrice)
{
   double sc=1;
   unsigned int pos=0;
   for (ivint iv=site.begin();iv!=site.end();iv++){
      const int base=*iv;
      if (base == 4){
         sc*=exp(-100);
      }
      else sc*=matrice[pos][base];
      pos++;
   }
   return sc;
}

unsigned int
shift(vint::const_iterator iseq,vvd & matrice, vint::const_iterator &seq_end, unsigned int extent)
{
   unsigned int shift=0;
//   if (args_info.motgen_given){
      double max=-100;
      for (unsigned int i=0;i<extent && iseq+i!=seq_end;i++){
         vint::const_iterator seqposi=iseq+i;
         double score=scoref(seqposi,matrice);
         if (score>max){
            max=score;
            shift=i;
         }
      }
//   } else {
//      double dist=1e6;
//      for (unsigned int i=0;i<extent && iseq+i!=seq_end;i++){
//         vint::const_iterator seqposi=iseq+i;
//         double score=scoref(seqposi,matrice);
//         if (fabs(2*i-extent)<dist && score> scorethr2){
//            dist=fabs(2*i-extent);
//            shift=i;
//         }
//      }
//   }

   return shift;
}

   int
basetoint(char base)
{
   int valueb=4;
   switch (base){
      case 'A':
         valueb=0;
         break;
      case 'a':
         valueb=0;
         break;
      case 'C':
         valueb=1;
         break;
      case 'c':
         valueb=1;
         break;
      case 'G':
         valueb=2;
         break;
      case 'g':
         valueb=2;
         break;
      case 'T':
         valueb=3;
         break;
      case 't':
         valueb=3;
         break;
   }
   return valueb;
}

Sequence
coordtoseq(Coordinate & coord)
{
//   cout << "coord to seq\n";
   Sequence seq;
   string dummystr;
   vint dummyvint;

   seq.seqsrealigned.insert(seq.seqsrealigned.begin(),nbspecies,dummystr);
   seq.imaps.insert(seq.imaps.begin(),nbspecies,dummyvint);
   seq.imapsinv.insert(seq.imapsinv.begin(),nbspecies,dummyvint);
   seq.species.insert(seq.species.begin(),nbspecies,0);
   seq.seqs.insert(seq.seqs.begin(),nbspecies,dummystr);
   seq.iseqs.insert(seq.iseqs.begin(),nbspecies,dummyvint);

   
   int test=0;
   for (ivcoord ivs=alignscoord.begin();ivs!=alignscoord.end();ivs++){
      Coordinate & ali=*ivs;
      if (ali.chrom==coord.chrom){
        if (ali.start<coord.start+1 && ali.stop>coord.start-1){
           
           // check if we need next alignment
           unsigned int isnextali=0;
           if (coord.stop>ali.stop){
              // check overlap between pieces of alignement, and that next alignment is sufficient
              if ( (ivs+1)->chrom == ivs->chrom && (ivs+1)->start < (ivs)->stop+2 && coord.stop < (ivs+1)->stop+1){
                 if (ivs != alignscoord.end()-1) isnextali=1;
              }
              else continue;
           }
   

           //if (species=="droso" && isnextali==1) continue; // pb with coordinates *** To be corrected ??

           Sequence trueali;
           ifstream fileseq(ali.name.c_str());
           fileseq >> trueali;
           fileseq.close();
           unsigned int pos=0;
           unsigned int truestart=0;
           unsigned int truestop=0;
           unsigned int counter=0;
           
           for (istring is=trueali.seqsrealigned[0].begin();is!=trueali.seqsrealigned[0].end();is++){
              if (*is!='-'){
                 if (coord.start-trueali.start==pos){
                    truestart=counter;
                 }
                 // if we need the next alignment, then we need to first stop at ali.stop
                 if (min(coord.stop,ali.stop)-trueali.start==pos){
                    truestop=counter;
                    break;
                 }
                 pos++;
              }
              counter++;
           }
           for (unsigned int i=0;i<nbspecies;i++){
              if (trueali.species[i]){
//                 cout << trueali.seqsrealigned[i] << "\n";
                 string seqrea=trueali.seqsrealigned[i].substr(truestart,truestop-truestart+1);
                 seq.seqsrealigned[i]=seqrea;
                 if (!isnextali){
                    seq.imaps[i]=alignedtomap(seqrea);
                    seq.imapsinv[i]=alignedtorevmap(seqrea);
                    string seqwogap=remgaps(seqrea);
                    if (i==0){
                       //cout << seqwogap << "\n\n";
                    }
                    if (seqwogap.size()>extraction_cutoff){
                       seq.species[i]=1;
                    } else {
                       seq.species[i]=0;
                    }
                    seq.seqs[i]=seqwogap;
                    seq.iseqs[i]=stringtoint(seqwogap);
                 }
              } else {
                 if (isnextali){
                       // fill absent species with gaps in case it is present in the next alignment
                       string dum(truestop-truestart+1,'-');
                       seq.seqsrealigned[i]=dum;
                 } else {
                    seq.species[i]=0;
                 }
              }
           }

           if (isnextali){
              Coordinate & nextali=*(ivs+1);
              ifstream fileseq(nextali.name.c_str());
              fileseq >> trueali;
              fileseq.close();
              // previous position
              int prevpos=ali.stop-nextali.start; //usually prevpos is -1
              pos=0;
              truestart=0;
              truestop=0;
              counter=0;
              for (istring is=trueali.seqsrealigned[0].begin();is!=trueali.seqsrealigned[0].end();is++){
                 if (*is!='-'){
                    if (pos-(prevpos+1)==0){
                       truestart=counter;
                    }
                    if (coord.stop-trueali.start==pos){
                       truestop=counter;
                       break;
                    }
                    pos++;
                 }
                 counter++;
              }
              for (unsigned int i=0;i<nbspecies;i++){
                 string seqrea;
                 if (trueali.species[i]){
                    seqrea=trueali.seqsrealigned[i].substr(truestart,truestop-truestart+1);
                 } else {
                    string dum(truestop-truestart+1,'-');
                    seqrea=dum;
                 }
                 seq.seqsrealigned[i].append(seqrea);
                 seq.imaps[i]=alignedtomap(seq.seqsrealigned[i]);
                 seq.imapsinv[i]=alignedtorevmap(seq.seqsrealigned[i]);
                 string seqwogap=remgaps(seq.seqsrealigned[i]);

                 if (seqwogap.size()>extraction_cutoff){
                    seq.species[i]=1;
                 } else {
                    seq.species[i]=0;
                 }
                 seq.seqs[i]=seqwogap;
                 seq.iseqs[i]=stringtoint(seqwogap);
              }
           }

           seq.name=coord.name;
           seq.start=coord.start;
           seq.stop=coord.stop;
           seq.chrom=coord.chrom;
           seq.nbN=compN(seq.iseqs[0]);
           seq.nbtb=seq.iseqs[0].size()-seq.nbN;
           test=1;
           
           vint spe;
           int i=0;
           for (ivint iv=seq.species.begin();iv!=seq.species.end();iv++){
              if (*iv==1) spe.push_back(i);
              i++;
           }
           if (iscons(spe)==0) test=2;
        }
      }
   }
   if (test==0) {
      cout << coord.name << ": sequence not found or overlap\n";
      cout.flush();
   }
   if (test==2) {
 //     cout << coord.name << ": sequence not conserved\n";
 //     cout.flush();
   }

   return seq;
}

//!! Here spe is a vint containing the species number (eg 0 1 4 6) and not if they are present (1 or 0)
   bool
iscons(vint & spe)
{
   vint boolspe(nbspecies,0);
   for (ivint iv=spe.begin();iv!=spe.end();iv++){
      boolspe[*iv]=1;
   }
   
   int nbfr=0;
   
   if (species=="droso"){
      if (boolspe[5]) nbfr++; 
      if (boolspe[6] || boolspe[7]) nbfr++;
      if (boolspe[8]) nbfr++;
      if (boolspe[9] || boolspe[10] || boolspe[11]) nbfr++;
      if (nbfr>1) return true;
   }
   else if (species=="eutherian"){
      if (boolspe[2] || boolspe[3] || boolspe[4] || boolspe[5] || boolspe[6] || boolspe[7]) nbfr++; // primates
      if (boolspe[8]) nbfr++; // horse
      if (boolspe[9]) nbfr++; // dog
      if (boolspe[10]) nbfr++; // wild boar
      if (boolspe[11]) nbfr++; // cow
      if (nbfr>1) return true;
   }
   return false;	
}

Sequence::Sequence()
{
   name="";
   score=0;
   sign=0;
   cons=0;
   vint dum(nbmots_for_score,0);
   motis=dum;
   signpermot=dum;
}

Coordinate::Coordinate()
{
   name="";
   chrom=-1;
   start=0;
   stop=0;
   strand=0;
   check=0;
}

void
Sequence::instances2instancescons()
{
   vint dum(nbmots_for_score,0);
   motis=dum;
   vinstseq vdum;
   vvinstseq vvinstspe(nbspecies,vdum);
   for (ivinstseq ivs=instances.begin();ivs!=instances.end();ivs++){
      vvinstspe[ivs->species].push_back(*ivs);
   }

   for (ivinstseq ivi=vvinstspe[0].begin();ivi!=vvinstspe[0].end();ivi++){

      vint vspe;
      vinstseq vtmp;
      vtmp.push_back(*ivi);

      for (int i=1;i<nbspecies;i++){

         int pdist=neighbext+1;
         ivinstseq bestivc;
         
         for (ivinstseq ivc=vvinstspe[i].begin();ivc!=vvinstspe[i].end();ivc++){
            if (ivc->motindex==ivi->motindex){
               int dist=fabs(ivc->pos-ivi->pos);
               if (dist<=neighbext && dist<pdist){
                  bestivc=ivc;
                  pdist=dist;
               }
            }
         }

         if (pdist<=neighbext){
            vtmp.push_back(*bestivc);
            vspe.push_back(bestivc->species);
         }
      }

      if (iscons(vspe)){
         instancescons.push_back(vtmp);
         motis[ivi->motindex]++;
      }

   }

   return;
}

Instanceseq::Instanceseq(unsigned int moti,int s, unsigned int p,unsigned int tp,double sc,int spe,string n)
{
   motindex=moti;
   sens=s;
   pos=p;
   truepos=tp;
   species=spe;
   score=sc;
   motname=n;
}
      
