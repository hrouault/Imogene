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
#ifndef Sequence_H
#define Sequence_H

#include <string>

#include "vectortypes.hpp"


extern vint lengthchrom;

class TSS
{
   public:
      int chrom;
      int coord;
      int sens;
      string gene;

      TSS(int position,char dir,string chr,string gname);
      TSS();
};

typedef vector<TSS> vTSS;
typedef vTSS::iterator ivTSS;
typedef istream_iterator<TSS> iisTSS;

void importTSS(vTSS & vt,ifstream & file);

istream & operator >>(istream &is,vTSS &vt);

class Instanceseq
{
   public:
      int motindex;
      string motname;
      int sens;
      int pos;//w/ gaps, and imaps[spe][pos] is pos on iseq
      int truepos;//wo gaps
      double score;
      int species;
      int shift; //shift compare to ref species

      string seq;
      vint iseq;

      Instanceseq(unsigned int moti,int s, unsigned int p,unsigned int tp,double sc,int spe,string n);
};


typedef vector<Instanceseq> vinstseq;
typedef vinstseq::iterator ivinstseq;
typedef vector<vinstseq> vvinstseq;
typedef vvinstseq::iterator ivvinstseq;

ostream & operator <<(ostream &os,Instanceseq &ist);
ostream & operator <<(ostream &os, vinstseq &vist);
ostream & operator <<(ostream &os, vvinstseq &vvist);

class Sequence
{
   public:
      string name;
      vstring sitenames;
      vstring speciesname;//drosoname
      string finame;
      int chrom;
      int start;
      int stop;
      vstring seqs;
      vstring seqsrealigned;
      vvint iseqs;
      vvint imaps; //relative pos w/o gaps (eg A-T, T is 2nd). i.e, given an aligned sequence (w/ gaps), gives position in iseq (w/o gaps)
      vvint imapsinv; //relative pos including gaps (eg A-T, T is 3rd). i.e, given an iseq (w/o gaps), gives position in an aligned seq (w/ gaps)
      vint species;//drosos;
      int nmot; // nb mots NON CONS
      int nbmot; // nb mots CONS
      vinstseq instances;
      vvinstseq instancescons;
      int nbN;
      int nbtb;
      // This is for score use
      double score;
      int sign;
      vint motis; //number of conserved motifs per mot index
      vint signpermot; //at a given cut-off, is 1 if sequence is >0, -1 else.

      bool cons; //0 if not cons, 1 if cons
      Sequence();
      void instances2instancescons();
};

typedef vector<Sequence> vseq;
typedef vseq::iterator ivseq;
typedef vseq::const_iterator civseq;
istream & operator >>(istream &is,Sequence & seq);
bool operator<(const Sequence & seqscr1,const Sequence & seqscr2);
ostream& operator <<(ostream &os,const vseq &vs);
ostream& operator <<(ostream &os,const Sequence & s);

class Chromosome
{
   public:
      string name;
      string seq;
};

typedef vector<Chromosome> vchrom;
typedef vchrom::iterator ivchrom;
typedef istream_iterator<Chromosome> iischrom;
bool operator<(const Chromosome & chr1,const Chromosome & chr2);

istream & operator >>(istream &is,Chromosome &chrom);

class Coordinate
{
   public:
      string name;
      int chrom;
      int start;
      int stop;

      int strand;

      vstring neargenes;
      int check;
      Coordinate();
};
typedef vector<Coordinate> vcoord;
typedef vcoord::iterator ivcoord;
typedef istream_iterator<Coordinate> iiscoord;


bool operator<(const Coordinate & coord1,const Coordinate & coord2);
istream & operator >>(istream &is,Coordinate &coord);
ostream & operator <<(ostream &os,Coordinate &coord);


int compN(vint & bs);

string vinttostring(vint & iseq);
string inttostring(int iseq);

vint stringtoint(string & seq);

string remgaps(string & seq);

vint alignedtomap(string seq);

vint alignedtorevmap(string seq);

vint reversecomp(vint & istr);

vvd reversecomp(vvd & matrice);
   
vseq loadseqs(const char * folder);

vstring loadfilenames(const char * folder);

vseq loadsequences(ifstream & list);

vseq loadsequencesints(ifstream & list);

Sequence loadseqconserv(string & filename);

vseq loadsequencesconserv(ifstream & list);

vseq loadsequencesconservonly(ifstream & list);

double scoref(vint::const_iterator &iseq, vvd &matrice);
double scoref(vint site, vvd &matrice);
double compprob(vint site, vvd &matrice);
double scorefshift(vint &iseq, vvd &matrice);

unsigned int shift(vint::const_iterator iseq,vvd & matrice, vint::const_iterator &seq_end, unsigned int extent);

int basetoint(char base);

string chromfromint(int chr);

int intfromchrom(string chrname);

vcoord loadcoordconservwstrand(ifstream & list);
vcoord loadcoordfromTSS(ifstream & list);
vcoord loadcoordconserv(ifstream & list);

Sequence coordtoseq(Coordinate & coord);

extern vcoord alignscoord;

bool iscons(vint & spe);

#endif // Sequence_H
