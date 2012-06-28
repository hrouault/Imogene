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


class TSS
{
    public:
        int chrom;
        int coord;
        int sens;
        string gene; // usual name
        string genevar; // CG (flybase) or NM_ (UCSC) name

        TSS(int position, char dir, string chr, string gname, string gnamevar);
        TSS();
};

typedef vector<TSS> vTSS;
typedef vTSS::iterator ivTSS;
typedef vTSS::const_iterator civTSS;
typedef istream_iterator<TSS> iisTSS;

void importTSS(vTSS & vt, ifstream & file);

istream & operator >>(istream & is, vTSS & vt);

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

        Instanceseq(unsigned int moti, int s, unsigned int p, unsigned int tp, double sc, int spe, string n);
};

typedef vector<Instanceseq> vinstseq;
typedef vinstseq::iterator ivinstseq;
typedef vector<vinstseq> vvinstseq;
typedef vvinstseq::iterator ivvinstseq;

ostream & operator <<(ostream & os, Instanceseq & ist);
ostream & operator <<(ostream & os, vinstseq & vist);

class Sequence
{
    public:
        string name;
        string finame;
        int chrom;
        int start;
        int stop;
        vstring seqs;
        vstring seqsrealigned;
        vvint iseqs;
        vvint imaps; //relative pos w/o gaps (eg A-T, T is 2nd). i.e, given an aligned sequence (w/ gaps), gives position in iseq (w/o gaps)
        vvint imapsinv; //relative pos including gaps (eg A-T, T is 3rd). i.e, given an iseq (w/o gaps), gives position in an aligned seq (w/ gaps)
        vint species;//0 or 1 per species
        int nbtot; // formerly nmot nb mots NON CONS
        int nbcons; // formerly nbmot nb mots CONS
        vinstseq instances;
        vvinstseq instancescons;
        int nbN;
        unsigned int nbtb;

        Sequence();
        void instances2instancescons();
};

typedef vector<Sequence> vseq;
typedef vseq::iterator ivseq;
typedef vseq::const_iterator civseq;
istream & operator >>(istream & is, Sequence & seq);
bool operator<(const Sequence & seqscr1, const Sequence & seqscr2);
ostream & operator <<(ostream & os, const vseq & vs);
ostream & operator <<(ostream & os, const Sequence & s);


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

extern vcoord alignscoord;

bool operator<(const Coordinate & coord1, const Coordinate & coord2);
istream & operator >>(istream & is, Coordinate & coord);
ostream & operator <<(ostream & os, Coordinate & coord);


string vinttostring(vint & iseq);
string inttostring(int iseq);
vint stringtoint(string & seq);
vint stringtointmaskrepeats(string & seq);

void maskrepeats(vseq & seqs);
string remgaps(string & seq);
vint alignedtomap(string seq);
vint alignedtorevmap(string seq);
vint reversecomp(vint & istr);
vvd reversecomp(vvd & matrice);
int compN(vint & bs);

double scoref(vint::const_iterator & iseq, vvd & matrice);
double scorefhamming(vint site1, vint site2);
double scoref(vint site, vvd & matrice);
unsigned int shift(vint::const_iterator iseq, vvd & matrice, vint::const_iterator & seq_end, unsigned int extent);
unsigned int shifthamming(vint::const_iterator iseq, vint seqmel, vint::const_iterator &seq_end, unsigned int extent);

string chromfromint(int chr);
int intfromchrom(string chrname);

vseq loadseqs(const char * folder);
vstring loadfilenames(const char * folder);
Sequence loadseqconserv(string & filename);
vcoord loadcoordconserv(ifstream & list);
Sequence coordtoseq(Coordinate & coord);

bool iscons(vint & spe);

extern string sequence_datapath;

#endif // Sequence_H
