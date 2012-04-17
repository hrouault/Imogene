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
#ifndef Motif_H
#define Motif_H

#include <string>
#include <vector>
#include <gsl/gsl_vector.h>

#include "vectortypes.hpp"
#include "const.hpp"
#include "sequence.hpp"

class Motif;

class Motalign
{
    public:
        vvint alignseq;
        vivint matchespos;
        vint matches;
        ivint seq_start;
        ivint seq_stop;
        int strand;

        Motalign(unsigned int pos, Sequence & seq, Motif & mot, int sens);
        Motalign();
        bool iscons();
        void print();
        void mask();
};

typedef vector<Motalign> vma;
typedef vma::iterator ivma;
typedef vma::const_iterator icvma;

bool operator ==(const vma & seqs1, const vma & seqs2);

class Instance
{
    public:
        int motindex;
        int chrom;
        int coord;
        int sens;
        double score;
        string site;
        int isassigned; // non zero if has been asigned to a CRM

        Instance(int chr, int pos, int sen, int moti);
        Instance(int chr, int pos, int sen, int moti, double sco, string s);
};

bool operator<(const Instance & inst1, const Instance & inst2);

typedef vector<Instance> vinst;
typedef vinst::iterator ivinst;
typedef vinst::const_iterator civinst;
typedef vector<vinst> vvinst;

ostream & operator <<(ostream & os, const Instance & inst);

class Motif
{
    public:
        int pos;
        string bsinit;
        vvd matrice;
        vvd matprec;
        vvd matfreq;

        int nbtot;
        int nbcons;
        int nbconsback;
        double lambdaback;
        double lambdatrain;
        double pvalue;
        double meanpoiss;

        double scorepoiss;

        int * distmot;

        vvd matricerevcomp;
        vvd matprecrevcomp;
        vma seqs;
        bool check; // for distinfo
        vvinst instances;

        Sequence refinstances;
        vinst refinstances_short; // basic infos about motifs instances for scangen

        string name;
        unsigned int index;

        double motscorethr2;
        double motscorethr;
        double motscorethrcons;
        unsigned int motwidth;

        unsigned int tottest; // number of true bases in the background

        Motif();
        void matinit(double scth);
        void matinitforscanmots(Sequence & seq);
        void compprec();
        void compprec_MCMC();
        void compprec_inde();
        void pvaluecomp();
        void updatebacksites(Sequence & seq);
        void display(ostream & streamfile);
        int nbmatchcons(Sequence & seq);
        void findinstancesnmask(Sequence & seq);
        void findinstances(Sequence & seq);
        void findinstances(vseq & vs);
        void calcmeanpoiss();
        void calcscorepoiss();
        int dispmots(Sequence & seq, int motindex);
        void corrprec();
        int statemot(Sequence & seq, int pos, int num, double & scoremot);
        int nbmatchnmaskforsvg(Sequence & seq, unsigned int moti);
        void setscorethr2meaninfo();
};

typedef vector<Motif> vmot;
typedef vector<vmot> vvmot;
typedef vmot::iterator ivmot;

class GroupInstance
{
    public:
        int chrom;
        int start;
        int stop;
        vint nbmots;
        vinst instances;
        vTSS TSSs;
        TSS besttss;
        int distbesttss;
        int goodpheno;
        int discarded;
        double score;
        int totmots;

        GroupInstance();
        GroupInstance(int sta, int sto, int chr);
        void compbestannot();
        void compscore(vmot & lmots, unsigned int nbmots_score);
        void compscoreweight(vmot & lmots, unsigned int nbmots_score);
        int distance(const GroupInstance & gi);
        void isdiscarded();
};

typedef vector<GroupInstance> vginst;
typedef vginst::iterator ivginst;
typedef vginst::const_iterator civginst;
typedef vector<vginst> vvginst;

bool operator<(const GroupInstance & ginst1, const GroupInstance & ginst2);
ostream & operator <<(ostream & os, const GroupInstance & ginst);


class Combination
{
    public:
        vint motis;
        unsigned int FP, TP;
        double FPR, TPR;
        int e10up, e10down, e10neg;
        int e18up, e18down, e18neg;
        double score;
        vstring posnames;
        vseq posseqs;
};

typedef vector<Combination> vcombi;
typedef vcombi::iterator ivcombi;
typedef vcombi::const_iterator civcombi;

bool operator<(const Combination & comb1, const Combination & comb2);
ostream & operator <<(ostream & os, const Combination & comb);
ostream & operator <<(ostream & os, const vcombi & vcomb);


vd colmean(unsigned int pos, Motif * mot);
vd colopti(unsigned int pos, Motif * mot);
vd colinde(unsigned int pos, Motif * mot);
double likelyhood(vd x, void * params);
double loglikelyhood(vd x, void * params);
double loglikely(const gsl_vector * v, void * params);
void loadmots(const char * filename, vmot & mots);

void displayhist(vginst & vgi, ostream & ostr);
void displayhist_set(vginst & vgi, vstring geneset, ostream & ostr);

int compalpha();

void freqtolog(vvd & mat);
void countfreq(vvd & mat);
void countbases(Motif & mot, Sequence & bds);
void getmatrices(ifstream & file, Motif & mot);

Motif comprefmot(Motif & mot);
vvd mattofreq(vvd & mat);
void  matfreqdisp(vvd & matrice);
void displaymat(vvd & mat);

void scanseqforinstances(Sequence & seq, vmot & mots);
void scanseqsforinstances(vseq & align, Motif & mot);
void scanseqsforinstances(vseq & align, vmot & mots);
void scanseqforinstancesnmask(Sequence & seq, vmot & mots);
void scanseqsforinstancesnmask(vseq & align, vmot & mots);
void scanseqforconsinstances(Sequence & seq, vmot & mots);


#endif // Motif_H
