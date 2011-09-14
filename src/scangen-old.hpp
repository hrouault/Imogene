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
#ifndef Motsearch_H 
#define Motsearch_H 

#include "imogene-genmot_cmdline.h"
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

using namespace std;

extern gengetopt_args_info args_info;


int basetoint(char base);
vint reversecomp(vint & istr);
vvd reversecomp(vvd & matrice);

class svg
{
public:
	int xsize;
	int ysize;
	int xoffset;
	int yoffset;
	int pos;
	
	svg();
};

class Resultdist
{
   public:
      int nbmotres;
      double mean;
      double min;
};


int nbmatchmat(Sequence & seq,double * matrice);
ostream& operator << (ostream &os,const vvd &matrice);
ostream& operator <<(ostream &os,const vint &bs);

vvd arraytomat();

extern vseq regtests;
extern vseq regints;
extern vvginst groupedinst;

double scoref(vint::const_iterator &iseq, vvd &matrice);
bool matchneighb(int pos, Sequence & seq,vvd & mat);

void outputresults(unsigned int nbmots_score);//,unsigned int nb_file);

void scanseqsforsvg(vseq & align);
void scanseqsforinstances(vseq & align,vmot & mots);
void scanseqsforinstances(vseq & align,Motif & mot);


void dispinit(ofstream & outf);
void dispseqwmots (Sequence & seq, vmot & mots, ofstream & outf);
void dispclose(ofstream & outf);

vseq loadseqs();
void disptex(vseq & seqs);
#endif /* Motsearch_H */
