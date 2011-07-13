#ifndef Motsearch_H 
#define Motsearch_H 

#include "cmdline.h"
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
