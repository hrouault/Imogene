#ifndef Matcons_H
#define Matcons_H

#include<vector>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_rng.h>

#include "vectortypes.hpp"
#include "motif.hpp"

using namespace std;

class noeud
{
   public:
      int esp1;
      int esp2;
      int noe;
      double prox1;
      double prox2;

      noeud(int e1,int e2, int n, double p1, double p2);
};

typedef vector<noeud> vnoe;
typedef vnoe::iterator ivnoe;

extern vnoe treedist;


class icprob
{
   public:
      double ic;
      double prob;
};
typedef vector<icprob> vic;
typedef vector<vic> vvic;
typedef vic::iterator ivic;
typedef vvic::iterator ivvic;

class column
{
   public:
      vint bases;

      column(Motalign & mots, int c);
      void complete();
      void compnoeud(noeud & noe);
};


typedef vector<gsl_matrix *> vpgslmat;
typedef vpgslmat::iterator ivpgslmat;

void phylotest(Motif & mot);
void evolvesite(Motif & mot);
void evolvesite(vmot & mots);
void evolvebase(Motif & mot);
void fitdistperbase(Motif & mot);
void fitdistpersite(Motif & mot);
void fitkmeanspersite(Motif & mot,const char * filename);
void fitdistkmeanspersite(Motif & mot,const char * filename);
void fitscorepersite(Motif & mot);
void inittreedist();

int speciestonum(string name);//species2num
string numtospecies(int num);//num2species

int func (double t, const double y[], double f[], void *params);
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);
double proba_fixation_rel(double ratio);

int instant_rates (const gsl_vector * w, gsl_matrix * rates);

#endif // Matcons_H
