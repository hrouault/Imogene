#ifndef Motif_H
#define Motif_H

#include <string>
#include <vector>
#include <gsl/gsl_vector.h>

#include "sequence.hpp"
#include "const.hpp"

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

      Motalign(unsigned int pos,Sequence & seq,Motif & mot,int sens);
      Motalign();
      bool iscons();
      void print();
      void mask();
};

typedef vector<Motalign> vma;
typedef vma::iterator ivma;

class Instance
{
   public:
      int motindex;
      int chrom;
      int coord;
      int sens;
      double score;
      string site;

      Instance(int chr,int pos, int sen,int moti);
      Instance(int chr,int pos, int sen,int moti,double sco,string s);
};

bool operator<(const Instance & inst1,const Instance & inst2);

typedef vector<Instance> vinst;
typedef vinst::iterator ivinst;
typedef vector<vinst> vvinst;

class TFBS
{

   public:
      vint site;
      
      double score; // score w current motif
      double prob; // probability as observed in sequences
      double probref,probini,probmix; // probability to observe such a motif under the initial PWM, the PWM with refinement, or the Kmeans mixture
      int num; // number of instances found in a given set of sequences
      int numrand; // number of instances found in a(n average of) random sample(s) of TFBS from the PWM.

      vd scores; // one score for each motif


      TFBS();
};

typedef vector<TFBS> vtfbs;
typedef vtfbs::iterator ivtfbs;
typedef vtfbs::const_iterator civtfbs;

bool operator<(const TFBS & bs1,const TFBS & bs2);
ostream & operator <<(ostream &os,const TFBS & bs);
ostream & operator <<(ostream &os,const vtfbs & vbs);

class Motif
{
   public:
      int pos;
      string seqinit;
      string bsinit;
      vvd matrice;
      vvd matprec;
      vvd matfreq;
      vvd matenergy;
      int nbmot;
      vvd matricerevcomp;
      vvd matprecrevcomp;
      vvd matwocons;
      vvd matwoconsrevcomp;
      double lambda;
      double lambdatrain;
      int ntrain;
      double pvalue;
      double meanpoiss;
       double meanval;
      int nbmatchback;
      unsigned int nbmatch;
      int distmot[distwidth];
      double scorepoiss;
      vma seqs;
      bool check;
      vvinst instances;
      
      Sequence refinstances; 
      vinst refinstances_short;
  
      string name;
      string id; //for jaspar matrices
      vint indexes; //for clustered motifs

      double motscorethr2;
      double motscorethr;
      double motscorethrcons;
      unsigned int motwidth;

      double optthr; //thr for optimal discernemnt
      double optauc; // best ROC area
      double optTP; // TP at optimal thr
      double  optFP; // FP at optimal thr
      unsigned int index;

      Motif();
      void matinit(double scth);
      void matinitforscanmots(Sequence & seq);
      void matinithamming(double scth,unsigned int numhamm);
      void compprec();
      void compprec_MCMC();
      void compprec_test();
      void comprefinstances(vseq & regs,unsigned int nspe);
      void comprefinstancescons(unsigned int nspe);
      void comprefmot();
      void pvaluecomp();
      void calclambdaposneg(vseq & vscore);
      void calclambda();
      void calclambdaback();
      void lambdacomp();
      void lambdacompcons();
      void display(ostream & streamfile);
      void displaywname(ostream & streamfile);
      bool isdiff();
      Motif copy();
      int nbmatchmat(const Sequence & seq);
      int nbmatchcons(Sequence & seq);
      double scorematchcons (Sequence & seq);
      int nbmatchconsnmask(Sequence & seq);
      int nbmatchconsbest1kb(Sequence & seq);
      int nbmatchnmask(Sequence & seq,unsigned int moti);
      int nbmatchwomask(Sequence & seq,unsigned int moti);
      void findinstancesnmask(Sequence &seq);
      void findinstances(Sequence &seq);
      void findinstances(vseq &vs);
      void matchmatincr(Sequence & seq);
      void calcmeanpoiss();
      void calcscorepoiss();
      void printmatrice(ostream & streamfile);
      int dispmots (Sequence & seq, int motindex);
      void corrprec();
      int statemot (Sequence & seq,int pos,int num, double & scoremot);
      int nbmatchnmaskforsvg (Sequence & seq,unsigned int moti);
      void setscorethr2meaninfo();
      vvvd correlations(); //returns width matrices of correlation
      void enumeratewthres(int n,vint & site,vtfbs & sites);
      vint drawsite(double scorethr);// draws a site with matfreq probs

      void testpwmcv(ofstream & testcv,ofstream & motmeldb);
      void testpwmcvtimer(ofstream & testcv);
      void testpseudocount(ofstream & outf);
};

typedef vector<Motif> vmot;
typedef vector<vmot> vvmot;
typedef vmot::iterator ivmot;
      
void enumerateallsitesncomp(int n,vint & site,vtfbs & sites,Motif & motini,Motif & refmot,vmot & motmix,vd & kprobs,ofstream & outf);
void enumeratesitesforkmeans(int n,vint & site,vtfbs & sites,Motif & motini,Motif & refmot,vvmot & vmotmix,vvd & vkprobs,ofstream & outf);

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
      int goodpheno;
      int discarded;
      double score;
      int totmots;

      GroupInstance(int sta,int sto,int chr);
      void compbestannot();
      void compscore(vmot & lmots,unsigned int nbmots_score);
      void compscoreweight(vmot & lmots,unsigned int nbmots_score);
      int distance(const GroupInstance & gi);
      void isdiscarded();
};

typedef vector<GroupInstance> vginst;
typedef vginst::iterator ivginst;
typedef vginst::const_iterator civginst;
typedef vector<vginst> vvginst;

bool operator<(const GroupInstance & ginst1,const GroupInstance & ginst2);


class Combination
{
   public:
      vint motis;
      unsigned int FP,TP;
      double FPR,TPR;
      int e10up,e10down,e10neg;
      int e18up,e18down,e18neg;
      double score;
      vstring posnames;
      vseq posseqs;
};

typedef vector<Combination> vcombi;
typedef vcombi::iterator ivcombi;
typedef vcombi::const_iterator civcombi;

bool operator<(const Combination & comb1,const Combination & comb2);
ostream &operator <<(ostream &os,const Combination &comb);
ostream &operator <<(ostream &os,const vcombi &vcomb);


vd colmean(unsigned int pos,Motif * mot);
vd colmean_test_window_for_mean(unsigned int pos,Motif * mot);
vvd colmean_RW(unsigned int pos,Motif * mot);
vd colmean_dirichlet(unsigned int pos,Motif * mot);
vd colmean_dirichlet_metropolis(unsigned int pos,Motif * mot);
vd colmean_t_test(unsigned int pos,Motif * mot);
vd colmean_test_numiter(unsigned int pos,Motif * mot);
vd colmean_test_numiter_w_dirichlet(unsigned int pos,Motif * mot);
vd colmean_test_numiter_w_dirichlet_metropolis(unsigned int pos,Motif * mot);
void colmean_test_autocorr_w_dirichlet(unsigned int pos,Motif * mot);
vd colmean_test_likelihood_w_dirichlet(unsigned int pos,Motif * mot);
vd colmean_test_rejection_w_dirichlet_metropolis(unsigned int pos,Motif * mot);
vd colmean_test_rejection_n_autocorr(unsigned int pos,Motif * mot);
vd colmean_testCV(unsigned int pos,Motif * mot);
void colmean_testCV_w_dirichlet(unsigned int pos,Motif * mot);
vd colopti(unsigned int pos,Motif * mot);
double likelyhood(vd x, void *params);
double loglikelyhood(vd x, void *params);
double loglikely(const gsl_vector *v, void *params);
void loadmots ( const char * filename, vmot & mots );
void loadmotswnames ( const char * filename, vmot & mots );
void loadjaspardb ( vmot & mots );
   
void displayhist(vginst & vgi,ostream & ostr);
void displayhist_set(vginst & vgi, vstring geneset,ostream & ostr);

int compalpha();

void freqtolog(vvd & mat);
void countfreq(vvd & mat);
void countbases(Motif & mot,Sequence & bds);
void getmatrices(ifstream & file, Motif & mot);

Motif comprefmot(Motif & mot);
vvd mattofreq(vvd & mat);
void  matfreqdisp(vvd& matrice);
void displaymat(vvd & mat);

Motalign mahamming(unsigned int pos, Sequence & seq, Motif & mot,int sens,unsigned int numhamm);

#endif // Motif_H
