/*
 * =====================================================================================
 *
 *       Filename:  const.cpp
 *
 *    Description:  Constants used in the program
 *
 *        Version:  1.0
 *        Created:  02.07.2008 11:09:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  H. Rouault (HR), rouault@lps.ens.fr
 *        Company:  Laboratoire de physique statistique, ENS
 *
 * =====================================================================================
 */

#include <string>

using namespace std;


unsigned int width;

double scorethr1;
double scorethr2;
double scorethr;
double scorethrcons;

double conca,conct,concc,concg;
unsigned int evolutionary_model;
string species;
unsigned int nbspecies;

unsigned int neighbext=20;

unsigned int scanwidth=1000;
unsigned int scanstep=50;

unsigned int nbmots_for_score=20;

double alpha=0.176; // beta exponent for A,T
double beta=0.2/0.3*alpha; // beta exponent for C,G

unsigned int nbchrom;
