#ifndef Random_H
#define Random_H

#include <gsl/gsl_rng.h>

extern gsl_rng * gslran;

void rnginit();
void rngtest();


#endif // Random_H
