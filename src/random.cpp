#include <gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<iostream>
#include<fstream>

#include "vectortypes.hpp"
#include "const.hpp"
using namespace std;

gsl_rng * gslran;



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  rnginit
 *  Description:  Initialize the random number generator (inspired by gsl documentation)
 * =====================================================================================
 */
   void
rnginit()
{
   const gsl_rng_type * T;
   gsl_rng_env_setup();
   T = gsl_rng_default;
   gslran = gsl_rng_alloc (T);
   // Shouldn't always use the same random numbers?
   long seed=time(NULL) * getpid();
   gsl_rng_set(gslran,seed);
   // *** Output seed to a file if using a variable seed!!
   // see gsl_rng_fwrite()
}

// *** I erased rngtest(), is it ok?
