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
   // Shouldn't always use the same random numbers? ^^^ Agreed
   long seed=time(NULL) * getpid();
   gsl_rng_set(gslran,seed);
   // *** Output seed to a file if using a variable seed!!
   // see gsl_rng_fwrite()
}

// *** I erased rngtest(), is it ok? ^^^ No pb.
