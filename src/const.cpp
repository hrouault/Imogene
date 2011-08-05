/*    
 * Copyright (C) 2003-2011 Hervé Rouault
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
 *
 * =====================================================================================
 *
 *       Filename:  const.cpp
 *
 *    Description:  Constants used in the program
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
