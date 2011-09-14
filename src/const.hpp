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


#ifndef Const_H
#define Const_H

#include <string>


using namespace std;


extern unsigned int width;

const unsigned int distwidth=20;


const unsigned int coopext=250;

const unsigned int nbiter=3;

extern string species;
extern int evolutionary_model;
//const unsigned int nbspecies=4;
extern int nbspecies;

extern double scorethr1; // for distinfo use
extern double scorethr2; //main threshold for scangen
extern double scorethr; // init threshold for scangen
extern double scorethrcons; // for aligned species

extern unsigned int neighbext;

const unsigned int nbvalidated=10;//14; !! beware of the number of inputs !!

const double pvalthr=0.05;

const double valinf=0.081;
const double valsup=0.757;

extern double conca;//=0.268; // !! for the mouse
extern double conct;//=0.268;
extern double concc;//=0.232;
extern double concg;//=0.232;

/*
const double conca=0.3;
const double conct=0.3;
const double concc=0.2;
const double concg=0.2;
 */

extern double alpha; // beta exponent for A,T
extern double beta; // beta exponent for C,G

const double lookalikethr=40.0;

const unsigned int spacer=250;

extern unsigned int scanwidth;
extern unsigned int scanstep;

extern unsigned int nbchrom;
const unsigned int annotextent=10000;

extern unsigned int nbmots_for_score;

// sequences that are too short are considered as irrelevant in the sequence extraction process
extern unsigned int extraction_cutoff; 

#endif // Const_H
