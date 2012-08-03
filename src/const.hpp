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


extern unsigned int width; // motif width

// for genmot
extern const unsigned int nbiter; // max iterations for convergence

// I declare the following directly in the header file. It shouldn't be
// a prolbem since this is a const and so the multiply defined var problem
// doesn't hold here
const unsigned int distwidth = 20; // max number of instances for background motif distribution

extern unsigned int neighbext; // extent of conserved motifs search

extern string species; // droso or eutherian
extern int evolutionary_model; // 1 (felsen) or 2 (halpern)
extern unsigned int nbspecies;

extern double scorethr1; // for distinfo use
extern double scorethr2; //main threshold for scangen
extern double scorethr; // init threshold for scangen
extern double scorethrcons; // for aligned species

// background frequencies
extern double conca;
extern double conct;
extern double concc;
extern double concg;

// dirichlet exponents
extern double alpha; // beta exponent for A,T
extern double beta; // beta exponent for C,G

// for scangen
extern unsigned int nbchrom;
extern unsigned int scanwidth; // enhancers width
extern unsigned int annotextent; // extent of gene search for CRMs annotation
extern unsigned int nbmots_for_score;

// sequences that are too short are considered as irrelevant in the sequence extraction process
extern unsigned int extraction_cutoff;

// Show progress on the command line
extern bool progress;

#endif // Const_H
