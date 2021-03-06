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

#include "const.hpp"

using namespace std;

const unsigned int nbiter = 3;
// const unsigned int distwidth = 20; // max number of instances for background
                                      // motif distribution
                                      // This variable is now declared in the
                                      // header file const.hpp for reasons
                                      // explained there. 
unsigned int width = 10;

unsigned int neighbext = 20;

double scorethr1 = 13;
double scorethr2 = 13;
double scorethr = 11.5;
double scorethrcons = 11.5;

double conca, conct, concc, concg;
int evolutionary_model;
string species;
unsigned int nbspecies;

unsigned int scanwidth = 1000;
unsigned int scanstep = 50;

unsigned int annotextent = 10000;

unsigned int nbmots_for_score = 5;

double alpha; // beta exponent for A,T
double beta; // beta exponent for C,G

unsigned int nbchrom;

unsigned int extraction_cutoff = 20;

bool progress = false;
