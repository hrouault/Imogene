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

        noeud(int e1, int e2, int n, double p1, double p2);
};

typedef vector<noeud> vnoe;
typedef vnoe::iterator ivnoe;

extern vnoe treedist;

typedef vector<gsl_matrix *> vpgslmat;
typedef vpgslmat::iterator ivpgslmat;

void inittreedist();

int speciestonum(string name);//species2num
string numtospecies(int num);//num2species

int func(double t, const double y[], double f[], void * params);
int jac(double t, const double y[], double * dfdy, double dfdt[], void * params);
double proba_fixation_rel(double ratio);

int instant_rates(const gsl_vector * w, gsl_matrix * rates);

#endif // Matcons_H
