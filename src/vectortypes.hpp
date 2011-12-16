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
#ifndef Vectortypes_H
#define Vectortypes_H

#include <vector>
#include <iterator>

using namespace std;


typedef vector<int> vint;
typedef vector<vint> vvint;
typedef vector<double> vd;
typedef vector<vd> vvd;
typedef vector<vvd> vvvd;
typedef istream_iterator<double> iidouble;
typedef vector<string> vstring;
typedef string::iterator istring;

typedef vint::iterator ivint;
typedef vector<ivint> vivint;
typedef vvint::iterator ivvint;
typedef vd::iterator ivd;
typedef vvd::iterator ivvd;
typedef vint::const_iterator civint;
typedef vvint::const_iterator civvint;
typedef vd::const_iterator civd;
typedef vvd::const_iterator civvd;

ostream & operator <<(ostream & os, const vvd & matrice);
ostream & operator <<(ostream & os, const vint & bs);
ostream & operator <<(ostream & os, const vvint & vbs);
ostream & operator <<(ostream & os, const vd & vect);

typedef vstring::iterator ivstring;
typedef istream_iterator<string> iisstring;

vd operator+(const vd & vec1, const vd & vec2);
vvd operator +(const vvd & mat1, const vvd & mat2);
vd operator -(const vd & col1, const vd & col2);
vvd operator -(const vvd & mat1, const vvd & mat2);
vd abs(const vd & col1);
vvd abs(const vvd & mat1);
double sum(const vd & col1);
double sum(const vvd & mat1);
istream & operator >>(istream & is, vvd & matrice);
istream & operator >>(istream & is, int * distmot);
int min(int n1, int n2);

#endif // Vectortypes_H
