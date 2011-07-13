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

ostream& operator <<(ostream &os,const vvd &matrice);
ostream& operator <<(ostream &os,const vint &bs);
ostream& operator <<(ostream &os,const vvint &vbs);
ostream& operator <<(ostream &os,const vd &vect);

typedef vstring::iterator ivstring;
typedef istream_iterator<string> iisstring;

vd operator+(const vd & vec1,const vd & vec2);
vvd operator +(const vvd & mat1,const vvd & mat2);
vd operator -(const vd & col1,const vd & col2);
vvd operator -(const vvd & mat1,const vvd & mat2);
vd abs(const vd & col1);
vvd abs(const vvd & mat1);
double sum(const vd & col1);
double sum(const vvd & mat1);
istream& operator >>(istream &is,vvd &matrice);
istream& operator >>(istream &is,int * distmot);
int min(int n1, int n2);

#endif // Vectortypes_H
