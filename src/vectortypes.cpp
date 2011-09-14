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
#include <iostream>
#include <cmath>

#include "vectortypes.hpp"
#include "const.hpp"


   ostream&
operator <<(ostream &os,const vvd &matrice)
{
   for (vvd::const_iterator imat=matrice.begin();imat!=matrice.end();imat++){
      os << (*imat)[0] << ",";
   }
   os << "\t";
   for (vvd::const_iterator imat=matrice.begin();imat!=matrice.end();imat++){
      os << (*imat)[1] << ",";
   }
   os << "\t";
   for (vvd::const_iterator imat=matrice.begin();imat!=matrice.end();imat++){
      os << (*imat)[2] << ",";
   }
   os << "\t";
   for (vvd::const_iterator imat=matrice.begin();imat!=matrice.end();imat++){
      os << (*imat)[3] << ",";
   }
   os << "\t";
   return os;
}

   ostream&
operator <<(ostream &os,const vd &vect)
{
   for (civd ivect=vect.begin();ivect!=vect.end();ivect++){
      os << *ivect << ",";
   }
   os << "\n";
   return os;
}
   ostream&
operator <<(ostream &os,const vint &bs)
{
   for (vint::const_iterator ibs=bs.begin();ibs!=bs.end();ibs++){
      os << *ibs ;
   }
   return os;
}
   ostream&
operator <<(ostream &os,const vvint &vbs)
{
   for (civvint ivbs=vbs.begin();ivbs!=vbs.end();ivbs++){
//      for (civint ibs=ivbs->begin();ibs!=ivbs->end();ibs++){
//         os << *ibs ;
//      }
      os << *ivbs;
      os << "\n";
   }
   return os;
}

vd operator+(const vd & vec1,const vd & vec2) {
   vd vec;
   civd iv2=vec2.begin();
   for (civd iv1=vec1.begin();iv1!=vec1.end();iv1++){
      vec.push_back(*iv1+*iv2);
      iv2++;
   }
   return vec;
}

vvd operator +(const vvd & mat1,const vvd & mat2)
{
   vvd mat;
   civvd iv2=mat2.begin();
   for (civvd iv=mat1.begin();iv!=mat1.end();iv++){
      mat.push_back(*iv+*iv2);
      iv2++;
   }
   return mat;
}

vd operator -(const vd & col1,const vd & col2)
{
   vd col;
   civd iv2=col2.begin();
   for (civd iv=col1.begin();iv!=col1.end();iv++){
      col.push_back(*iv-*iv2);
      iv2++;
   }
   return col;
}

vvd operator -(const vvd & mat1,const vvd & mat2)
{
   vvd mat;
   civvd iv2=mat2.begin();
   for (civvd iv=mat1.begin();iv!=mat1.end();iv++){
      mat.push_back(*iv-*iv2);
      iv2++;
   }
   return mat;
}

vd abs(const vd & col1)
{
   vd col;
   for (civd iv=col1.begin();iv!=col1.end();iv++){
      col.push_back(fabs(*iv));
   }
   return col;
}

vvd abs(const vvd & mat1)
{
   vvd mat;
   for (civvd iv=mat1.begin();iv!=mat1.end();iv++){
      mat.push_back(abs(*iv));
   }
   return mat;
}

double sum(const vd & col1)
{
   double res=0;
   for (civd iv=col1.begin();iv!=col1.end();iv++){
      res+=*iv;
   }
   return res;
}

double sum(const vvd & mat1)
{
   double res=0;
   for (civvd iv=mat1.begin();iv!=mat1.end();iv++){
      res+=sum(*iv);
   }
   return res;
}

   istream&
operator >>(istream &is,vvd &matrice)
{
   char dummy;
   for (vvd::iterator imat=matrice.begin();imat!=matrice.end();imat++){
      is >> (*imat)[0] >> dummy;
   }
   for (vvd::iterator imat=matrice.begin();imat!=matrice.end();imat++){
      is >> (*imat)[1] >> dummy;
   }
   for (vvd::iterator imat=matrice.begin();imat!=matrice.end();imat++){
      is >> (*imat)[2] >> dummy;
   }
   for (vvd::iterator imat=matrice.begin();imat!=matrice.end();imat++){
      is >> (*imat)[3] >> dummy;
   }
   return is;
}

   istream&
operator >>(istream &is,int * distmot)
{
   char dummy;
   for (unsigned int i=0;i<distwidth-1;i++){
      is >> distmot[i] >> dummy;
   }
   is >> distmot[distwidth-1];
   return is;
}

   int
min(int n1, int n2)
{
   if (n1>n2) return n2;
   else return n1;
   return 0;
}
