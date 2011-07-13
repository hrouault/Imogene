/*
 * =====================================================================================
 *
 *       Filename:  diff.cpp
 *
 *    Description:  Differences between lists
 *
 *        Version:  1.0
 *        Created:  04.07.2008 17:13:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  H. Rouault (HR), rouault@lps.ens.fr
 *        Company:  Laboratoire de physique statistique, ENS
 *
 * =====================================================================================
 */


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

using namespace std;

/*
 * =====================================================================================
 *        Class:  Enhancer
 *  Description:  
 * =====================================================================================
 */
class Enhancer
{
   public:
      /* ====================  LIFECYCLE     ======================================= */
//      Enhancer ();                             /* constructor */

      /* ====================  OPERATORS     ======================================= */
      unsigned int distance(const Enhancer & enh);

      /* ====================  OPERATIONS    ======================================= */

      /* ====================  ACCESS        ======================================= */

      /* ====================  INQUIRY       ======================================= */

      /* ====================  DATA MEMBERS  ======================================= */
      string chrom;
      unsigned int start;
      unsigned int stop;
   protected:

   private:

}; /* -----  end of class  Enhancer  ----- */

typedef vector<Enhancer> venh;
typedef istream_iterator<Enhancer> iisenh;


istream &
operator >>(istream &is,Enhancer & enh)
{
   double dumscore;
   is >> dumscore;

//   cout << "dumscore : " << dumscore << endl;

   char dum[5000];
   is.getline(dum,100,':');
   enh.chrom=dum;
//   cout << "chrom : " << dum << endl;
   is.getline(dum,100,'.');
   string coord=dum;
//   cout << "start : " << coord << endl;
   stringstream sscoord(coord);
   sscoord >> enh.start;
   is.getline(dum,100,'.');
   is.getline(dum,100,' ');
   coord=dum;
   stringstream ssstop(coord);
   ssstop >> enh.stop;
//   cout << "stop : " << coord << endl;
   is.getline(dum,5000);
//   is.getline(enh.start,100,"..");
//   is.getline(enh.stop,100);

   return is;
}

unsigned int
Enhancer::distance(const Enhancer & enh)
{
   if (chrom!=enh.chrom) return 10000000;
   else return abs((int)start-(int)enh.start);
}

unsigned int
compcommon(venh & ve1,venh & ve2, unsigned int nbln)
{
   unsigned int dist=0;
   for (unsigned int i=0;i<nbln;i++){
      for (unsigned int j=0;j<nbln;j++){
         if (ve1[i].distance(ve2[j])<1000){
            dist++;
            break;
         }
      }
   }
   return dist;
}

   int
main ( int argc, char *argv[] )
{
   ifstream file1(argv[1]);
   ifstream file2(argv[2]);

   venh ve1;

   back_insert_iterator<venh> dest(ve1);
   copy(iisenh(file1),iisenh(),dest);

   venh ve2;
   back_insert_iterator<venh> dest2(ve2);
   copy(iisenh(file2),iisenh(),dest2);

   unsigned int size1=ve1.size();
   unsigned int size2=ve2.size();
   unsigned int sizeinf=size1;
   if (size2<size1) sizeinf=size2;

//   cout << sizeinf << endl;

   for (unsigned int i=0;i<sizeinf;i++){
      unsigned int nbcom=compcommon(ve1,ve2,i);
      cout << (double)compcommon(ve1,ve2,i)/(double)i << "\n";
   }
   
   return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */


