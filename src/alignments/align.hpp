/*
 * =====================================================================================
 *
 *       Filename:  align.hpp
 *
 *    Description:  Take as input a file containing the list of absolute path to ensembl emf alignment files
 *                   and generate fasta formated alignments readable by Imogene 
 *
 *        Version:  1.0
 *        Created:  06.07.2011 15:10:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef CLASS_HPP
#define CLASS_HPP



// ---------------------------------------------------------------------------
// ---------------------------------CLASS-------------------------------------
// ---------------------------------------------------------------------------

// Info about a sequence
class seq
{
   public:
      string name, chr; // species, chromosome
      int start,stop,strand,length,colnum; // start, end, chromosome length, column number in the alignment file		
      string rawseq; // sequence from the alignment (with '-'s)

      void SetNull(string initname) // initializes to null
      {
         name = initname;
         chr = "undef";
         rawseq = "...";
         start = -1;
         stop = -1;
         strand = 0;
         length = -1;
         colnum = 100;
      };
};

typedef vector<seq> vseq;
typedef vseq::iterator ivseq;
typedef vector<vseq> vvseq;
typedef vvseq::iterator ivvseq;
typedef vector<int> vint;
typedef vint::iterator ivint;


int nametoint(string strName); // Converts species name to int
string shortname(string strName);
string inttoname(int iname); // Converts int to species name
void SetNullVector(vseq& vSeq); // reset a vector of sequences to NULL
ostream& operator <<(ostream &os,vseq & vSeq); // Outputs infos about different sequences
ostream& operator <<(ostream &os,seq & cSeq); // Outputs infos about one sequence
string strClean(string& strInput, string strErase); // Erases the string strErase in the string strInput
string rawtoseq(string strRaw); // This takes an aligned string and returns the real sequence
string whichchrom(string chr);
void readinfo (ifstream& inf, vseq& seqs, string line); //get header info
string reversecomp(string rawseq);
void getseqs(string filename, vvseq & falign);
void writefiles(vvseq & falign);
#endif // CLASS_HPP
