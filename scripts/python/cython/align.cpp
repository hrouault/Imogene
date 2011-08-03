// We use the different I/O operations
#include <iostream> // for screen I/O
#include <fstream> // for file I/O
#include <sstream> // for string I/O
// We use strings
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

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


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// Converts species name to int
int nametoint(string strName)
{
   int iname;
   if (strName == "mus_musculus") iname = 0;
   else if (strName == "pan_troglodytes") iname = 1;
   else if (strName == "pongo_abelii") iname = 2;
   else if (strName == "macaca_mulatta") iname = 3;
   else if (strName == "homo_sapiens") iname = 4;
   else if (strName == "rattus_norvegicus") iname = 5;
   else if (strName == "equus_caballus") iname = 6;
   else if (strName == "canis_familiaris") iname = 7;
   else if (strName == "bos_taurus") iname = 8;
   else if (strName == "sus_scrofa") iname = 9;
   else if (strName == "gorilla_gorilla") iname = 10;
   else if (strName == "callithrix_jacchus") iname = 11;
   else {
      cerr << strName << endl;
      cerr << "Error : Name not found" << endl;
      exit(1);
   }

   return iname;
}

string shortname(string strName)
{
	string name;
	if (strName == "mus_musculus") name = "MusMus";
	else if (strName == "pan_troglodytes") name = "PanTro";
	else if (strName == "pongo_abelii") name = "PonAbe";
	else if (strName == "macaca_mulatta") name = "MacMul";
	else if (strName == "homo_sapiens") name = "HomSap";
	else if (strName == "rattus_norvegicus") name = "RatNor";
	else if (strName == "equus_caballus") name = "EquCab";
	else if (strName == "canis_familiaris") name = "CanFam";
	else if (strName == "bos_taurus") name = "BosTau";
	else if (strName == "sus_scrofa") name = "SusScr";
   else if (strName == "gorilla_gorilla") name = "GorGor";
   else if (strName == "callithrix_jacchus") name = "CalJac";
	else {
      cerr << strName << endl;
		cerr << "Error : Name not found" << endl;
		exit(1);
	}
	
	return name;
}

// Converts int to species name
string inttoname(int iname)
{
	string strName("");
	if (iname == 0) strName = "mus_musculus";
	else if (iname == 1) strName = "pan_troglodytes";
	else if (iname == 2) strName = "pongo_abelii";
	else if (iname == 3 )strName = "macaca_mulatta";
	else if (iname == 4) strName = "homo_sapiens" ;
	else if (iname == 5) strName = "rattus_norvegicus";
	else if (iname == 6) strName = "equus_caballus";
	else if (iname == 7) strName = "canis_familiaris";
	else if (iname == 8) strName = "bos_taurus";
	else if (iname == 9) strName = "sus_scrofa";
   else if (iname == 10) strName = "gorilla_gorilla";
   else if (iname == 11) strName = "callithrix_jacchus";
	else 
	{
		cerr << "Error in inttoname : Out of bound" << endl;
		exit(1);
	};
	return strName;
};

// reset a vector of sequences to NULL
void SetNullVector(vseq& vSeq)
{
	ivseq iter;
	int iiter=0;
	for(iter = vSeq.begin(); iter != vSeq.end(); iter++)
	{
		(*iter).SetNull(inttoname(iiter));
		iiter++;
	};
};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// Outputs infos about different sequences
ostream& operator <<(ostream &os,vseq & vSeq)
{
	
	for (ivseq iter=vSeq.begin(); iter<vSeq.end(); iter++)
	{
		os << endl;
		os << (*iter).name << " " ;
		os << (*iter).chr << " " ;
		//os << (*iter).strChrNum << " " ;
		os << (*iter).start << " " ;
		os << (*iter).stop << " " ;
		os << (*iter).strand  << " " ;
		os << (*iter).length << " " ;
		os << "Col#" << (*iter).colnum << " " ;
	}
	os << endl;
	return os;
	
};

// Outputs infos about one sequence
ostream& operator <<(ostream &os,seq & cSeq)
{
	
	os << cSeq.name << " Chr." ;
	os << cSeq.chr << " ";//" #" ;
	//	os << cSeq.strChrNum << " " ;
	os << cSeq.start << ":" ;
	os << cSeq.stop << " " ;
	os << cSeq.strand  << " " ;
	os << cSeq.length << " Col." ;
	os << cSeq.colnum ;
	os << endl;
	return os;
};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// Erases the string strErase in the string strInput
string strClean(string& strInput, string strErase)
{				
	while (strInput.find(strErase) != string::npos)
	{
		strInput.erase(strInput.find(strErase),strErase.size());
	};
	
	return strInput;
};

// This takes an aligned string and returns the real sequence
 string rawtoseq(string strRaw)
{
	return strClean(strRaw, "-");
}

string whichchrom(string chr)
{
   if (chr == "1") return "chr1";
   else if (chr == "2") return "chr2";
   else if (chr == "3") return "chr3";
   else if (chr == "4") return "chr4";
   else if (chr == "5") return "chr5";
   else if (chr == "6") return "chr6";
   else if (chr == "7") return "chr7";
   else if (chr == "8") return "chr8";
   else if (chr == "9") return "chr9";
   else if (chr == "10") return "chr10";
   else if (chr == "11") return "chr11";
   else if (chr == "12") return "chr12";
   else if (chr == "13") return "chr13";
   else if (chr == "14") return "chr14";
   else if (chr == "15") return "chr15";
   else if (chr == "16") return "chr16";
   else if (chr == "17") return "chr17";
   else if (chr == "18") return "chr18";
   else if (chr == "19") return "chr19";
   else if (chr == "X") return "chrX";
   else if (chr == "Y") return "chrY";
   else return "others";
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void readinfo (ifstream& inf, vseq& seqs, string line) //get header info
{
	int iThisSpecies=0;
	int nColTemp =0 ;
   seq cSeq;

   for (int iii=0;iii<12;iii++) //init cSeq to NULL
   {
      cSeq.SetNull(inttoname(iii));
      seqs.push_back(cSeq);
   }

   while (line.find("DATA") == string::npos || line.find("SEQ") != string::npos)
   {
      if (line.find("ancestral_sequences") == string::npos && line.find("TREE") == string::npos)
      {
         string dum;
         istringstream ins;
         ins.str(line);
         cSeq.colnum = nColTemp;
         ins >> dum; // SEQ
         ins >> dum; //name 
         cSeq.name = dum;
         iThisSpecies = nametoint(cSeq.name);
         ins >> dum;
         cSeq.chr = dum ; 
         // transforms a string into an integer:
         ins >> dum;
         istringstream strin(dum);
         strin >> cSeq.start;
         ins >> dum;
         strin.clear();
         strin.str(dum);
         strin >> cSeq.stop;
         ins >> dum;
         strin.clear();
         strin.str(dum);
         strin >> cSeq.strand;
         ins >> dum;	
         strClean(dum,"(chr_length="); strClean(dum, ")");
         strin.clear();
         strin.str(dum);
         strin >> cSeq.length;
         cSeq.rawseq = "";
         seqs.at(iThisSpecies) = cSeq;
      };
      nColTemp++;
      getline(inf,line);
   }
}

string reversecomp(string rawseq)
{
   string revseq;
   for (string::iterator is=rawseq.end()-1; is>=rawseq.begin(); is--){
      if (*is=='A') revseq+='T';
      else if (*is=='T') revseq+='A';
      else if (*is=='C') revseq+='G';
      else if (*is=='G') revseq+='C';
      else revseq+=*is;
   }
   return revseq;
}

void getseqs(string filename, vvseq & falign)
{
   ifstream inf;
   inf.open(filename.c_str());
   string line;
   getline(inf,line);

   while(!inf.eof()){
      vseq seqs;
      while((line.find("SEQ") == string::npos)){
         getline(inf,line);
      } 
      readinfo(inf,seqs,line);
      //here we are on the 'DATA' line
      getline(inf,line);
      bool cond= (seqs[nametoint("mus_musculus")].rawseq != "...") && (whichchrom(seqs[nametoint("mus_musculus")].chr)!="others");
      //only keep mus alignments
      while (line.find("//") == string::npos){ 
         if (cond)
         {
            for (ivseq iseqs=seqs.begin(); iseqs<seqs.end(); iseqs++)
            {
               if ((*iseqs).chr != "undef")
               {
                  char base(line[(*iseqs).colnum]);
                  // The following is done in the main program
//                  if (base == 'a') base='A';
//                  if (base == 't') base='T';
//                  if (base == 'c') base='C';
//                  if (base == 'g') base='G';
//                  if (base == 'n') base='N';
                  (*iseqs).rawseq += base; 
               }
            }
         }
         getline(inf,line); 
      }

      if (cond)
      {
         falign.push_back(seqs);
         //cout << seqs[nametoint("Mus_musculus")].chr << "\t" << seqs[nametoint("Mus_musculus")].start << "\t" << seqs[nametoint("Mus_musculus")].stop <<"\t" << seqs[nametoint("Mus_musculus")].strand << endl;
         //cout << seqs[nametoint("Mus_musculus")].rawseq << endl;
//         seq s=seqs[nametoint("Mus_musculus")]; 
//         int strand;
//         strand=s.strand;   
//         if (strand==-1){
//            cout << s.chr << " "<< s.start << " " << s.stop << " " << endl;
//            cout << reversecomp(s.rawseq) << endl;
//            exit(9);
//         }
      }
      //we are ont the '//' line
      getline(inf,line);
   }
   //cout << endl;
   
   inf.close();
   return;
}

void writefiles(vvseq & falign)
{
   for (ivvseq ivv=falign.begin();ivv!=falign.end();ivv++){
		string filename("");
		//filename = "/home/santolin/files/mus/epo/align/";
		filename = "";
      filename += whichchrom((*ivv)[nametoint("mus_musculus")].chr);
      filename+= "/";
      stringstream os;
      os << "[ ! -d " << filename << " ] && mkdir " << filename << "";
      system(os.str().c_str());
		stringstream out;
      out << (*ivv)[nametoint("mus_musculus")].start;
      out << "-";
      out << (*ivv)[nametoint("mus_musculus")].stop;
      filename+= out.str();
      filename+= ".fa";

		ofstream outf(filename.c_str());
      double strand;
      strand=(*ivv)[nametoint("mus_musculus")].strand;
		for (ivseq iseqs = (*ivv).begin(); iseqs<(*ivv).end(); iseqs++)
		{
			seq s;
			s = *iseqs;
			if (s.rawseq != "...")
			{
            if (s.name=="mus_musculus"){
               if (whichchrom(s.chr)!="others")
               {
                  outf << ">" << shortname(s.name) << " chr" <<  s.chr << " " << s.start << " " << s.stop << endl;
               }
               else
               {
                  outf << ">" << shortname(s.name) << " " <<  s.chr << " " << s.start << " " << s.stop << endl;
               }
            }
            else {
				   outf << ">" << shortname(s.name) << endl;
            }
            if (strand==-1) outf << reversecomp(s.rawseq) << endl;
            else outf << s.rawseq << endl;
			}
		}
	}
}

// ---------------------------------------------------------------------------
// --------------------------------------------------------------------------

int main (int argc, char * const argv[]) 
{
	// Filenames of the alignment files
   ifstream inf("files.dat");
	string file;
   getline(inf,file); 

   while (!inf.eof())
   {
      cout << file << endl;
      //for each file we get all the alignments in a vvseq
      vvseq falign;
      cout << "Reading file..." << endl;
      getseqs(file,falign);
      cout << "Writing align files..." << endl;
      writefiles(falign);
      getline(inf,file);
   } 
	return(0);
}



