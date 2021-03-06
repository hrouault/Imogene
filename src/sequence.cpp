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
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <dirent.h>
#include <cstring>
#include <errno.h>

#include "vectortypes.hpp"
#include "sequence.hpp"
#include "tree.hpp"

using namespace std;

vcoord alignscoord;


bool
operator<(const Coordinate & coord1, const Coordinate & coord2)
{
    if (coord1.chrom == coord2.chrom) return (coord1.start <= coord2.start);
    else return coord1.chrom < coord2.chrom;
}

int
compN(vint & bs)
{
    int N = 0;
    for (vint::const_iterator ibs = bs.begin(); ibs != bs.end(); ibs++) {
        if (*ibs == 4) N++;
    }
    return N;
}

string
inttostring(int iseq)
{
    string seq;
    if (iseq == 0) {
        seq = "A";
    } else if (iseq == 1) {
        seq = "C";
    } else if (iseq == 2) {
        seq = "G";
    } else if (iseq == 3) {
        seq = "T";
    } else if (iseq == 4) {
        seq = "N";
    } else if (iseq == 5) {
        seq = "-";
    } else {
        seq = "?";
    }
    return seq;
}

string
vinttostring(vint & iseq)
{
    vint::const_iterator istr;
    string seq;
    for (istr = iseq.begin(); istr != iseq.end(); istr++) {
        seq += inttostring(*istr);
    }
    return seq;
}

vint
stringtoint(string & seq)
{
    string::const_iterator istr;
    vint iseq;
    for (istr = seq.begin(); istr != seq.end(); istr++) {
        if (*istr == 'A' || *istr == 'a') {
            iseq.push_back(0);
        } else if (*istr == 'C' || *istr == 'c') {
            iseq.push_back(1);
        } else if (*istr == 'G' || *istr == 'g') {
            iseq.push_back(2);
        } else if (*istr == 'T' || *istr == 't') {
            iseq.push_back(3);
        } else if (*istr == '-') {
            iseq.push_back(5);
        } else {
            iseq.push_back(4);
        }
    }
    return iseq;
}

// *** lowercase = REPEATS, N = CDS
vint
stringtointmaskrepeats(string & seq)
{
    string::const_iterator istr;
    vint iseq;
    for (istr = seq.begin(); istr != seq.end(); istr++) {
        if (*istr == 'A') {
            iseq.push_back(0);
        } else if (*istr == 'C') {
            iseq.push_back(1);
        } else if (*istr == 'G') {
            iseq.push_back(2);
        } else if (*istr == 'T') {
            iseq.push_back(3);
        } else if (*istr == '-') {
            iseq.push_back(5);
        } else {
            iseq.push_back(4);
        }
    }
    return iseq;
}

void
maskrepeats(vseq & seqs)
{
    for (ivseq ivs = seqs.begin(); ivs != seqs.end(); ivs++) {
        for (unsigned int i = 0; i < nbspecies; i++) {
            if (ivs->species[i]) {
                ivs->iseqs[i] = stringtointmaskrepeats(ivs->seqs[i]);
            }
        }
    }
    return;
}

string
remgaps(string & seq)
{
    string::const_iterator istr;
    string oseq;
    for (istr = seq.begin(); istr != seq.end(); istr++) {
        if (*istr != '-') {
            oseq.push_back(*istr);
        }
    }
    return oseq;
}

TSS::TSS()
{
}

TSS::TSS(int position, char dir, string chr, string gname, string gnamevar)
{
    if (chr == "2L") {
        chrom = 0;
    } else if (chr == "2R") {
        chrom = 1;
    } else if (chr == "3L") {
        chrom = 2;
    } else if (chr == "3R") {
        chrom = 3;
    } else if (chr == "4") {
        chrom = 4;
    } else if (chr == "X") {
        chrom = 5;
    }
    coord = position;
    if (dir == '+') {
        sens = 1;
    } else sens = -1;
    gene = gname;
    genevar = gnamevar;
}

string
chromfromint(int chr)
{
    if (species == "droso") {
        if (chr == 0) {
            return "2L";
        } else if (chr == 1) {
            return "2R";
        } else if (chr == 2) {
            return "3L";
        } else if (chr == 3) {
            return "3R";
        } else if (chr == 4) {
            return "4";
        } else if (chr == 5) {
            return "X";
        }
    } else if (species == "eutherian") {
        if (chr == 0) {
            return "1";
        } else if (chr == 1) {
            return "2";
        } else if (chr == 2) {
            return "3";
        } else if (chr == 3) {
            return "4";
        } else if (chr == 4) {
            return "5";
        } else if (chr == 5) {
            return "6";
        } else if (chr == 6) {
            return "7";
        } else if (chr == 7) {
            return "8";
        } else if (chr == 8) {
            return "9";
        } else if (chr == 9) {
            return "10";
        } else if (chr == 10) {
            return "11";
        } else if (chr == 11) {
            return "12";
        } else if (chr == 12) {
            return "13";
        } else if (chr == 13) {
            return "14";
        } else if (chr == 14) {
            return "15";
        } else if (chr == 15) {
            return "16";
        } else if (chr == 16) {
            return "17";
        } else if (chr == 17) {
            return "18";
        } else if (chr == 18) {
            return "19";
        } else if (chr == 19) {
            return "X";
        } else if (chr == 20) {
            return "Y";
        }
    }
    return "unknown";
}

int
intfromchrom(string chrname)
{
    if (species == "droso") {
        if (chrname == "2L" || chrname == "chr2L") {
            return 0;
        } else if (chrname == "2R" || chrname == "chr2R") {
            return 1;
        } else if (chrname == "3L" || chrname == "chr3L") {
            return 2;
        } else if (chrname == "3R" || chrname == "chr3R") {
            return 3;
        } else if (chrname == "4" || chrname == "chr4") {
            return 4;
        } else if (chrname == "X" || chrname == "chrX") {
            return 5;
        }
    } else if (species == "eutherian") {
        if (chrname == "1" || chrname == "chr1") {
            return 0;
        } else if (chrname == "2" || chrname == "chr2") {
            return 1;
        } else if (chrname == "3" || chrname == "chr3") {
            return 2;
        } else if (chrname == "4" || chrname == "chr4") {
            return 3;
        } else if (chrname == "5" || chrname == "chr5") {
            return 4;
        } else if (chrname == "6" || chrname == "chr6") {
            return 5;
        } else if (chrname == "7" || chrname == "chr7") {
            return 6;
        } else if (chrname == "8" || chrname == "chr8") {
            return 7;
        } else if (chrname == "9" || chrname == "chr9") {
            return 8;
        } else if (chrname == "10" || chrname == "chr10") {
            return 9;
        } else if (chrname == "11" || chrname == "chr11") {
            return 10;
        } else if (chrname == "12" || chrname == "chr12") {
            return 11;
        } else if (chrname == "13" || chrname == "chr13") {
            return 12;
        } else if (chrname == "14" || chrname == "chr14") {
            return 13;
        } else if (chrname == "15" || chrname == "chr15") {
            return 14;
        } else if (chrname == "16" || chrname == "chr16") {
            return 15;
        } else if (chrname == "17" || chrname == "chr17") {
            return 16;
        } else if (chrname == "18" || chrname == "chr18") {
            return 17;
        } else if (chrname == "19" || chrname == "chr19") {
            return 18;
        } else if (chrname == "X" || chrname == "chrX") {
            return 19;
        } else if (chrname == "Y" || chrname == "chrY") {
            return 20;
        }
    }
    return -1;
}

istream &
operator >>(istream & is, TSS & tss)
{
    int pos;
    char dir;
    string gname;
    string gnamevar;
    string chrom;
    is >> pos;
    is >> dir;
    is >> chrom;
    is >> gname;
    if (species == "droso"){
      is >> gnamevar;
    } else {
       gnamevar = gname;
    }
    if (dir == '+') {
        tss.sens = 1;
    } else if (dir == '-') {
        tss.sens = -1;
    } else if (dir == '.') {
        tss.sens = 0;
    } else if (!is.eof()) {
        cout << "error during TSS importation!! pos : " << pos << "\n";
    }
    tss.chrom = intfromchrom(chrom);
    tss.coord = pos;
    tss.gene = gname;
    tss.genevar = gnamevar;
    return is;
}

void
importTSS(vTSS & vt, ifstream & file)
{
    back_insert_iterator<vTSS> dest(vt);
    copy(iisTSS(file), iisTSS(), dest);
    //   cout << vt.size() << endl;
}

istream &
operator >>(istream & is, Sequence & seq)
{
    string dummystr;
    vint dummyvint;
    seq.seqsrealigned.insert(seq.seqsrealigned.begin(), nbspecies, dummystr);
    seq.imaps.insert(seq.imaps.begin(), nbspecies, dummyvint);
    seq.imapsinv.insert(seq.imapsinv.begin(), nbspecies, dummyvint);
    seq.species.insert(seq.species.begin(), nbspecies, 0);
    seq.seqs.insert(seq.seqs.begin(), nbspecies, dummystr);
    seq.iseqs.insert(seq.iseqs.begin(), nbspecies, dummyvint);
    string fseqline;
    getline(is, fseqline);
    stringstream firstline(fseqline);
    string dros;
    firstline >> dros;
    string chrom;
    firstline >> chrom;
    seq.chrom = intfromchrom(chrom);
    //	cout << chrom << " ";
    firstline >> seq.start;
    // cout << seq.start << " ";
    firstline >> seq.stop;
    // cout << seq.stop << "\n";
    while (!is.eof()) {
        //cout << fseqline << endl;
        int numdro = speciestonum(fseqline.substr(1, 6));
        //      cout << numdro << "\n";
        getline(is, fseqline);
        seq.seqsrealigned[numdro] = fseqline;
        seq.imaps[numdro] = alignedtomap(fseqline);
        seq.imapsinv[numdro] = alignedtorevmap(fseqline);
        string seqwogap = remgaps(fseqline);
        seq.seqs[numdro] = seqwogap;
        seq.iseqs[numdro] = stringtoint(seqwogap);
        unsigned int nbtb = seq.iseqs[numdro].size() - compN(seq.iseqs[numdro]);
        if (nbtb > width + neighbext) { //>5){
            seq.species[numdro] = 1;
        } else {
            seq.species[numdro] = 0;
        }
        getline(is, fseqline);
    }
    seq.nbN = compN(seq.iseqs[0]);
    seq.nbtb = seq.iseqs[0].size() - seq.nbN;
    return is;
}

// *** Comply with bed4 standard
istream &
operator >>(istream & is, Coordinate & coord)
{
    string chromname;
    is >> chromname;
    if (chromname == "") return is;
    coord.chrom = intfromchrom(chromname);
    is >> coord.start;
    is >> coord.stop;
    is >> coord.name;
    return is;
}

ostream &
operator <<(ostream & os, const Sequence & s)
{
    os << s.name << "\t";
    os << "chrom: " << chromfromint(s.chrom) << "\t";
    os << "start: " << s.start << "\t";
    os << "stop: " << s.stop << "\t";
    os << "length: " << s.iseqs[0].size() << "\t";
    os << endl;
    return os;
}

ostream &
operator <<(ostream & os, const vseq & vs)
{
    for (civseq ivs = vs.begin(); ivs != vs.end(); ivs++) {
        os << *ivs;
    }
    return os;
}

ostream &
operator <<(ostream & os, Coordinate & coord)
{
    os << coord.name << "\t";
    os << chromfromint(coord.chrom) << "\t";
    os << coord.start << "\t";
    os << coord.stop << "\n";
    //   os << chromfromint(coord.chrom) << "\t";
    //   os << coord.start << "\t";
    //   os << coord.stop << "\t";
    //   os << coord.name << "\n";
    return os;
}

ostream &
operator <<(ostream & os, Instanceseq & ist)
{
    // os << "Motif #" << ist.motindex+1 << "\t";
    os <<  ist.motname << "\t";
    os <<  ist.seq << "\t";
    os << "sens " << ist.sens << "\t";
    os << "pos (w/gaps) " << ist.pos << "\t";
    os << "pos " << ist.truepos << "\t";
    os << "score " << ist.score << "\t";
    os << "species " << ist.species ;
    os << " (" << numtospecies(ist.species) << ")" << "\t";
    os << "\n";
    return os;
}

ostream &
operator <<(ostream & os, vinstseq & vist)
{
    for (ivinstseq ivs = vist.begin(); ivs != vist.end(); ivs++) {
        os << *ivs;
    }
    return os;
}

vint
alignedtomap(string seq)
{
    string::const_iterator istr;
    int i = 0;
    vint map;
    for (istr = seq.begin(); istr != seq.end(); istr++) {
        map.push_back(i);
        if (*istr != '-') {
            i++;
        }
    }
    return map;
}

vint
alignedtorevmap(string seq)
{
    string::const_iterator istr;
    int i = 0;
    vint map;
    for (istr = seq.begin(); istr != seq.end(); istr++) {
        if (*istr != '-') {
            map.push_back(i);
        }
        i++;
    }
    return map;
}

vint
reversecomp(vint & istr)
{
    vint revistr;
    for (vint::reverse_iterator is = istr.rbegin(); is != istr.rend(); is++) {
        if (*is == 0) {
            revistr.push_back(3);
        } else if (*is == 1) {
            revistr.push_back(2);
        } else if (*is == 2) {
            revistr.push_back(1);
        } else if (*is == 3) {
            revistr.push_back(0);
        } else revistr.push_back(*is);
    }
    return revistr;
}

vvd
reversecomp(vvd & matrice)
{
    vvd matrev;
    for (vvd::reverse_iterator imat = matrice.rbegin(); imat != matrice.rend(); imat++) {
        vd line;
        line.push_back((*imat)[3]);
        line.push_back((*imat)[2]);
        line.push_back((*imat)[1]);
        line.push_back((*imat)[0]);
        matrev.push_back(line);
    }
    return matrev;
}

Sequence loadseqconserv(string & filename)
{
    ifstream fseq;
    fseq.open(filename.c_str());
    if (fseq.fail()) {
        cerr << "Sequence file opening failed: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    Sequence seq;
    fseq >> seq;
    seq.name = filename.c_str(); // for display purpose
    fseq.close();
    return seq;
}

//Loads seqs from a folder containing .fa aligned sequences
vseq
loadseqs(const char * folder)
{
    vseq seqs;
    DIR * dp;
    struct dirent * ep;
    dp = opendir(folder);
    if (dp != NULL) {
        while ((ep = readdir(dp))) {
            string file = string(folder);
            file += "/";
            file += ep->d_name;
            if (file.substr(file.size() - 3, file.size()) == ".fa") {
                Sequence seq = loadseqconserv(file);
                // get rid of root path
                seq.name = string(ep->d_name);
                // get rid of .fa
                seq.name = seq.name.substr(0 , seq.name.size() - 3);
                seqs.push_back(seq);
            }
        }
        (void) closedir(dp);
    } else {
        cerr << "Couldn't open the " << folder << " directory: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    return seqs;
}

//Loads fasta filenames from a folder
vstring
loadfilenames(const char * folder)
{
    vstring filenames;
    DIR * dp;
    struct dirent * ep;
    dp = opendir(folder);
    if (dp != NULL) {
        while ((ep = readdir(dp))) {
            string file = string(folder);
            file += "/";
            file += ep->d_name;
            if (file.substr(file.size() - 3, file.size()) == ".fa") {
                filenames.push_back(file);
            }
        }
        (void) closedir(dp);
    } else {
        cerr << "Couldn't open the " << folder << " directory: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    return filenames;
}

vcoord
loadcoordconserv(ifstream & list)
{
    vcoord vcds;
    string dum;
    getline(list, dum);
    while (!list.eof()) {
        Coordinate coord;
        stringstream line(dum);
        string chrom;
        line >> chrom;
        coord.chrom = intfromchrom(chrom);
        line >> coord.start;
        line >> coord.stop;
        line >> dum;
        coord.name = dum;
        if (coord.chrom != -1) {
            vcds.push_back(coord);
        }
        getline(list, dum);
    }
    return vcds;
}


// the matrix is a scoring PWM (matprec)
double
scoref(vint::const_iterator & iseq, vvd & matrice)
{
    double sc = 0;
    vint::const_iterator ibs = iseq;
    for (vvd::const_iterator imat = matrice.begin(); imat != matrice.end(); imat++) {
        const int base = *ibs;
        if (base == 4) {
            sc -= 100;
        } else sc += (*imat)[base];
        ibs++;
    }
    return sc;
}

double
scoref(vint site, vvd & matrice)
{
    double sc = 0;
    unsigned int pos = 0;
    for (ivint iv = site.begin(); iv != site.end(); iv++) {
        const int base = *iv;
        if (base == 4) {
            sc -= 100;
        } else sc += matrice[pos][base];
        pos++;
    }
    return sc;
}

unsigned int
shift(vint::const_iterator iseq, vvd & matrice, vint::const_iterator & seq_end, unsigned int extent)
{
    unsigned int shift = 0;
    //   if (args_info.motgen_given){
    double max = -100;
    for (unsigned int i = 0; i < extent && iseq + i != seq_end; i++) {
        vint::const_iterator seqposi = iseq + i;
        double score = scoref(seqposi, matrice);
        if (score > max) {
            max = score;
            shift = i;
        }
    }
    //   } else {
    //      double dist=1e6;
    //      for (unsigned int i=0;i<extent && iseq+i!=seq_end;i++){
    //         vint::const_iterator seqposi=iseq+i;
    //         double score=scoref(seqposi,matrice);
    //         if (fabs(2*i-extent)<dist && score> scorethr2){
    //            dist=fabs(2*i-extent);
    //            shift=i;
    //         }
    //      }
    //   }
    return shift;
}

string sequence_datapath;

Sequence
coordtoseq(Coordinate & coord)
{
    //   cout << "coord to seq\n";
    Sequence seq;
    string dummystr;
    vint dummyvint;
    seq.seqsrealigned.insert(seq.seqsrealigned.begin(), nbspecies, dummystr);
    seq.imaps.insert(seq.imaps.begin(), nbspecies, dummyvint);
    seq.imapsinv.insert(seq.imapsinv.begin(), nbspecies, dummyvint);
    seq.species.insert(seq.species.begin(), nbspecies, 0);
    seq.seqs.insert(seq.seqs.begin(), nbspecies, dummystr);
    seq.iseqs.insert(seq.iseqs.begin(), nbspecies, dummyvint);
    int test = 0;
    for (ivcoord ivs = alignscoord.begin(); ivs != alignscoord.end(); ivs++) {
        Coordinate & ali = *ivs;
        if (ali.chrom == coord.chrom) {
            if (ali.start < coord.start + 1 && ali.stop > coord.start - 1) {
                // check if we need next alignment
                unsigned int isnextali = 0;
                if (coord.stop > ali.stop) {
                    // check overlap between pieces of alignement, and that next alignment is sufficient
                    if ((ivs + 1)->chrom == ivs->chrom && (ivs + 1)->start < (ivs)->stop + 2 && coord.stop < (ivs + 1)->stop + 1) {
                        if (ivs != alignscoord.end() - 1) isnextali = 1;
                    } else continue;
                }
                //if (species=="droso" && isnextali==1) continue; // pb with coordinates *** To be corrected ??
                Sequence trueali;
                string alifn;
                if (species == "droso") {
                    alifn = sequence_datapath + "/droso/" + ali.name.c_str();
                } else if (species == "eutherian") {
                    alifn = sequence_datapath + "/eutherian/" + ali.name.c_str();
                }
                ifstream fileseq(alifn.c_str());
                if (fileseq.fail()) {
                    cerr << "Cannot open alignment sequence file for reading: " << strerror(errno) << endl;
                    exit(EXIT_FAILURE);
                }
                fileseq >> trueali;
                fileseq.close();
                int pos = 0;
                int truestart = 0;
                int truestop = 0;
                unsigned int counter = 0;
                for (istring is = trueali.seqsrealigned[0].begin(); is != trueali.seqsrealigned[0].end(); is++) {
                    if (*is != '-') {
                        if (coord.start - trueali.start == pos) {
                            truestart = counter;
                        }
                        // if we need the next alignment, then we need to first stop at ali.stop
                        if (min(coord.stop, ali.stop) - trueali.start == pos) {
                            truestop = counter;
                            break;
                        }
                        pos++;
                    }
                    counter++;
                }
                for (unsigned int i = 0; i < nbspecies; i++) {
                    if (trueali.species[i]) {
                        //                 cout << trueali.seqsrealigned[i] << "\n";
                        string seqrea = trueali.seqsrealigned[i].substr(truestart, truestop - truestart + 1);
                        seq.seqsrealigned[i] = seqrea;
                        if (!isnextali) {
                            seq.imaps[i] = alignedtomap(seqrea);
                            seq.imapsinv[i] = alignedtorevmap(seqrea);
                            string seqwogap = remgaps(seqrea);
                            if (i == 0) {
                                //cout << seqwogap << "\n\n";
                            }
                            if (seqwogap.size() > extraction_cutoff) {
                                seq.species[i] = 1;
                            } else {
                                seq.species[i] = 0;
                            }
                            seq.seqs[i] = seqwogap;
                            seq.iseqs[i] = stringtoint(seqwogap);
                        }
                    } else {
                        if (isnextali) {
                            // fill absent species with gaps in case it is present in the next alignment
                            string dum(truestop - truestart + 1, '-');
                            seq.seqsrealigned[i] = dum;
                        } else {
                            seq.species[i] = 0;
                        }
                    }
                }
                if (isnextali) {
                    Coordinate & nextali = *(ivs + 1);
                    string alifn;
                    if (species == "droso") {
                        alifn = sequence_datapath + "/droso/" + nextali.name.c_str();
                    } else if (species == "eutherian") {
                        alifn = sequence_datapath + "/eutherian/" + nextali.name.c_str();
                    }
                    ifstream fileseq2(alifn.c_str());
                    if (fileseq2.fail()) {
                        cerr << "Cannot open alignment sequence file for reading: " << strerror(errno) << endl;
                        exit(EXIT_FAILURE);
                    }
                    fileseq2 >> trueali;
                    fileseq2.close();
                    // previous position
                    int prevpos = ali.stop - nextali.start; //usually prevpos is -1
                    pos = 0;
                    truestart = 0;
                    truestop = 0;
                    counter = 0;
                    for (istring is = trueali.seqsrealigned[0].begin(); is != trueali.seqsrealigned[0].end(); is++) {
                        if (*is != '-') {
                            if (pos - (prevpos + 1) == 0) {
                                truestart = counter;
                            }
                            if (coord.stop - trueali.start == pos) {
                                truestop = counter;
                                break;
                            }
                            pos++;
                        }
                        counter++;
                    }
                    for (unsigned int i = 0; i < nbspecies; i++) {
                        string seqrea;
                        if (trueali.species[i]) {
                            seqrea = trueali.seqsrealigned[i].substr(truestart, truestop - truestart + 1);
                        } else {
                            string dum(truestop - truestart + 1, '-');
                            seqrea = dum;
                        }
                        seq.seqsrealigned[i].append(seqrea);
                        seq.imaps[i] = alignedtomap(seq.seqsrealigned[i]);
                        seq.imapsinv[i] = alignedtorevmap(seq.seqsrealigned[i]);
                        string seqwogap = remgaps(seq.seqsrealigned[i]);
                        if (seqwogap.size() > extraction_cutoff) {
                            seq.species[i] = 1;
                        } else {
                            seq.species[i] = 0;
                        }
                        seq.seqs[i] = seqwogap;
                        seq.iseqs[i] = stringtoint(seqwogap);
                    }
                }
                seq.name = coord.name;
                seq.start = coord.start;
                seq.stop = coord.stop;
                seq.chrom = coord.chrom;
                seq.nbN = compN(seq.iseqs[0]);
                seq.nbtb = seq.iseqs[0].size() - seq.nbN;
                test = 1;
                vint spe;
                int i = 0;
                for (ivint iv = seq.species.begin(); iv != seq.species.end(); iv++) {
                    if (*iv == 1) spe.push_back(i);
                    i++;
                }
                if (iscons(spe) == 0) test = 2;
            }
        }
    }
    if (test == 0) {
        cout << coord.name << ": sequence not found or overlap\n";
        cout.flush();
    }
    if (test == 2) {
        //     cout << coord.name << ": sequence not conserved\n";
        //     cout.flush();
    }
    return seq;
}

//!! Here spe is a vint containing the species number (eg 0 1 4 6) and not if they are present (1 or 0)
bool
iscons(vint & spe)
{
    vint boolspe(nbspecies, 0);
    for (ivint iv = spe.begin(); iv != spe.end(); iv++) {
        boolspe[*iv] = 1;
    }
    int nbfr = 0;
    if (species == "droso") {
        if (boolspe[5]) nbfr++;
        if (boolspe[6] || boolspe[7]) nbfr++;
        if (boolspe[8]) nbfr++;
        if (boolspe[9] || boolspe[10] || boolspe[11]) nbfr++;
        if (nbfr > 1) return true;
    } else if (species == "eutherian") {
        if (boolspe[2] || boolspe[3] || boolspe[4] || boolspe[5] || boolspe[6] || boolspe[7]) nbfr++; // primates
        if (boolspe[8]) nbfr++; // horse
        if (boolspe[9]) nbfr++; // dog
        if (boolspe[10]) nbfr++; // wild boar
        if (boolspe[11]) nbfr++; // cow
        if (nbfr > 1) return true;
    }
    return false;
}

Sequence::Sequence()
{
    name = "";
    chrom = -1;
    start = 0;
    stop = 0;
    nbtot = 0;
    nbcons = 0;
    nbN = 0;
    nbtb = 0;
}

Coordinate::Coordinate()
{
    name = "";
    chrom = -1;
    start = 0;
    stop = 0;
    strand = 0;
    check = 0;
}

void
Sequence::instances2instancescons()
{
    vint dum(nbmots_for_score, 0);
    vinstseq vdum;
    vvinstseq vvinstspe(nbspecies, vdum);
    for (ivinstseq ivs = instances.begin(); ivs != instances.end(); ivs++) {
        vvinstspe[ivs->species].push_back(*ivs);
    }
    for (ivinstseq ivi = vvinstspe[0].begin(); ivi != vvinstspe[0].end(); ivi++) {
        vint vspe;
        vinstseq vtmp;
        vtmp.push_back(*ivi);
        for (unsigned int i = 1; i < nbspecies; i++) {
            unsigned int pdist = neighbext + 1;
            ivinstseq bestivc;
            for (ivinstseq ivc = vvinstspe[i].begin(); ivc != vvinstspe[i].end(); ivc++) {
                if (ivc->motindex == ivi->motindex) {
                    unsigned int dist = abs(ivc->pos - ivi->pos);
                    if (dist <= neighbext && dist < pdist) {
                        bestivc = ivc;
                        pdist = dist;
                    }
                }
            }
            if (pdist <= neighbext) {
                vtmp.push_back(*bestivc);
                vspe.push_back(bestivc->species);
            }
        }
        if (iscons(vspe)) {
            instancescons.push_back(vtmp);
        }
    }
    return;
}

Instanceseq::Instanceseq(unsigned int moti, int s, unsigned int p, unsigned int tp, double sc, int spe, string n)
{
    motindex = moti;
    sens = s;
    pos = p;
    truepos = tp;
    species = spe;
    score = sc;
    motname = n;
}
