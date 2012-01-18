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

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <sys/stat.h>

using namespace std;

#include "display_cmdline.h"
#include "display.hpp"
#include "const.hpp"
#include "motif.hpp"
#include "sequence.hpp"
#include "tree.hpp"

display_args_info display_args;

void
dispweblogo(vmot & mots)
{
    unsigned int index = 1;
    for (ivmot ivm = mots.begin(); ivm != mots.end(); ivm++) {
        if (ivm->check) {
            stringstream ss;
            ss << "weblogo-display ";
            ss << "Motif";
            ss << index << " ";
            ss << concc << " ";
            string sep1 = "";
            for (ivvd ivv = ivm->matfreq.begin(); ivv != ivm->matfreq.end(); ivv++) {
                ss << sep1;
                string sep2 = "";
                for (ivd iv = ivv->begin(); iv != ivv->end(); iv++) {
                    ss << sep2;
                    ss << *iv;
                    sep2 = ",";
                }
                sep1 = ":";
            }
            if (display_args.png_given) ss << " --png ";
            if (display_args.pdf_given) ss << " --pdf";
            system(ss.str().c_str());
            index++;
        }
    }
}

void
dispjaspar(vmot & mots)
{
    unsigned int index = 1;
    for (ivmot ivm = mots.begin(); ivm != mots.end(); ivm++) {
        if (ivm->check) {
            cout << "Motif";
            cout << index << endl;

            stringstream lineA;
            stringstream lineC;
            stringstream lineG;
            stringstream lineT;
            lineA << "A [";
            lineC << "C [";
            lineG << "G [";
            lineT << "T [";

            for (ivvd ivv = ivm->matfreq.begin(); ivv != ivm->matfreq.end(); ivv++) {
                string sep2 = "";
                lineA << (int)((*ivv)[0] * 100) << " ";
                lineC << (int)((*ivv)[1] * 100) << " ";
                lineG << (int)((*ivv)[2] * 100) << " ";
                lineT << (int)((*ivv)[3] * 100) << " ";
            }
            lineA << "]\n";
            lineC << "]\n";
            lineG << "]\n";
            lineT << "]\n";

            cout << lineA.str();
            cout << lineC.str();
            cout << lineG.str();
            cout << lineT.str();
            index++;
        }
    }
}

void
texify(string & str)
{
    size_t found = str.find("_");
    while (found != string::npos) {
        str.insert(found, "\\");
        found = str.find("_", found + 3);
    }
    found = str.find("#");
    while (found != string::npos) {
        str.insert(found, "\\");
        found = str.find("#", found + 3);
    }
}

// for TFBS color
string
colfromint(int i)
{
    if (i == 0) return "red";
    if (i == 1) return "blue";
    if (i == 2) return "green";
    if (i == 3) return "yellow";
    if (i == 4) return "brown";
    return "black";
}

void
svginit(ofstream & svgfile, Sequence & seq)
{
    int nbgroup = seq.seqsrealigned[0].size() / 2000 + 1;
    int nbspetot = 0;
    for (unsigned int i = 0; i < nbspecies; i++) {
        if (seq.species[i]) nbspetot++;
    }
    svgfile << "<?xml version=\"1.0\" encoding=\"utf-8\" standalone=\"no\" ?>\n";
    svgfile << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
    svgfile << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
    svgfile << endl;
    svgfile << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" viewBox=\"0 0 816 " << (nbgroup * (1.5 * nbspetot + 2) + 1.5) * 12 << "\" version=\"1.1\">" << endl;
    svgfile << endl;
    svgfile << "    <style type=\"text/css\" >" << endl;
    svgfile << "      <![CDATA[" << endl;
    svgfile << endl;
    svgfile << "text {color:#222;background:#fff;font-family:\"Helvetica Neue\", Arial, Helvetica, sans-serif;}" << endl;
    //   svgfile << "text.score {font-size: 6;}" << endl;
    svgfile << ".sequence {" << endl;
    svgfile << "   stroke: black;" << endl;
    svgfile << "   stroke-width: 0.2em;" << endl;
    svgfile << "}" << endl;
    svgfile << endl;
    svgfile << ".gap {" << endl;
    svgfile << "   stroke: black;" << endl;
    svgfile << "   stroke-width: 0.05em;" << endl;
    svgfile << "}" << endl;
    svgfile << ".mot1 {" << endl;
    svgfile << "   stroke:  #730046;" << endl;
    svgfile << "}" << endl;
    svgfile << ".mot2 {" << endl;
    svgfile << "   stroke: #BFBB11;" << endl;
    svgfile << "}" << endl;
    svgfile << ".mot3 {" << endl;
    svgfile << "   stroke: #FFC200;" << endl;
    svgfile << "}" << endl;
    svgfile << ".mot4 {" << endl;
    svgfile << "   stroke: #E88801;" << endl;
    svgfile << "}" << endl;
    svgfile << ".conserved {" << endl;
    svgfile << "   stroke-width: 0.3em;" << endl;
    svgfile << "}" << endl;
    svgfile << ".notconserved {" << endl;
    svgfile << "   stroke-width: 0.1em;" << endl;
    svgfile << "}" << endl;
    svgfile << ".score {" << endl;
    svgfile << "   font-size: 0.7em;" << endl;
    svgfile << "}" << endl;
    svgfile << "" << endl;
    svgfile << "      ]]>" << endl;
    svgfile << "    </style>" << endl;
    svgfile << endl;
    svgfile << "<g font-size=\"12\">" << endl;
}

int
realindex(int spenb, Sequence & seq)
{
    unsigned int indspe = 0;
    for (unsigned int i = 0; i < spenb; i++) {
        if (seq.species[i]) indspe++;
    }
    return indspe;
}

void
svgdisplay(ofstream & svgfile, Sequence & seq)
{
    double yoffset = 1.5;
    double xbegin = 6;
    int nbspetot = 0;
    for (unsigned int i = 0; i < nbspecies; i++) {
        if (seq.species[i]) nbspetot++;
    }
    int pos = 0;
    int groupnb = 0;
    int sizeseq = seq.seqsrealigned[0].size();
    int stop = min(2000, sizeseq);
    int start = 0;
    double xscale = 0.03;
    while (start < sizeseq) {
        int linenb = 0;
        for (unsigned int i = 0; i < nbspecies; i++) {
            if (seq.species[i]) {
                double yline = yoffset + 1.5 * linenb;
                linenb++;
                string dro;
                dro = numtospecies(i);
                svgfile << "<text class=\"species\" x=\"0.2em\" y=\"" << yline + 0.4 << "em\">" << dro << "</text>\n";
                svgfile << "<line class=\"gap\" x1=\"" << xbegin << "em\" y1=\"" << yline << "em\" x2=\"" << (stop - start)*xscale + xbegin << "em\" y2=\"" << yline << "em\"/>\n";
                vvint coords;
                vint curcoord;
                int seqorali = 0;
                unsigned int p = 0;
                for (istring is = seq.seqsrealigned[i].begin() + start; is != seq.seqsrealigned[i].begin() + stop; is++) {
                    if (*is == '-') {
                        if (seqorali == 1) {
                            curcoord.push_back(p);
                            coords.push_back(curcoord);
                            curcoord.clear();
                            seqorali = 0;
                        }
                    } else {
                        if (seqorali == 0) {
                            curcoord.push_back(p);
                            seqorali = 1;
                        }
                        if (seqorali == 1 && is == seq.seqsrealigned[i].begin() + stop - 1) {
                            curcoord.push_back(p);
                            coords.push_back(curcoord);
                            curcoord.clear();
                        }
                    }
                    p++;
                }
                for (ivvint ivv = coords.begin(); ivv != coords.end(); ivv++) {
                    svgfile << "<line class=\"sequence\" x1=\"" << xbegin + 0.03 * (*ivv)[0] << \
                            "em\" y1=\"" << yline << "em\" x2=\"" << xbegin + 0.03 * (*ivv)[1] << "em\" y2=\"" << yline << "em\"/>\n";
                }
            }
        }
        for (ivinstseq ivi = seq.instances.begin(); ivi != seq.instances.end(); ivi++) {
            int moti = (*ivi).motindex;
            if (moti < 8 && (*ivi).pos < stop && (*ivi).pos > start) {
                double xmot = xbegin + xscale * ((*ivi).pos - start);
                double yline = yoffset + 1.5 * realindex((*ivi).species, seq);
                string width;
                string motclass;
                if ((*ivi).score > scorethr2)
                    motclass = "conserved";
                else
                    motclass = "notconserved";
                svgfile << "<line class=\"" + motclass + " mot" << moti + 1 << "\" x1=\"" << xmot << "em\" y1=\"" << yline - 0.5 << "em\" x2=\"" << xmot << "em\" y2=\"" << yline + 0.5 << "em\"/>\n";
                svgfile << "<text x=\"" << xmot + 0.2 << "em\" y=\"" << (yline + 0.8) << "em\"><tspan class=\"score\">" << fixed << setprecision(1) << (*ivi).score << "</tspan></text>\n";
            }
        }
        groupnb++;
        stop = min(2000 * (groupnb + 1), sizeseq);
        start = 2000 * groupnb;
        yoffset += nbspetot * 1.5 + 2;
    }
}

void
svgclose(ofstream & svgfile)
{
    svgfile << "</g>" << endl;
    svgfile << "</svg>\n";
}

void
scanseqforsvg(Sequence & align, vmot & mots)
{
    int xsize(0);
    int ysize(0);
    string filename("");
    filename += align.name;
    filename += ".svg";
    ofstream svgfile(filename.c_str());
    if (svgfile.fail()) {
        cerr << "Cannot open file for svg recording: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    svginit(svgfile, align);
    svgdisplay(svgfile, align);
    svgclose(svgfile);
}

// display TFBS on reference sequences
void
dispseqwmots(Sequence & seq, vmot & mots, ofstream & outf)
{
    //HEADER
    string name = numtospecies(0);
    string texname = seq.name;
    outf << "$>$" << name <<  "\t";
    outf << texname << "\t" ;
    outf << chromfromint(seq.chrom) << "\t";
    outf << seq.start << "\t" ;
    outf << seq.stop << " \\\\\n";
    outf << " \\\\\n";
    //DEFINING STATES, COLORS AND SCORES
    unsigned int sizeseq = seq.seqsrealigned[0].size();
    vint vdum(sizeseq, 0);
    vvint vvstate(nbspecies, vdum); // 0 nothing, 1 normal tfbs, 2 cons tfbs
    vvint vvcol(nbspecies, vdum); // 0 red, 1 blue, 2 green, 3 yellow
    vd vs(sizeseq, 0.);
    vvd scores(nbspecies, vs);
    //NOT CONS
    for (ivinstseq ivs = seq.instances.begin(); ivs != seq.instances.end(); ivs++) {
        int spe = ivs->species;
        int motwidth = mots[ivs->motindex].motwidth;
        int ipos = seq.imaps[spe][ivs->pos];
        int stop = seq.imapsinv[spe][ipos + motwidth - 1];
        for (int i = ivs->pos; i < stop + 1; i++) {
            vvstate[spe][i] = 1;
            vvcol[spe][i] = ivs->motindex;
            scores[spe][i] = ivs->score;
            //if (seq.seqsrealigned[spe][i]=='-') vvcol[spe][i]=-1;
        }
    }
    //CONS
    for (ivvinstseq ivvs = seq.instancescons.begin(); ivvs != seq.instancescons.end(); ivvs++) {
        for (ivinstseq ivs = (*ivvs).begin(); ivs != (*ivvs).end(); ivs++) {
            int spe = ivs->species;
            int motwidth = mots[ivs->motindex].motwidth;
            int ipos = seq.imaps[spe][ivs->pos];
            int stop = seq.imapsinv[spe][ipos + motwidth - 1];
            for (int i = ivs->pos; i < stop + 1; i++) {
                vvstate[spe][i] = 2;
                scores[spe][i] = ivs->score;
            }
        }
    }
    //DISPLAY
    unsigned int start = 0;
    unsigned int stop = min(60, sizeseq);
    outf << "\\texttt{";
    while (start < sizeseq) {
        for (unsigned int spe = 0; spe < 1; spe++) {
            if (seq.species[spe]) {
                for (unsigned int i = start; i < stop; i++) {
                    if (vvstate[spe][i] == 0) {
                        outf << seq.seqsrealigned[spe][i];
                    } else if (vvstate[spe][i] == 1) {
                        outf << "\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                        outf << seq.seqsrealigned[spe][i];
                        outf << "}";
                    } else if (vvstate[spe][i] == 2) {
                        outf << "\\textit{\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                        outf << seq.seqsrealigned[spe][i];
                        outf << "}}";
                    }
                    if (i > 0 && i % 10 == 0) {
                        outf << "\t";
                    }
                }
                outf << "\\\\\n";
                int pstate = 0;
                int ppos = start;
                for (unsigned int i = start; i < stop; i++) {
                    if (vvstate[spe][i] > 0  && vvstate[spe][i] != pstate) {
                        int deca;
                        deca = i - ppos;
                        outf << "\\hspace*{" << deca << "\\charwidth}";
                        if (vvstate[spe][i] == 2) outf << "\\textit{";
                        outf << "\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                        outf << setprecision(1) << fixed <<  scores[spe][i] ;
                        outf << "}";
                        if (vvstate[spe][i] == 2) outf << "}";
                        // to count the number of digits:
                        stringstream score;
                        score << setprecision(1) << fixed <<  scores[spe][i] ;
                        for (unsigned int ii = 1; ii <= score.str().size(); ii++) {
                            if (i > 0 && i % 10 == 0) {
                                outf << "\\hspace*{1\\charwidth}";
                            }
                            i++;
                        }
                        ppos = i;
                    }
                    pstate = vvstate[spe][i];
                    if (i > 0 && i % 10 == 0) {
                        outf << "\\hspace*{1\\charwidth}";
                    }
                }
                outf << "\\\\\n";
            }
        }
        start = stop;
        stop = min(stop + 60, sizeseq);
    }
    outf << "}\n";
    outf << "\\\\\n";
}

void
dispmotifs_html(ofstream & outf, vmot & mots)
{
    outf << "<h2> Motifs</h2>" << endl;
    outf << "<table>" << endl;
    outf << "<tr>" << endl;
    outf << "<th>Color</th><th>Rank</th> <th>Logo</th><th>P-value</th><th>Over-representation</th>" << endl;
    outf << "</tr>" << endl;
    unsigned int i = 1;
    for (ivmot iv = mots.begin(); iv != mots.end(); iv++) {
        outf << "<tr>" << endl;
        outf << "<td class=\"mot" << i << "\"></td><td>" << i << "</td> <td><a href=\"Motif" << i;
        outf << ".pdf\"><img src=\"Motif" << i << ".png\" alt=\"Motif " << i;
        outf << "\" /></a></td><td>" << (*iv).pvalue << "</td><td>" << (*iv).lambdatrain / (*iv).lambdaback << "</td>" << endl;
        outf << "</tr>" << endl;
        i++;
    }
    outf << "</table>" << endl;
}

// display TFBS on reference sequences
void
dispseqwmots_html(Sequence & seq, vmot & mots, ofstream & outf)
{
    //HEADER
    outf << "<h3>";
    string name = numtospecies(0);
    outf << ">" << name <<  " ";
    outf << seq.name << " " ;
    outf << chromfromint(seq.chrom) << " ";
    outf << seq.start << " " ;
    outf << seq.stop << "</h3>" << endl;
    //DEFINING STATES, COLORS AND SCORES
    unsigned int sizeseq = seq.seqsrealigned[0].size();
    vint vdum(sizeseq, 0);
    vvint vvstate(nbspecies, vdum); // 0 nothing, 1 normal tfbs, 2 cons tfbs
    vvint vvcol(nbspecies, vdum); // 0 red, 1 blue, 2 green, 3 yellow
    vd vs(sizeseq, 0.);
    vvd scores(nbspecies, vs);
    //NOT CONS
    for (ivinstseq ivs = seq.instances.begin(); ivs != seq.instances.end(); ivs++) {
        int spe = ivs->species;
        int motwidth = mots[ivs->motindex].motwidth;
        int ipos = seq.imaps[spe][ivs->pos];
        int stop = seq.imapsinv[spe][ipos + motwidth - 1];
        for (int i = ivs->pos; i < stop + 1; i++) {
            vvstate[spe][i] = 1;
            vvcol[spe][i] = ivs->motindex;
            scores[spe][i] = ivs->score;
            //if (seq.seqsrealigned[spe][i]=='-') vvcol[spe][i]=-1;
        }
    }
    //CONS
    for (ivvinstseq ivvs = seq.instancescons.begin(); ivvs != seq.instancescons.end(); ivvs++) {
        for (ivinstseq ivs = (*ivvs).begin(); ivs != (*ivvs).end(); ivs++) {
            int spe = ivs->species;
            int motwidth = mots[ivs->motindex].motwidth;
            int ipos = seq.imaps[spe][ivs->pos];
            int stop = seq.imapsinv[spe][ipos + motwidth - 1];
            for (int i = ivs->pos; i < stop + 1; i++) {
                vvstate[spe][i] = 2;
                scores[spe][i] = ivs->score;
            }
        }
    }
    //DISPLAY
    unsigned int start = 0;
    unsigned int stop2 = min(60, sizeseq);
    outf << "<pre>" << endl;
    int modulo = 0;
    int check = 0;
    while (start < sizeseq) {
        for (unsigned int i = start; i < stop2; i++) {
            if (i >= sizeseq) break;
            if (vvstate[0][i] == 0) {
                if (check == 1) {
                    outf << "</span>";
                    check = 0;
                }
            } else if (vvstate[0][i] == 1) {
                if (check == 0) {
                    check = 1;
                    outf << "<span class=\"mot" << vvcol[0][i] + 1 << "\">";
                }
            } else if (vvstate[0][i] == 2) {
                if (check == 0) {
                    check = 1;
                    outf << "<span class=\"conserved mot" << vvcol[0][i] + 1 << "\">";
                }
            }
            if (seq.seqsrealigned[0][i] != '-') {
                outf << seq.seqsrealigned[0][i];
                if (i > 0 && (i + 1 - modulo) % 10 == 0) {
                    outf << "  ";
                }
            } else {
                stop2++;
                modulo++;
            }
        }
        int pstate = 0;
        int ppos = start;
        //      for (unsigned int i=start;i<stop;i++){
        //         if (vvstate[0][i]>0  && vvstate[0][i]!=pstate){
        //            int deca;
        //            deca=i-ppos;
        //            outf << deca ;
        //            if (vvstate[0][i]==2) outf << "<span class=emph>";
        //            outf << "mot" << vvcol[0][i];
        //            outf << setprecision(1) << fixed <<  scores[0][i] ;
        //            if (vvstate[0][i]==2) outf << "}";
        //            // to count the number of digits:
        //            stringstream score;
        //            score << setprecision(1) << fixed <<  scores[0][i] ;
        //            for (unsigned int ii=1;ii<=score.str().size();ii++){
        //               if ( i>0 && i%10 == 0 ) {
        //                  outf << " ";
        //               }
        //               i++;
        //            }
        //            ppos=i;
        //         }
        //         pstate=vvstate[0][i];
        //         if ( i>0 && (i+1)%10 == 0 ) {
        //            outf << " ";
        //         }
        //      }
        start = stop2;
        stop2 = min(stop2 + 60, sizeseq);
        outf << endl;
    }
    outf << "</pre>" << endl;
}

// display TFBS on aligned sequences
void
dispseqwmotswgaps(Sequence & seq, vmot & mots, ofstream & outf)
{
    //HEADER
    string name = numtospecies(0);
    string texname = seq.name;
    outf << "$>$" << name <<  "\t";
    outf << texname << "\t" ;
    outf << "chr" << chromfromint(seq.chrom) << "\t";
    outf << seq.start << "\t" ;
    outf << seq.stop << " \\\\\n";
    outf << " \\\\\n";
    //DEFINING STATES, COLORS AND SCORES
    unsigned int sizeseq = seq.seqsrealigned[0].size();
    vint vdum(sizeseq, 0);
    vvint vvstate(nbspecies, vdum); // 0 nothing, 1 normal tfbs, 2 cons tfbs
    vvint vvcol(nbspecies, vdum); // 0 red, 1 blue, 2 green, 3 yellow
    vd vs(sizeseq, 0.);
    vvd scores(nbspecies, vs);
    //NOT CONS
    for (ivinstseq ivs = seq.instances.begin(); ivs != seq.instances.end(); ivs++) {
        int spe = ivs->species;
        int motwidth = mots[ivs->motindex].motwidth;
        int ipos = seq.imaps[spe][ivs->pos];
        int stop = seq.imapsinv[spe][ipos + motwidth - 1];
        for (int i = ivs->pos; i < stop + 1; i++) {
            vvstate[spe][i] = 1;
            vvcol[spe][i] = ivs->motindex;
            scores[spe][i] = ivs->score;
            //if (seq.seqsrealigned[spe][i]=='-') vvcol[spe][i]=-1;
        }
    }
    //CONS
    for (ivvinstseq ivvs = seq.instancescons.begin(); ivvs != seq.instancescons.end(); ivvs++) {
        for (ivinstseq ivs = (*ivvs).begin(); ivs != (*ivvs).end(); ivs++) {
            int spe = ivs->species;
            int motwidth = mots[ivs->motindex].motwidth;
            int ipos = seq.imaps[spe][ivs->pos];
            int stop = seq.imapsinv[spe][ipos + motwidth - 1];
            for (int i = ivs->pos; i < stop + 1; i++) {
                vvstate[spe][i] = 2;
                scores[spe][i] = ivs->score;
            }
        }
    }
    //DISPLAY
    unsigned int start = 0;
    unsigned int stop = min(60, sizeseq);
    outf << "\\texttt{";
    while (start < sizeseq) {
        for (unsigned int spe = 0; spe < nbspecies; spe++) {
            if (seq.species[spe]) {
                outf << numtospecies(spe) << "\t";
                for (unsigned int i = start; i < stop; i++) {
                    if (vvstate[spe][i] == 0) {
                        outf << seq.seqsrealigned[spe][i];
                    } else if (vvstate[spe][i] == 1) {
                        outf << "\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                        outf << seq.seqsrealigned[spe][i];
                        outf << "}";
                    } else if (vvstate[spe][i] == 2) {
                        outf << "\\textit{\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                        outf << seq.seqsrealigned[spe][i];
                        outf << "}}";
                    }
                    if (i > 0 && i % 10 == 0) {
                        outf << "\t";
                    }
                }
                outf << "\\\\\n";
                int pstate = 0;
                int ppos = start;
                outf << "\\hspace*{" << 7 << "\\charwidth}";// species name + space = 7 characters
                for (unsigned int i = start; i < stop; i++) {
                    if (vvstate[spe][i] > 0  && vvstate[spe][i] != pstate) {
                        int deca;
                        deca = i - ppos;
                        outf << "\\hspace*{" << deca << "\\charwidth}";
                        if (vvstate[spe][i] == 2) outf << "\\textit{";
                        outf << "\\textcolor{" << colfromint(vvcol[spe][i]) << "}{";
                        outf << setprecision(1) << fixed <<  scores[spe][i] ;
                        outf << "}";
                        if (vvstate[spe][i] == 2) outf << "}";
                        // to count the number of digits:
                        stringstream score;
                        score << setprecision(1) << fixed <<  scores[spe][i] ;
                        for (unsigned int ii = 1; ii <= score.str().size(); ii++) {
                            if (i > 0 && i % 10 == 0) {
                                outf << "\\hspace*{1\\charwidth}";
                            }
                            i++;
                        }
                        ppos = i;
                    }
                    pstate = vvstate[spe][i];
                    if (i > 0 && i % 10 == 0) {
                        outf << "\\hspace*{1\\charwidth}";
                    }
                }
                outf << "\\\\\n";
            }
        }
        outf << "\\\\\n";
        start = stop;
        stop = min(stop + 60, sizeseq);
    }
    outf << "}\n";
}

// general tex header
// TODO: more portable for international languages?
void
disptexinit(ofstream & outf)
{
    outf << "\\documentclass[11pt,twoside,reqno,a4paper]{article}\n" <<
         "\\usepackage[french]{babel}\n" <<
         "\\usepackage{color}\n" <<
         "\\usepackage[utf8]{inputenc}\n" <<
         "\\usepackage{geometry}\n" <<
         "\\geometry{a4paper}\n" <<
         "\\usepackage{graphicx}\n" <<
         "\\usepackage{array,multirow}\n" <<
         "\\begin{document}\n" <<
         "\\noindent\n" <<
         "\\newlength{\\charwidth}" <<
         "\\settowidth{\\charwidth}{\\texttt{A}}";
}

// general html header
// TODO: more portable for international languages?
void
disphtmlinit(ofstream & outf, string title)
{
    outf << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\" ";
    outf << "\"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\">" << endl;
    outf << endl;
    outf << "<html xmlns=\"http://www.w3.org/1999/xhtml\">" << endl;
    outf << endl;
    outf << "<head>" << endl;
    outf << "<link rel=\"stylesheet\" href=\"css/screen.css\"";
    outf << " type=\"text/css\" media=\"screen, projection\" />" << endl;
    outf << "<link rel=\"stylesheet\" href=\"css/print.css\"";
    outf << " type=\"text/css\" media=\"print\" />" << endl;
    outf << "<!--[if lt IE 8]>" << endl;
    outf << "<link rel=\"stylesheet\" href=\"css/ie.css\"";
    outf << " type=\"text/css\" media=\"screen, projection\" />" << endl;
    outf << "<![endif]-->" << endl;
    outf << "<link rel=\"stylesheet\" href=\"css/custom.css\"";
    outf << " type=\"text/css\" media=\"print, projection, screen\" />" << endl;
    outf << endl;
    outf << "<meta http-equiv=\"Content-Type\" content=\"text/html;charset=utf-8\" />" << endl;
    outf << "<title>" << title << "</title>" << endl;
    outf << "</head>" << endl;
    outf << endl;
    outf << "<body>" << endl;
    outf << "<div class=\"container\">" << endl;
    outf << "<div class=\"span-24\">" << endl;
    outf << "<h1>" << title << "</h1>" << endl;
}

void
disphtmlclose(ofstream & outf)
{
    outf << "</div>" << endl;
    outf << "</div>" << endl;
    outf << "</body>" << endl;
    outf << endl;
    outf << "</html>" << endl;
}

void
disptexclose(ofstream & outf)
{
    outf << "\\end{document}";
}

void
disptex(vseq & seqs, vmot & mots)
{
    string filename;
    filename = "results.tex";
    ofstream outf(filename.c_str());
    if (outf.fail()) {
        cerr << "Cannot open file for tex recording: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    disptexinit(outf);
    for (ivseq ivs = seqs.begin(); ivs != seqs.end(); ivs++) {
        cout << "Scanning " << ivs->name << endl;
        dispseqwmots(*ivs, mots, outf);
    }
    disptexclose(outf);
    outf.close();
}

void
copycss(string cssname)
{
    ifstream fin((display_datapath + "/css/" + cssname).c_str());
    ofstream fout(("css/" + cssname).c_str(), fstream::trunc);
    if (fin.fail() || fout.fail()) {
        cerr << "Cannot open " << cssname << ": " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    fout << fin.rdbuf();
}

void
disphtml_genmot(vseq & seqs, vmot & mots)
{
    // Generate svg files
    for (ivseq ivs = seqs.begin(); ivs != seqs.end(); ivs++) {
        scanseqforsvg(*ivs, mots);
    }
    string filename;
    filename = "results_genmot.html";
    ofstream outf(filename.c_str());
    if (outf.fail()) {
        cerr << "Cannot open file for html recording: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    int retcode = mkdir("css" , S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    if (retcode) {
        cerr << "Cannot create directory: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    int retsym = symlink((display_datapath + "/css").c_str(), "css");
    if (retsym && errno != EEXIST) {
        cerr << "Cannot create symlink: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    copycss("screen.css");
    copycss("print.css");
    copycss("ie.css");
    copycss("custom.css");
    disphtmlinit(outf, "Imogene genmot output");
    dispmotifs_html(outf, mots);
    outf << "<h2>Motif instances in the training set</h2>";
    for (ivseq ivs = seqs.begin(); ivs != seqs.end(); ivs++) {
        cout << "Scanning " << ivs->name << endl;
        dispseqwmots_html(*ivs, mots, outf);
    }
    outf << "<h2>Motif presence in alignments</h2>";
    for (ivseq ivs = seqs.begin(); ivs != seqs.end(); ivs++) {
        outf << "<h3>" << (*ivs).name << "</h3>" << endl;
        string filename;
        filename += (*ivs).name;
        filename += ".svg";
        outf << "<object data=\"" + filename + "\" width=\"950px\" type=\"image/svg+xml\"></object>" << endl;
    }
    disphtmlclose(outf);
    outf.close();
}

void
disphtml_scangen()
{
    string filename;
    filename = "results_scangen.html";
    ofstream outf(filename.c_str());
    if (outf.fail()) {
        cerr << "Cannot open file for html recording: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    int retsym = symlink((display_datapath + "/css").c_str(), "css");
    if (retsym && errno != EEXIST) {
        cerr << "Cannot create symlink: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    disphtmlinit(outf, "Imogene scangen output");
    outf << "<h2>Genome-wide predicted enhancers</h2>" << endl;
    outf << "<p>TSS stands for Transcription start sites</p>" << endl;
    ifstream enhf(display_args.enhancers_arg);
    if (enhf.fail()) {
        cerr << "Cannot open enhancer file: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    outf << "<table>" << endl;
    outf << "<tr>" << endl;
    outf << "<th>Score</th><th>Coordinate</th> <th>Closest TSS</th><th>distance to closest TSS</th><th>TSSs</th>" << endl;
    outf << "</tr>" << endl;
    unsigned int i = 0;
    string strline;
    string strfield;
    getline(enhf, strline);
    while (enhf) {
        stringstream line(strline);
        getline(line, strfield, ' ');
        stringstream field(strfield);
        double score;
        field >> score;
        string coordinate;
        getline(line, coordinate, ' ');
        string closesttss;
        getline(line, closesttss, ' ');
        string disttss;
        getline(line, disttss, ' ');
        int lendisttss = disttss.length();
        string tsses;
        getline(line, tsses, ' ');
        outf << "<tr><td>" << score << "</td><td>" << coordinate << "</td><td>";
        outf << closesttss << "</td><td>" << disttss.substr(1, lendisttss - 2) << "</td><td>" << tsses << "</td></tr>" << endl;
        if (i > 200) {
            break;
        }
        i++;
        getline(enhf, strline);
    }
    if (i == 0) { //No motifs
        cerr << "No motifs" << endl;
        exit(EXIT_FAILURE);
    }
    disphtmlclose(outf);
    outf.close();
}

void
disptexwgaps(vseq & align, vmot & mots)
{
    string folder("");
    for (ivseq ivs = align.begin(); ivs != align.end(); ivs++) {
        stringstream file;
        file << folder;
        file << ivs->name;
        file << ".tex";
        ofstream outf(file.str().c_str());
        disptexinit(outf);
        cout << "Scanning " << (*ivs).name << endl;
        dispseqwmotswgaps(*ivs, mots, outf);
        disptexclose(outf);
        outf.close();
    }
}

void
display_args_init()
{
    if (display_args.genmot_mode_counter) {
        if (!strcmp(display_args.species_arg, "droso")) {
            species = "droso";
            nbspecies = 12;
            conca = 0.3;
        } else if (!strcmp(display_args.species_arg, "eutherian")) {
            species = "eutherian";
            nbspecies = 12;
            conca = 0.263;
        }
        concc = 0.5 - conca;
        conct = conca;
        concg = concc;
        scorethr2 = display_args.threshold_arg * log(2);
        //scorethr=scorethr2-1;
        //scorethrcons=scorethr2-1;
        nbmots_for_score = display_args.nbmots_arg;
    }
}

void
loadenhancers(const char * filename)
{
    if (!filename) {
        cerr << "Please give a enhancer file. Exiting..." << endl;
        exit(EXIT_FAILURE);
    }
    ifstream fmotifs;
    fmotifs.open(filename);
    if (fmotifs.fail()) {
        cerr << "Cannot open enhancer file: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    string dum;
    fmotifs >> dum;
}

string display_datapath;

/**
 * ===  FUNCTION  ======================================================================
 *         Name:  cmd_display
 *  Description:  Display results
 * =====================================================================================
 */
int
cmd_display(int argc, char ** argv)
{
    if (display_cmdline_parser(argc, argv, & display_args) != 0)
        exit(EXIT_FAILURE);
    display_args_init();
    const char * imo_display_datapath = getenv("IMOGENE_DATA");
    if (imo_display_datapath == NULL) {
        display_datapath = DATA_PATH;
    } else {
        display_datapath = imo_display_datapath;
    }
    if (display_args.genmot_mode_counter) {
        cout << "Thresholds: thr2=" << scorethr2 / log(2) << " thr=" << scorethr / log(2) << " thrcons=" << scorethrcons / log(2) << endl;
        cout << "Loading alignments " << endl;
        vseq align;
        align = loadseqs(display_args.align_arg);
        cout << "Nb sequences to scan: " << align.size() << endl;
        cout << "Loading Motifs" << endl;
        vmot mots;
        loadmots(display_args.motifs_arg, mots);
        if (nbmots_for_score < mots.size()) mots.erase(mots.begin() + nbmots_for_score, mots.end());
        for (ivmot iv = mots.begin(); iv != mots.end(); iv++) {
            // use same score on all species for detection
            iv->motscorethrcons = iv->motscorethr2;
            // avoid problems with _ and # characters for latex
            texify(iv->name);
        }
        cout << "Loaded " << mots.size() << " motifs." << endl;
        cout << "Scanning sequences for instances..." << endl;
        scanseqsforinstancesnmask(align, mots);
        cout << "Defining conserved instances..." << endl;
        for (ivseq ivs = align.begin(); ivs != align.end(); ivs++) {
            ivs->instances2instancescons();
        }
        if (display_args.tex_ref_given) {
            cout << "Creating tex file for reference species... " << endl;
            disptex(align, mots);
        } else if (display_args.html_ref_given) {
            cout << "Creating html file for reference species... " << endl;
            disphtml_genmot(align, mots);
        } else if (display_args.tex_align_given) {
            cout << "Creating fasta/tex files... " << endl;
            disptexwgaps(align, mots);
        } else if (display_args.jaspar_given) {
            cout << "Outputting matrices... " << endl;
            dispjaspar(mots);
        } else if (display_args.logos_given) {
            cout << "Creating logos... " << endl;
            dispweblogo(mots);
        } else if (display_args.svg_given) {
            for (ivseq ivs = align.begin(); ivs != align.end(); ivs++) {
                scanseqforsvg(*ivs, mots);
            }
        } else {
            cout << "No output mode was given. Exiting." << endl;
        }
    } else {
        if (display_args.scangen_mode_counter) {
            disphtml_scangen();
        }
    }
    return EXIT_SUCCESS;
}		/* -----  end of function extract  ----- */
