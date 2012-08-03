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
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>

using namespace std;

#include "genmot_cmdline.h"
#include "genmot.hpp"
#include "const.hpp"
#include "random.hpp"
#include "motif.hpp"
#include "tree.hpp"
#include "distinfo.hpp"
#include "sequence.hpp"
#include "display.hpp"

genmot_args_info genmot_args;

vstring regtests;
vseq regints;

/**
 * ===  FUNCTION  =============================================================
 *         Name:  motscoreorder
 *  Description:  Compare motif overrepresentation
 * ============================================================================
 */
bool
motscoreorder(Motif mot1, Motif mot2)
{
    return mot1.pvalue < mot2.pvalue;
}		/* -----  end of function motscoreorder  ----- */

bool
motchi2order(Motif mot1, Motif mot2)
{
    return mot1.scorepoiss < mot2.scorepoiss;
}		/* -----  end of function motscoreorder  ----- */

// infinite norm
double
distcv(vvd & mat1, vvd & mat2)
{
    double max(0);
    int col(0);
    for (ivvd iv = mat1.begin(); iv != mat1.end(); iv++) {
        int line(0);
        for (ivd il = (*iv).begin(); il != (*iv).end(); il++) {
            double diff;
            diff = fabs((*il) - mat2[col][line]);
            if (diff > max) max = diff;
            line++;
        }
        col++;
    }
    return max;
}

void
motiftomat(vint & seq, Motif & mot)
{
    vd dum(4, 0);
    mot.matrice = vvd(width, dum);
    vvd::iterator imat = mot.matrice.begin();
    for (vint::iterator iseq = seq.begin(); iseq != seq.end(); iseq++) {
        int base = *iseq;
        double a, t, c, g;
        a = t = alpha / (1 + 2 * alpha + 2 * beta);
        c = g = beta / (1 + 2 * alpha + 2 * beta);
        if (base == 0) {
            a = (1 + alpha) / (1 + 2 * alpha + 2 * beta);
        } else if (base == 1) {
            c = (1 + beta) / (1 + 2 * alpha + 2 * beta);
        } else if (base == 2) {
            g = (1 + beta) / (1 + 2 * alpha + 2 * beta);
        } else if (base == 3) {
            t = (1 + alpha) / (1 + 2 * alpha + 2 * beta);
        }
        (*imat)[0] = log(a / conca);
        (*imat)[1] = log(c / concc);
        (*imat)[2] = log(g / concg);
        (*imat)[3] = log(t / conct);
        imat++;
    }
}

void
seqanalysis(Sequence & currseq, vmot & genmots)
{
    unsigned int i = 0;
    for (vint::iterator istr = currseq.iseqs[0].begin();
            istr != currseq.iseqs[0].end() - width + 1; istr++) {
        if (progress) {
            cout << "\r" << i + 1 << "/";
            cout << currseq.iseqs[0].size() - width + 1 ;
            cout.flush();
        }
        vint bs(istr, istr + width);
        if (compN(bs) > 0) continue;
        Motif currmot;
        vma pseqs;
        currmot.bsinit = vinttostring(bs);
        currmot.pos = i;
        motiftomat(bs, currmot);
        currmot.matricerevcomp = reversecomp(currmot.matrice);
        currmot.matprec = currmot.matrice;
        currmot.matprecrevcomp = currmot.matricerevcomp;
        pseqs = currmot.seqs;
        for (unsigned int nb = 1; nb <= nbiter; nb++) {
            if (nb > 2) currmot.matinit(scorethr2);
            else currmot.matinit(scorethr);
            if (currmot.nbcons < 1) break;
            
            if (!strcmp(genmot_args.method_arg, "max")) {
               currmot.compprec();
            } else if (!strcmp(genmot_args.method_arg, "mean")) {
               currmot.compprec_MCMC();
            } else if (!strcmp(genmot_args.method_arg, "inde")) {
               currmot.compprec_inde();
            }

            if (pseqs == currmot.seqs) break;
            pseqs = currmot.seqs;
            if (nb == 1) {
                currmot.corrprec();
                currmot.matprecrevcomp = reversecomp(currmot.matprec);
            }
        }
        currmot.matinit(scorethr2);
        if (currmot.nbcons > 2) {
            currmot.matprecrevcomp = reversecomp(currmot.matprec);
            currmot.matfreq = mattofreq(currmot.matprec);
            currmot.motscorethr = scorethr2;
            currmot.motwidth = width;
            genmots.push_back(currmot);
        }
        i++;
    }
    cout << endl;
}

void
genmot_args_init()
{
    if (!strcmp(genmot_args.species_arg, "droso")) {
        species = "droso";
        nbspecies = 12;
    } else if (!strcmp(genmot_args.species_arg, "eutherian")) {
        species = "eutherian";
        nbspecies = 12;
    }
    initconc();
    width = genmot_args.width_arg;
    scorethr2 = genmot_args.threshold_arg * log(2);
    scorethr = scorethr2 * (1 - 2.0 / width);
    scorethrcons = scorethr2 * (1 - 1.0 / width);
    evolutionary_model = genmot_args.evolutionary_model_arg;
    neighbext = genmot_args.neighbext_arg;
    if (genmot_args.progress_given) progress = true;
}

string genmot_datapath;

/**
 * ===  FUNCTION  =============================================================
 *         Name:  cmd_genmot
 *  Description:  Motif generation
 * ============================================================================
 */
int
cmd_genmot(int argc, char ** argv)
{
    if (genmot_cmdline_parser(argc, argv, & genmot_args) != 0)
        exit(EXIT_FAILURE);
    genmot_args_init();
    rnginit();
    const char * imo_genmot_datapath = getenv("IMOGENE_DATA");
    if (imo_genmot_datapath == NULL) {
        genmot_datapath = DATA_PATH;
    } else {
        genmot_datapath = imo_genmot_datapath;
    }
    sequence_datapath = genmot_datapath;
    // printconfig(); *** to be written so that one can rerun exacltly the same instance
    compalpha();
    //   printpriorsandthrs(); *** to be written
    cout << "alpha=" << alpha << endl;
    cout << "Loading background file names..." << endl;
    if (genmot_args.background_given){
        regtests = loadfilenames(genmot_args.background_arg);
    } else {
        regtests = loadfilenames((genmot_datapath + "/" + species + "/background").c_str());
    }
    cout << "Background set size : " << regtests.size() << endl;
    cout << "Loading training set..." << endl;
    regints = loadseqs(genmot_args.align_arg);
    cout << "Training set size : " << regints.size() << endl;
    cout << "Masking repeats in training set..." << endl;
    maskrepeats(regints);
    inittreedist();
    cout << "Generating motifs..." << endl;
    vmot genmots;
    for (vseq::iterator iseq = regints.begin(); iseq != regints.end(); iseq++) {
        cout << (*iseq).name << endl;
        seqanalysis(*iseq, genmots);
    }
    cout << "Loading background sequences and scoring motifs..." << endl;
    unsigned int counter = 1;
    for (ivstring ivs = regtests.begin(); ivs != regtests.end(); ivs++) {
        if (progress) {
            cout << "\r" << counter << "/" << regtests.size();
            cout.flush();
        }
        counter++;
        Sequence seq;
        string filename = *ivs;
        seq = loadseqconserv(filename);
        for (ivmot ivm = genmots.begin(); ivm != genmots.end(); ivm++) {
            ivm->updatebacksites(seq);
        }
    }
    cout << "\n";
    for (ivmot ivm = genmots.begin(); ivm != genmots.end(); ivm++) {
        ivm->pvaluecomp();
    }
    cout << "Sorting motifs..." << endl;
    // Sort on chi2
    sort(genmots.begin(), genmots.end(), motchi2order);
    // find 3rd quartile
    int shiftchi2 = genmots.size() * 3 / 4;
    genmots.erase(genmots.begin() + shiftchi2, genmots.end());
    // Sort
    sort(genmots.begin(), genmots.end(), motscoreorder);
    // Cluster...
    cout << "Clustering motifs (keeping only 20 best)..." << endl;
    compmotsdist(genmots);
    //
    cout << "Creating output file..." << endl;
    ofstream motmeldb("motifs.txt");
    if (motmeldb.fail()) {
        cerr << "Cannot write to motifs.txt file: " << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }
    for (ivmot ivm = genmots.begin(); ivm != genmots.end(); ivm++) {
        if (ivm->check)
            ivm->display(motmeldb);
    }
    motmeldb.close();
    gsl_rng_free(gslran);
    genmot_cmdline_parser_free(&genmot_args);
    return EXIT_SUCCESS;
}		/* -----  end of function genmot  ----- */


