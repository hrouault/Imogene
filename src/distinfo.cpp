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

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>

#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include	<gsl/gsl_multimin.h>

#include "distinfo_cmdline.h"
#include "const.hpp"
#include "vectortypes.hpp"
#include "motif.hpp"
#include "tree.hpp"

vd fback;

double
scalprod(vd & col1, vd & col2)
{
    ivd iv2 = col2.begin();
    double sum = 0;
    for (ivd iv = col1.begin(); iv != col1.end(); iv++) {
        sum += *iv * *iv2;
        iv2++;
    }
    return sum;
}		/* -----  end of function sumcol  ----- */

double
ftomin(double x, vvd mat)
{
    double sum = 0;
    for (ivvd iv = mat.begin(); iv != mat.end(); iv++) {
        vd cpcol;
        for (ivd ic = (*iv).begin(); ic != (*iv).end(); ic++) {
            cpcol.push_back(exp(-*ic * x));
        }
        sum += log(scalprod(fback, cpcol));
    }
    //cout << x << " " << sum+scorethr*x << endl;
    return sum + scorethr * x;
}

double
fwrap(double x, void * params)
{
    const vvd * mats = (const vvd *)params;
    const vvd & mat = *mats;
    return ftomin(x, mat);
}		/* -----  end of function f2wrap  ----- */

double
ftominp(double x, const vvd & mat)
{
    double sum = 0;
    for (civvd iv = mat.begin(); iv != mat.end(); iv++) {
        vd cpcol, cpcol2;
        for (civd ic = (*iv).begin(); ic != (*iv).end(); ic++) {
            cpcol.push_back(-*ic * exp(-*ic * x));
            cpcol2.push_back(exp(-*ic * x));
        }
        sum += scalprod(fback, cpcol) / scalprod(fback, cpcol2);
    }
    return sum + scorethr;
}
double
ftominpp(double x, const vvd & mat)
{
    double sumed = 0;
    for (civvd iv = mat.begin(); iv != mat.end(); iv++) {
        vd z, emoy, e2moy;
        for (civd ic = (*iv).begin(); ic != (*iv).end(); ic++) {
            z.push_back(exp(-*ic * x));
            emoy.push_back(*ic * exp(-*ic * x));
            e2moy.push_back((*ic) * (*ic)*exp(-*ic * x));
        }
        double sz = scalprod(fback, z);
        double semoy = scalprod(fback, emoy);
        double se2moy = scalprod(fback, e2moy);
        sumed += se2moy / sz - semoy * semoy / (sz * sz);
    }
    return sumed;
}		/* -----  end of function ftominpp  ----- */

double
ftomin2(vd x, const vvd & mat1, const vvd & mat2)
{
    double s = 0;
    civvd iv2 = mat2.begin();
    for (civvd iv = mat1.begin(); iv != mat1.end(); iv++) {
        vd ccol;
        civd ic2 = (*iv2).begin();
        for (civd ic = (*iv).begin(); ic != (*iv).end(); ic++) {
            ccol.push_back(exp(-*ic * x[0] - *ic2 * x[1]));
            ic2++;
        }
        s += log(scalprod(fback, ccol));
        iv2++;
    }
    s += scorethr1 * x[0] + scorethr2 * x[1];
    //   cout << "sum : " << s << endl;
    return s;
}		/* -----  end of function ftomin2  ----- */


double
f2wrap(const gsl_vector * v, void * params)
{
    const vvd ** mats = (const vvd **)params;
    const vvd & mat1 = *mats[0];
    const vvd & mat2 = *mats[1];
    vd x;
    x.push_back(gsl_vector_get(v, 0));
    x.push_back(gsl_vector_get(v, 1));
    return ftomin2(x, mat1, mat2);
}		/* -----  end of function f2wrap  ----- */


vd
ftomin2p(const vd x, const vvd & mat1, const vvd & mat2)
{
    double s1 = 0;
    double s2 = 0;
    civvd iv2 = mat2.begin();
    for (civvd iv = mat1.begin(); iv != mat1.end(); iv++) {
        vd ccol, ccolp;
        vd ccol2;
        civd ic2 = (*iv2).begin();
        for (civd ic = (*iv).begin(); ic != (*iv).end(); ic++) {
            ccol.push_back(-*ic * exp(-*ic * x[0] - *ic2 * x[1]));
            ccolp.push_back(exp(-*ic * x[0] - *ic2 * x[1]));
            ccol2.push_back(-*ic2 * exp(-*ic * x[0] - *ic2 * x[1]));
            ic2++;
        }
        double sp = scalprod(fback, ccolp);
        s1 += scalprod(fback, ccol) / sp;
        s2 += scalprod(fback, ccol2) / sp;
        iv2++;
    }
    vd res;
    res.push_back(s1);
    res.push_back(s2);
    return res;
}		/* -----  end of function ftomin2p  ----- */

vvd
ftomin2pp(const vd & x, const vvd & mat1, const vvd & mat2)
{
    double dxx = 0, dxxp = 0, dxpxp = 0;
    civvd iv2 = mat2.begin();
    for (civvd iv = mat1.begin(); iv != mat1.end(); iv++) {
        vd z, emoy, epmoy, e2moy, ep2moy, epemoy;
        civd ic2 = (*iv2).begin();
        for (civd ic = (*iv).begin(); ic != (*iv).end(); ic++) {
            z.push_back(exp(-*ic * x[0] - *ic2 * x[1]));
            emoy.push_back(-*ic * exp(-*ic * x[0] - *ic2 * x[1]));
            epmoy.push_back(-*ic2 * exp(-*ic * x[0] - *ic2 * x[1]));
            e2moy.push_back(*ic * (*ic)*exp(-*ic * x[0] - *ic2 * x[1]));
            ep2moy.push_back(*ic2 * (*ic2)*exp(-*ic * x[0] - *ic2 * x[1]));
            epemoy.push_back(*ic2 * (*ic)*exp(-*ic * x[0] - *ic2 * x[1]));
            ic2++;
        }
        double sz = scalprod(fback, z);
        double semoy = scalprod(fback, emoy);
        double sepmoy = scalprod(fback, epmoy);
        double se2moy = scalprod(fback, e2moy);
        double sep2moy = scalprod(fback, ep2moy);
        double sepemoy = scalprod(fback, epemoy);
        dxx += se2moy / sz - semoy * semoy / (sz * sz);
        dxxp += sepemoy / sz - (semoy * sepmoy) / (sz * sz);
        dxpxp += sep2moy / sz - (sepmoy * sepmoy) / (sz * sz);
        iv2++;
    }
    //   cout << dxx << " " << dxxp << " " << dxpxp << endl;
    vd res1;
    vd res2;
    res1.push_back(dxx);
    res1.push_back(dxxp);
    res2.push_back(dxxp);
    res2.push_back(dxpxp);
    vvd res;
    res.push_back(res1);
    res.push_back(res2);
    return res;
}		/* -----  end of function f2min2pp  ----- */

void
f2pwrap(const gsl_vector * v, void * params, gsl_vector * df)
{
    const vvd ** mats = (const vvd **)params;
    const vvd & mat1 = *mats[0];
    const vvd & mat2 = *mats[1];
    vd x;
    x.push_back(gsl_vector_get(v, 0));
    x.push_back(gsl_vector_get(v, 1));
    vd gradf = ftomin2p(x, mat1, mat2);
    gsl_vector_set(df, 0, gradf[0]);
    gsl_vector_set(df, 1, gradf[1]);
}		/* -----  end of function f2pwrap  ----- */

void
f2f2pwrap(const gsl_vector * v, void * params, double * f, gsl_vector * df)
{
    const vvd ** mats = (const vvd **)params;
    const vvd & mat1 = *mats[0];
    const vvd & mat2 = *mats[1];
    vd x;
    x.push_back(gsl_vector_get(v, 0));
    x.push_back(gsl_vector_get(v, 1));
    vd gradf = ftomin2p(x, mat1, mat2);
    *f = ftomin2(x, mat1, mat2);
    gsl_vector_set(df, 0, gradf[0]);
    gsl_vector_set(df, 1, gradf[1]);
}		/* -----  end of function f2pwrap  ----- */

vd
solvex2(const vvd & mat1, const vvd & mat2)
{
    const gsl_multimin_fdfminimizer_type * T;
    gsl_multimin_fdfminimizer * s;
    /* Position of the minimum (1,2), scale factors
       10,20, height 30. */
    const vvd * par[2] = { &mat1, &mat2 };
    gsl_vector * x;
    gsl_multimin_function_fdf my_func;
    my_func.n = 2;
    my_func.f = &f2wrap;
    my_func.df = &f2pwrap;
    my_func.fdf = &f2f2pwrap;
    my_func.params = &par;
    /* Starting point, x = (-1,-1) */
    x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, -1.0);
    gsl_vector_set(x, 1, -1.0);
    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc(T, 2);
    gsl_multimin_fdfminimizer_set(s, &my_func, x, 0.01, 1e-3);
    int iter = 0;
    int status = 0;
    do {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);
        double x0 = gsl_vector_get(s->x, 0);
        double x1 = gsl_vector_get(s->x, 0);
        if (x0 != x0 && x1 != x1) {
            vd dum;
            dum.push_back(-1000);
            dum.push_back(-1000);
            return dum;
        }
        if (status)
            break;
        //cout << gsl_vector_get(x,0) << gsl_vector_get(x,1) << s->f << endl;
        status = gsl_multimin_test_gradient(s->gradient, 1e-3);
        //
        //      if (status == GSL_SUCCESS)
        //         printf ("Minimum found at:\n");
        //
        //      printf ("%5d %.5f %.5f %10.5f\n", iter,
        //            gsl_vector_get (s->x, 0),
        //            gsl_vector_get (s->x, 1),
        //            s->f);
    } while (status == GSL_CONTINUE && iter < 1000);
    //cout << "minimum found at : " << gsl_vector_get(s->x, 0) << " " << gsl_vector_get (s->x, 1) << " " << s->f << endl;
    vd res;
    res.push_back(gsl_vector_get(s->x, 0));
    res.push_back(gsl_vector_get(s->x, 1));
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
    return res;
}

double
solvex(const vvd & mat)
{
    int status;
    int iter = 0, max_iter = 1000;
    const gsl_min_fminimizer_type * T;
    gsl_min_fminimizer * s;
    double m = -1.0;
    double a = -100.0, b = 100.0;
    gsl_function F;
    F.function = &fwrap;
    F.params = (void *)&mat;
    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc(T);
    gsl_min_fminimizer_set(s, &F, m, a, b);
    do {
        iter++;
        status = gsl_min_fminimizer_iterate(s);
        m = gsl_min_fminimizer_x_minimum(s);
        a = gsl_min_fminimizer_x_lower(s);
        b = gsl_min_fminimizer_x_upper(s);
        status
            = gsl_min_test_interval(a, b, 0.001, 0.0);
        // if (status == GSL_SUCCESS)
        //printf ("Converged\n");
    } while (status == GSL_CONTINUE && iter < max_iter);
    // cout << "converged to : " << m << endl;
    gsl_min_fminimizer_free(s);
    return m;
}

double
p(double x, const vvd & mat)
{
    double exppart = exp(ftomin(x, mat));
    double prefact = sqrt(2 * M_PI / ftominpp(x, mat));
    return prefact * exppart;
}		/* -----  end of function p  ----- */


double
p2(const vd & x, const vvd & mat1, const vvd & mat2)
{
    double exppart = exp(ftomin2(x, mat1, mat2));
    //   cout << "exppart : " << exppart << endl;
    vvd fp2 = ftomin2pp(x, mat1, mat2);
    double prefact = 2 * M_PI * sqrt(1 / (fp2[0][0] * fp2[1][1] - fp2[0][1] * fp2[1][0]));
    return exppart * prefact;
}		/* -----  end of function p2  ----- */

double
distmatinfo(const vvd & mat1, const vvd & mat2)
{
    double scoreinit = scorethr;
    scorethr = scorethr1;
    double x1 = solvex(mat1);
    scorethr = scorethr2;
    double x2 = solvex(mat2);
    vd x1x2 = solvex2(mat1, mat2);
    scorethr = scoreinit;
    if (x1x2[0] < -800 && x1x2[1] < -800) return 1e5;
    //   cout << x1 << " " << p(x1,mat1) << endl;
    //   cout << x2 << " " << p(x2,mat2) << endl;
    //   cout << x1x2 << endl;
    //   cout << "p1 : " << p(x1,mat1) << endl;
    //   cout << "p2 : " << p(x2,mat2) << endl;
    //   cout << "pboth : " << p2(x1x2,mat1,mat2) << endl;
    //cout << "distance : " << -log(2*p2(x1x2,mat1,mat2)/(p(x1,mat1)+p(x2,mat2))) << endl;
    return -log(2 * p2(x1x2, mat1, mat2) / (p(x1, mat1) + p(x2, mat2)));
}

unsigned int decamax;

double
distmat(const vvd & mat1, const vvd & mat2)
{
    double dist = 1e5;
    for (unsigned int deca = 0; deca < decamax; deca++) {
        vd dum(4, 0.0);
        vvd m1(mat1);
        m1.insert(m1.end(), deca, dum);
        vvd m2(mat2);
        m2.insert(m2.begin(), deca, dum);
        //      cout << m1<< " " << m2 << endl;
        if (sum(abs(m1 - m2)) > 1e-3) {
            double disttemp = distmatinfo(m1, m2);
            if (disttemp < dist) {
                dist = disttemp;
            }
        }
    }
    for (unsigned int deca = 1; deca < decamax; deca++) {
        vd dum(4, 0);
        vvd m1(mat1);
        m1.insert(m1.begin(), deca, dum);
        vvd m2(mat2);
        m2.insert(m2.end(), deca, dum);
        //      cout << m1<< " " << m2 << endl;
        if (sum(abs(m1 - m2)) > 1e-3) {
            double disttemp = distmatinfo(m1, m2);
            if (disttemp < dist) {
                dist = disttemp;
            }
        }
    }
    return dist;
}

double
distmot(const Motif & mot1, const Motif & mot2)
{
    scorethr1 = mot1.motscorethr;
    scorethr2 = mot2.motscorethr;
    decamax = max((double)mot1.motwidth / 2, (double)mot2.motwidth / 2);
    //   cout << "deca: " << decamax << endl;
    double dist = 1e5;
    double dist1 = distmat(mot1.matprec, mot2.matprec);
    double dist2 = distmat(mot1.matprec, mot2.matprecrevcomp);
    if (dist1 > dist2) {
        dist = dist2;
    } else {
        dist = dist1;
    }
    return dist;
}

void
compmotsthr(vmot & mots)
{
    vmot fimots;
    unsigned int sizemax = width;
    for (ivmot ivm = mots.begin(); ivm != mots.end(); ivm++) {
        if (ivm->bsinit.size() > sizemax) sizemax = ivm->bsinit.size();
    }
    width = sizemax;
    for (ivmot ivm = mots.begin(); ivm != mots.end(); ivm++) {
        double maxinfo(0);
        double maxelem(0);
        double meaninfo(0);
        int j(0);
        for (ivvd ivv = ivm->matprec.begin(); ivv != ivm->matprec.end(); ivv++) {
            int i(0);
            double maxcol(-10);
            for (ivd iv = ivv->begin(); iv != ivv->end(); iv++) {
                if (*iv > maxcol) maxcol = *iv;
                if (*iv > maxelem) maxelem = *iv;
                meaninfo += ivm->matfreq[j][i] * (*iv);
                i++;
            }
            j++;
            maxinfo += maxcol;
        }
        //cout << ivm->name << " " << maxinfo << " " << maxelem << " " << meaninfo << endl;
        ivm->motscorethr = min(0.95 * maxinfo, meaninfo);
        ivm->motwidth = ivm->bsinit.size();
        unsigned int btomax, bleft, bright;
        btomax = sizemax - ivm->bsinit.size();
        bleft = btomax / 2;
        bright = btomax - bleft;
        vd dum(4, 0.0);
        vvd m1(ivm->matprec);
        m1.insert(m1.begin(), bleft, dum);
        m1.insert(m1.end(), bright, dum);
        ivm->matprec = m1;
        ivm->matprecrevcomp = reversecomp(m1);
        fimots.push_back(*ivm);
    }
    mots = fimots;
    return;
}

distinfo_args_info distinfo_args;


// sets check flags on/off to filter identical motifs in a p-value sorted motifs list
// // *** BEWARE, hidden is a cutoff on vmot size
void
compmotsdist(vmot & mots)
{
    if (mots.size() == 0) return;
    fback.push_back(conca);
    fback.push_back(concc);
    fback.push_back(concg);
    fback.push_back(conct);
    unsigned int counter = 2;
    unsigned int countertrue = 1;
    ivmot ivstop = mots.end() - 1;
    for (ivmot ivm = mots.begin() + 1; ivm != mots.end(); ivm++) {
        if (progress) {
            cout << "\r" << counter << "/" << mots.size();
            cout.flush();
        }
        counter++;
        for (ivmot ivm2 = mots.begin(); ivm2 != ivm; ivm2++) {
            double distance = 0;
            if ((*ivm2).check) {
                if (sum(abs((*ivm).matprec - (*ivm2).matprec)) > 1e-3 && sum(abs((*ivm).matprec - (*ivm2).matprecrevcomp)) > 1e-3) {
                    distance = distmot(*ivm, *ivm2);
                    if (distinfo_args.displaydist_given){
                       cout << ivm->bsinit << ", " << ivm2->bsinit;
                       cout << " : " << distance << endl;
                    }
                }
                if (distance < 3.0) {
                    (*ivm).check = false;
                    break;
                }
            }
        }
        if (ivm->check == true) countertrue++;
        if (countertrue >= 20) { // *** CUTOFF TO AVOID LONG COMPUTATION
            cout << "\nMore than 20 motifs where found. Stopping clustering." << endl;
            ivstop = ivm;
            break;
        }
    }
    if (ivstop != mots.end() - 1) mots.erase(ivstop + 1, mots.end());
    cout << "\n";
    return;
}


double
compdistance(Motif & mot1, Motif & mot2)
{
    //   cout << "conca: " << conca << endl;
    conct = conca;
    concg = concc;
    vmot mots;
    mots.push_back(mot1);
    mots.push_back(mot2);
    compmotsthr(mots);
    mot1 = mots[0];
    mot2 = mots[1];
    
    fback.clear();
    fback.push_back(conca);
    fback.push_back(concc);
    fback.push_back(concg);
    fback.push_back(conct);
    double distance = 0;
    if (sum(abs(mot1.matprec - mot2.matprec)) > 1e-3 && sum(abs(mot1.matprec - mot2.matprecrevcomp)) > 1e-3) {
       distance = distmot(mot1, mot2);
    }
    else distance = 0;

    return distance;
}

void
distinfo(const char * motfile)
{
    initconc();
    if (distinfo_args.width_given) {
        width = distinfo_args.width_arg;
    }
    scorethr = distinfo_args.threshold_arg;
    vmot mots;
    cout << "Loading motifs..." << endl;
    loadmots(motfile, mots);
    cout << "Computing motifs thresholds..." << endl;
    compmotsthr(mots);
    cout << "Maximum size : " << width << endl;
    compmotsdist(mots);
    ofstream outf("bestuniq-pval.dat");
    for (ivmot ivm = mots.begin(); ivm != mots.end(); ivm++) {
        if (ivm->check) {
            ivm->display(outf);
        }
    }
    outf.close();
}				/* ----------  end of function main  ---------- */

/**
 * ===  FUNCTION  ======================================================================
 *         Name:  cmd_distinfo
 *  Description:  Species argument
 * =====================================================================================
 */
int
cmd_distinfo(int argc, char ** argv)
{
    if (distinfo_cmdline_parser(argc, argv, & distinfo_args) != 0)
        exit(EXIT_FAILURE);
    if (strcmp(distinfo_args.species_arg, "droso")) {
        species = "droso";
    } else if (strcmp(distinfo_args.species_arg, "eutherian")) {
        species = "eutherian";
    }
    distinfo(distinfo_args.motifs_arg);
    return EXIT_SUCCESS;
}		/* -----  end of function extract  ----- */

