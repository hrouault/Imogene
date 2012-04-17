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
#include	<cmath>
#include	<gsl/gsl_vector.h>
#include	<gsl/gsl_matrix.h>
#include	<gsl/gsl_blas.h>
#include	<gsl/gsl_odeiv.h>
#include <gsl/gsl_linalg.h> // for exponential of matrices

#include "tree.hpp"
#include "const.hpp"

using namespace std;

vnoe treedist;

double pa;
double pc;

const double kappa = 2.0;

const double integr_step = 0.01; //0.01;//=0.001; //!! 0.001 instead of 0.01
double fat, fac, fag, fta, ftc, ftg, fca, fct, fcg, fga, fgt, fgc;

noeud::noeud(int e1, int e2, int n, double p1, double p2)
{
    esp1 = e1;
    esp2 = e2;
    noe = n;
    if (evolutionary_model == 2) {
        prox1 = p1;
        prox2 = p2;
    } else if (evolutionary_model == 1) {
        double correction = 0.5 + 4.0 * conca * concc;
        //double correction=1.0;
        prox1 = exp(-p1 / correction);
        prox2 = exp(-p2 / correction);
    }
}

vpgslmat vtransi;
gsl_matrix * instrates;
gsl_vector * pnoe1out;
gsl_vector * pnoe2out;
gsl_matrix * id;
gsl_matrix * m1;
gsl_matrix * m2;
gsl_matrix * m3;
gsl_matrix * m4;
gsl_matrix * pij;
gsl_matrix * pijp;

gsl_vector * proba2;
unsigned int noemax;

// ** TODO the following functions use the not-enough-tested-yet
// exponential function from gsl, they could nicely replace
// the current way to compute on the tree.


// This function returns M in
// P(t)=M(t)*P(0)
// where M is exp(R*t), with R the rate matrix
vvd
transition_halpern(vd & w, double dist)
{
    gsl_matrix * rates = gsl_matrix_calloc(4, 4);
    pa = conca;
    pc = concc;
    double w0 = w[0];
    double w1 = w[1];
    double w2 = w[2];
    double w3 = w[3];
    if (w0 < 0 || w1 < 0 || w2 < 0 || w3 < 0) {
        cerr << "Problem in instant_rates_halpern" << endl;
        exit(EXIT_FAILURE);
    }

    double fat = proba_fixation_rel(w3 / w0);
    double fac = proba_fixation_rel(pa * w1 / pc / w0);
    double fag = proba_fixation_rel(pa * w2 / pc / w0);
    
    double fta = proba_fixation_rel(w0 / w3);
    double ftc = proba_fixation_rel(pa * w1 / pc / w3);
    double ftg = proba_fixation_rel(pa * w2 / pc / w3);
    
    double fca = proba_fixation_rel(pc * w0 / pa / w1);
    double fct = proba_fixation_rel(pc * w3 / pa / w1);
    double fcg = proba_fixation_rel(w2 / w1);
    
    double fga = proba_fixation_rel(pc * w0 / pa / w2);
    double fgt = proba_fixation_rel(pc * w3 / pa / w2);
    double fgc = proba_fixation_rel(w1 / w2);
    
    double prefact = 1.0 / (4 * kappa * pa * pc + 0.5);
    
    gsl_matrix_set(m1, 0, 0, -(pa * fat + pc * fac + pc * kappa * fag));
    gsl_matrix_set(m1, 0, 1, pa * fca);
    gsl_matrix_set(m1, 0, 2, pa * kappa * fga);
    gsl_matrix_set(m1, 0, 3, pa * fta);
    
    gsl_matrix_set(m1, 1, 0, pc * fac);
    gsl_matrix_set(m1, 1, 1, -(pa * fca + pa * kappa * fct + pc * fcg));
    gsl_matrix_set(m1, 1, 2, pc * fgc);
    gsl_matrix_set(m1, 1, 3, pc * kappa * ftc);
    
    gsl_matrix_set(m1, 2, 0, pc * kappa * fag);
    gsl_matrix_set(m1, 2, 1, pc * fcg);
    gsl_matrix_set(m1, 2, 2, -(pa * kappa * fga + pa * fgt + pc * fgc));
    gsl_matrix_set(m1, 2, 3, pc * ftg);
    
    gsl_matrix_set(m1, 3, 0, pa * fat);
    gsl_matrix_set(m1, 3, 1, pa * kappa * fct);
    gsl_matrix_set(m1, 3, 2, pa * fgt);
    gsl_matrix_set(m1, 3, 3, -(pa * fta + pc * kappa * ftc + pc * ftg));

    gsl_matrix_scale(m1, prefact * dist);
    // *** This is the exponentiation step
    gsl_linalg_exponential_ss(m1, rates, 1e-3);
    //Matrice d'évolution Runge Kutta 4
    // Mat(RG4) = Id + h*M + h^2/2*M^2 + h^3/6*M^3 + h^4/24*M^4
    //   gsl_matrix * mattemp;
    //
    //   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,0.5,m1,m1,0.0,m2);
    //   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1/3.0,m1,m2,0.0,m3);
    //   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1/4.0,m1,m3,0.0,m4);
    // rates=exp(m1), with precision of 1e-3
    //gsl_matrix_memcpy(rates,id);
    //gsl_matrix_add(rates,m1);
    //gsl_matrix_add(rates,m2);
    //gsl_matrix_add(rates,m3);
    //gsl_matrix_add(rates,m4);
    vd dum(4, 0.);
    vvd r(4, dum);
    for (unsigned int col = 0; col < 4; col++) {
        for (unsigned int row = 0; row < 4; row++) {
            r[col][row] = gsl_matrix_get(rates, row, col);
        }
    }
    return r;
}

vvd
transition_felsen(vd & w, double dist)
{
    vd dum(4, 0.);
    vvd M(4, dum);
    double correction = 0.5 + 4.0 * conca * concc;
    double prox = exp(-dist / correction);
    for (unsigned int col = 0; col < 4; col++) {
        for (unsigned int row = 0; row < 4; row++) {
            if (row == col) M[col][row] = prox + (1 - prox) * w[row];
            else M[col][row] = (1. - prox) * w[row];
        }
    }
    return M;
}
//probs is the initial base prob in the tree, freq its selection frequency, dist

vd
evolvedist_felsen_backwards(vd & probs, vd & freqs, double dist)
{
    vd pf(4, 0.);
    vvd M = transition_felsen(freqs, dist);
    double sum = 0;
    for (unsigned int row = 0; row < 4; row++) {
        for (unsigned int k = 0; k < 4; k++) {
            pf[row] += M[row][k] * probs[k];
        }
        sum += pf[row];
    }
    for (unsigned int row = 0; row < 4; row++) {
        pf[row] /= sum;
    }
    return pf;
}

vd
evolvedist_halpern(vd & probs, vd & freqs, double dist)
{
    vd pf(4, 0.);
    vvd M = transition_halpern(freqs, dist);

    for (unsigned int row = 0; row < 4; row++) {
        for (unsigned int k = 0; k < 4; k++) {
            pf[row] += M[k][row] * probs[k];
        }
    }
    return pf;
}

vd
evolvedist_felsen(vd & probs, vd & freqs, double dist)
{
    vd pf(4, 0.);
    vvd M = transition_felsen(freqs, dist);
    for (unsigned int row = 0; row < 4; row++) {
        for (unsigned int k = 0; k < 4; k++) {
            pf[row] += M[k][row] * probs[k];
        }
    }
    return pf;
}

vd
evolvedist(vd probs, vd freqs, double dist)
{
    vd dum;
    if (evolutionary_model == 1) return evolvedist_felsen(probs, freqs, dist);
    else if (evolutionary_model == 2) return evolvedist_halpern(probs, freqs, dist);
    return dum;
}

void   // Currently used phylogenetic tree for drosophilae :  Heger and Pontig, 2007
inittreedist()
{
    pa = conca;
    pc = concc;
    treedist.clear();
    if (species == "droso") {
        treedist.push_back(noeud(1, 2, 12, 0.017, 0.021)); // 12
        treedist.push_back(noeud(0, 12, 13, 0.0518, 0.028)); // 13
        treedist.push_back(noeud(3, 4, 14, 0.0885, 0.0780)); // 14
        treedist.push_back(noeud(13, 14, 15, 0.0623, 0.0308)); // 15
        treedist.push_back(noeud(5, 15, 16, 0.82, 0.58)); // 16
        treedist.push_back(noeud(6, 7, 17, 0.0056, 0.0140)); // 17
        treedist.push_back(noeud(16, 17, 18, 0.11, 0.4876)); // 18
        treedist.push_back(noeud(9, 10, 19, 0.3118, 0.358)); // 19
        treedist.push_back(noeud(11, 19, 20, 0.3598, 0.056)); // 20
        treedist.push_back(noeud(8, 18, 21, 0.6513, 0.078)); // 21
        treedist.push_back(noeud(20, 21, 22, 0.2092, 0.0150)); // 22
        noemax = 22;
    } else if (species == "eutherian") {
        // Arbre de ensembl epo 12 eutharian
        //
         // (
         // (
         // (
         // (
         // (
         // (
         // (
         // Pan_troglodytes:0.0067,
         // Homo_sapiens:0.0067):0.0022,
         // Gorilla_gorilla:0.0088):0.0097,
         // Pongo_abelii:0.0183):0.0143,
         // Macaca_mulatta:0.0375):0.0220,
         // Callithrix_jacchus:0.0661):0.0891,
         // (
         // Mus_musculus:0.0845,
         // Rattus_norvegicus:0.0916):0.2720):0.0206,
         // (
         // (
         // Equus_caballus:0.1094,
         // Canis_familiaris:0.1523):0.0107,
         // (
         // Sus_scrofa:0.0790,
         // Bos_taurus:0.1689):0.0202):0.0329);
        treedist.push_back(noeud(0, 1, 12, 0.0845, 0.0916)); // 12: mus and rat
        treedist.push_back(noeud(2, 3, 13, 0.0067, 0.0067)); // 13: pan and hom
        treedist.push_back(noeud(13, 4, 14, 0.0022, 0.0088)); // 14: 13 and gor
        treedist.push_back(noeud(14, 5, 15, 0.0097, 0.0183)); // 15: 14 and pon
        treedist.push_back(noeud(15, 6, 16, 0.0143, 0.0375)); // 16: 15 and mac
        treedist.push_back(noeud(16, 7, 17, 0.0220, 0.0661)); // 17: 16 and cal
//        treedist.push_back(noeud(8, 15, 16, 0.1477, 0.0796)); // 16 can and 15 !!!!!! OLD TREE
        treedist.push_back(noeud(12, 17, 18, 0.2720, 0.0891)); // 18: 12 and 17
        treedist.push_back(noeud(8, 9, 19, 0.1094, 0.1523)); // 19: equ and can
        treedist.push_back(noeud(10, 11, 20, 0.0790, 0.1689)); // 20: sus and bos
        treedist.push_back(noeud(19, 20, 21, 0.0107, 0.0202)); // 21: 19 and 20
        treedist.push_back(noeud(18, 21, 22, 0.0206, 0.0329)); // 22: 18 and 21
        noemax = 22;
    }
    vtransi.clear();
    if (evolutionary_model == 2) {
        for (unsigned int i = 0; i < noemax + 1; i++) {
            gsl_matrix * m = gsl_matrix_alloc(4, 4);
            vtransi.push_back(m);
        }
        instrates = gsl_matrix_calloc(4, 4);
        pnoe1out = gsl_vector_calloc(4);
        pnoe2out = gsl_vector_calloc(4);
        id = gsl_matrix_calloc(4, 4);
        for (unsigned int i = 0; i < 4; i++) {
            gsl_matrix_set(id, i, i, 1);
        }
        m1 = gsl_matrix_calloc(4, 4);
        m2 = gsl_matrix_calloc(4, 4);
        m3 = gsl_matrix_calloc(4, 4);
        m4 = gsl_matrix_calloc(4, 4);
        pij = gsl_matrix_calloc(4, 4);
        pijp = gsl_matrix_calloc(4, 4);
    } else if (evolutionary_model == 1) {
        proba2 = gsl_vector_alloc(4);
    }
}

int
speciestonum(string name)
{
    if (species == "droso") {
        if (name == "DroMel") return 0;
        else if (name == "DroSim") return 1;
        else if (name == "DroSec") return 2;
        else if (name == "DroYak") return 3;
        else if (name == "DroEre") return 4;
        else if (name == "DroAna") return 5;
        else if (name == "DroPse") return 6;
        else if (name == "DroPer") return 7;
        else if (name == "DroWil") return 8;
        else if (name == "DroVir") return 9;
        else if (name == "DroMoj") return 10;
        else if (name == "DroGri") return 11;
    } else if (species == "eutherian") {
        if (name == "MusMus") return 0;
        else if (name == "RatNor") return 1;
        else if (name == "PanTro") return 2;
        else if (name == "HomSap") return 3;
        else if (name == "GorGor") return 4;
        else if (name == "PonAbe") return 5;
        else if (name == "MacMul") return 6;
        else if (name == "CalJac") return 7;
        else if (name == "EquCab") return 8;
        else if (name == "CanFam") return 9;
        else if (name == "SusScr") return 10;
        else if (name == "BosTau") return 11;
    }
    return -1;
}

string numtospecies(int num)
{
    if (species == "droso") {
        if (num == 0) return "DroMel";
        else if (num == 1) return "DroSim";
        else if (num == 2) return "DroSec";
        else if (num == 3) return "DroYak";
        else if (num == 4) return "DroEre";
        else if (num == 5) return "DroAna";
        else if (num == 6) return "DroPse";
        else if (num == 7) return "DroPer";
        else if (num == 8) return "DroWil";
        else if (num == 9) return "DroVir";
        else if (num == 10) return "DroMoj";
        else if (num == 11) return "DroGri";
        else return "No name";
    } else if (species == "eutherian") {
        if (num == 0) return "MusMus";
        else if (num == 1) return "RatNor";
        else if (num == 2) return "PanTro";
        else if (num == 3) return "HomSap";
        else if (num == 4) return "GorGor";
        else if (num == 5) return "PonAbe";
        else if (num == 6) return "MacMul";
        else if (num == 7) return "CalJac";
        else if (num == 8) return "EquCab";
        else if (num == 9) return "CanFam";
        else if (num == 10) return "SusScr";
        else if (num == 11) return "BosTau";
        else return "No name";
    }
    return ("");
}

int
instant_rates(const gsl_vector * w, gsl_matrix * rates)
{
    double w0 = gsl_vector_get(w, 0);
    double w1 = gsl_vector_get(w, 1);
    double w2 = gsl_vector_get(w, 2);
    double w3 = 1 - w0 - w1 - w2;
    if (w0 < 0 || w1 < 0 || w2 < 0 || w3 < 0) {
        return -1;
    }
    double fat = proba_fixation_rel(w3 / w0);
    double fac = proba_fixation_rel(pa * w1 / pc / w0);
    double fag = proba_fixation_rel(pa * w2 / pc / w0);
    
    double fta = proba_fixation_rel(w0 / w3);
    double ftc = proba_fixation_rel(pa * w1 / pc / w3);
    double ftg = proba_fixation_rel(pa * w2 / pc / w3);
    
    double fca = proba_fixation_rel(pc * w0 / pa / w1);
    double fct = proba_fixation_rel(pc * w3 / pa / w1);
    double fcg = proba_fixation_rel(w2 / w1);
    
    double fga = proba_fixation_rel(pc * w0 / pa / w2);
    double fgt = proba_fixation_rel(pc * w3 / pa / w2);
    double fgc = proba_fixation_rel(w1 / w2);
    
    double prefact = 1.0 / (4 * kappa * pa * pc + 0.5);
    
    gsl_matrix_set(m1, 0, 0, -(pa * fat + pc * fac + pc * kappa * fag));
    gsl_matrix_set(m1, 0, 1, pa * fca);
    gsl_matrix_set(m1, 0, 2, pa * kappa * fga);
    gsl_matrix_set(m1, 0, 3, pa * fta);
    
    gsl_matrix_set(m1, 1, 0, pc * fac);
    gsl_matrix_set(m1, 1, 1, -(pa * fca + pa * kappa * fct + pc * fcg));
    gsl_matrix_set(m1, 1, 2, pc * fgc);
    gsl_matrix_set(m1, 1, 3, pc * kappa * ftc);
    
    gsl_matrix_set(m1, 2, 0, pc * kappa * fag);
    gsl_matrix_set(m1, 2, 1, pc * fcg);
    gsl_matrix_set(m1, 2, 2, -(pa * kappa * fga + pa * fgt + pc * fgc));
    gsl_matrix_set(m1, 2, 3, pc * ftg);
    
    gsl_matrix_set(m1, 3, 0, pa * fat);
    gsl_matrix_set(m1, 3, 1, pa * kappa * fct);
    gsl_matrix_set(m1, 3, 2, pa * fgt);
    gsl_matrix_set(m1, 3, 3, -(pa * fta + pc * kappa * ftc + pc * ftg));
    
    gsl_matrix_scale(m1, prefact * integr_step);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.5, m1, m1, 0.0, m2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1 / 3.0, m1, m2, 0.0, m3);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1 / 4.0, m1, m3, 0.0, m4);
    
    gsl_matrix_memcpy(rates, id);
    gsl_matrix_add(rates, m1);
    gsl_matrix_add(rates, m2);
    gsl_matrix_add(rates, m3);
    gsl_matrix_add(rates, m4);
    
    return 0;
}

int
func(double t, const double y[], double f[],
     void * params)
{
    gsl_matrix * rates = (gsl_matrix *)params;
    gsl_vector_const_view yview = gsl_vector_const_view_array(y, 4);
    gsl_vector_view fview = gsl_vector_view_array(f, 4);
    gsl_blas_dgemv(CblasNoTrans, 1.0, rates, &yview.vector, 0.0, &fview.vector);
    return GSL_SUCCESS;
}

int
jac(double t, const double y[], double * dfdy,
    double dfdt[], void * params)
{
    gsl_matrix_view dfdy_mat
        = gsl_matrix_view_array(dfdy, 4, 4);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix * rates = (gsl_matrix *)params;
    gsl_matrix_memcpy(m, rates);
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    dfdt[2] = 0.0;
    dfdt[3] = 0.0;
    return GSL_SUCCESS;
}

double
proba_fixation_rel(double ratio)
{
    if (ratio == 1.0) ratio = 1.0 + 1e-6;
    return ratio * log(ratio) / (ratio - 1.0);
}

void
initprobatree(const unsigned int pos, Motalign & ma, gsl_matrix * probatree)
{
    for (unsigned int i = 0; i < nbspecies; i++) {
        if (ma.matches[i]) {
            gsl_matrix_set(probatree, 0, i, 0);
            gsl_matrix_set(probatree, 1, i, 0);
            gsl_matrix_set(probatree, 2, i, 0);
            gsl_matrix_set(probatree, 3, i, 0);
            gsl_matrix_set(probatree, ma.alignseq[i][pos], i, 1);
        } else {
            gsl_matrix_set(probatree, 0, i, -1);
        }
    }
}


double
loglikely_column(const unsigned int pos, Motalign & ma, vpgslmat & vtrans, const gsl_vector * w)
{
    gsl_matrix * probatree = gsl_matrix_alloc(4, noemax + 1);
    gsl_vector * wfull = gsl_vector_alloc(4);
    initprobatree(pos, ma, probatree);
    gsl_vector_set(wfull, 0, gsl_vector_get(w, 0));
    gsl_vector_set(wfull, 1, gsl_vector_get(w, 1));
    gsl_vector_set(wfull, 2, gsl_vector_get(w, 2));
    gsl_vector_set(wfull, 3, 1 - gsl_vector_get(w, 0) - gsl_vector_get(w, 1) - gsl_vector_get(w, 2));
    for (ivnoe iv = treedist.begin(); iv != treedist.end(); iv++) {
        int n1 = iv->esp1;
        int n2 = iv->esp2;
        if (evolutionary_model == 2) {
            if (gsl_matrix_get(probatree, 0, n1) < -0.5 && gsl_matrix_get(probatree, 0, n2) < -0.5) {
                gsl_matrix_set(probatree, 0, iv->noe, -1.0);
            } else if (gsl_matrix_get(probatree, 0, n2) < -0.5) {
                gsl_vector_view pnoe1 = gsl_matrix_column(probatree, (*iv).esp1);
                gsl_vector_view pnoe = gsl_matrix_column(probatree, (*iv).noe);
                gsl_blas_dgemv(CblasTrans, 1.0, vtrans[2 * ((*iv).noe - nbspecies)], &pnoe1.vector, 0.0, &pnoe.vector);
            } else if (gsl_matrix_get(probatree, 0, n1) < -0.5) {
                gsl_vector_view pnoe2 = gsl_matrix_column(probatree, (*iv).esp2);
                gsl_vector_view pnoe = gsl_matrix_column(probatree, (*iv).noe);
                gsl_blas_dgemv(CblasTrans, 1.0, vtrans[2 * ((*iv).noe - nbspecies) + 1], &pnoe2.vector, 0.0, &pnoe.vector);
            } else {
                gsl_vector_view pnoe1 = gsl_matrix_column(probatree, (*iv).esp1);
                gsl_vector_view pnoe2 = gsl_matrix_column(probatree, (*iv).esp2);
                gsl_vector_view pnoe = gsl_matrix_column(probatree, (*iv).noe);
                gsl_blas_dgemv(CblasTrans, 1.0, vtrans[2 * ((*iv).noe - nbspecies)], &pnoe1.vector, 0.0, pnoe1out);
                gsl_blas_dgemv(CblasTrans, 1.0, vtrans[2 * ((*iv).noe - nbspecies) + 1], &pnoe2.vector, 0.0, pnoe2out);
                gsl_vector_memcpy(&pnoe.vector, pnoe1out);
                gsl_vector_mul(&pnoe.vector, pnoe2out);
            }
        } else if (evolutionary_model == 1) {
            double prox1 = iv->prox1;
            double prox2 = iv->prox2;
            if (gsl_matrix_get(probatree, 0, n1) < -0.5 && gsl_matrix_get(probatree, 0, n2) < -0.5) {
                gsl_matrix_set(probatree, 0, iv->noe, -1.0);
            } else if (gsl_matrix_get(probatree, 0, n2) < -0.5) {
                gsl_vector_view pnoe1 = gsl_matrix_column(probatree, (*iv).esp1);
                gsl_vector_view pnoe = gsl_matrix_column(probatree, (*iv).noe);
                double sum;
                gsl_blas_ddot(wfull, &pnoe1.vector, &sum);
                sum *= (1 - prox1);
                gsl_vector_set_all(&pnoe.vector, sum);
                gsl_blas_daxpy(prox1, &pnoe1.vector, &pnoe.vector);
            } else if (gsl_matrix_get(probatree, 0, n1) < -0.5) {
                gsl_vector_view pnoe2 = gsl_matrix_column(probatree, (*iv).esp2);
                gsl_vector_view pnoe = gsl_matrix_column(probatree, (*iv).noe);
                double sum;
                gsl_blas_ddot(wfull, &pnoe2.vector, &sum);
                sum *= (1 - prox2);
                gsl_vector_set_all(&pnoe.vector, sum);
                gsl_blas_daxpy(prox2, &pnoe2.vector, &pnoe.vector);
            } else {
                gsl_vector_view pnoe1 = gsl_matrix_column(probatree, iv->esp1);
                gsl_vector_view pnoe2 = gsl_matrix_column(probatree, iv->esp2);
                gsl_vector_view pnoe = gsl_matrix_column(probatree, iv->noe);
                double sum1, sum2;
                gsl_blas_ddot(wfull, &pnoe1.vector, &sum1);
                sum1 *= (1 - prox1);
                gsl_blas_ddot(wfull, &pnoe2.vector, &sum2);
                sum2 *= (1 - prox2);
                gsl_vector_set_all(&pnoe.vector, sum1);
                gsl_vector_set_all(proba2, sum2);
                gsl_blas_daxpy(prox1, &pnoe1.vector, &pnoe.vector);
                gsl_blas_daxpy(prox2, &pnoe2.vector, proba2);
                gsl_vector_mul(&pnoe.vector, proba2);
            }
        }
    }
    if (gsl_matrix_get(probatree, 0, noemax) < -0.5) cout << "error!!!\n";
    double w0 = gsl_vector_get(w, 0);
    double w1 = gsl_vector_get(w, 1);
    double w2 = gsl_vector_get(w, 2);
    double w3 = 1 - w0 - w1 - w2;
    //   cout << log(w0*gsl_matrix_get(probatree,0,noemax)+
    //         w1*gsl_matrix_get(probatree,1,noemax)+
    //         w2*gsl_matrix_get(probatree,2,noemax)+
    //         w3*gsl_matrix_get(probatree,3,noemax)) << endl;
    double res = log(w0 * gsl_matrix_get(probatree, 0, noemax) +
                     w1 * gsl_matrix_get(probatree, 1, noemax) +
                     w2 * gsl_matrix_get(probatree, 2, noemax) +
                     w3 * gsl_matrix_get(probatree, 3, noemax));
    gsl_matrix_free(probatree);
    gsl_vector_free(wfull);
    return res;
}

double
loglikely(const gsl_vector * w, void * params)
{
    void ** par = (void **) params;
    Motif & mot = *((Motif *)(par[0]));
    const unsigned int pos = *((const unsigned int *)(par[1]));
    gsl_matrix * pmattemp;
    if (evolutionary_model == 2) {
        if (instant_rates(w, instrates)) return 1e10;
        gsl_matrix_memcpy(pij, id);
        if (species == "droso") {
            for (unsigned int i = 1; i < 117; i++) {
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, instrates, pij, 0.0, pijp);
                pmattemp = pij;
                pij = pijp;
                pijp = pmattemp;
                //      gsl_matrix_memcpy(pij,pijp);
                if (i == 1) {
                    gsl_matrix_memcpy(vtransi[10], pij);
                } else if (i == 2) {
                    gsl_matrix_memcpy(vtransi[0], pij);
                    gsl_matrix_memcpy(vtransi[11], pij);
                } else if (i == 3) {
                    gsl_matrix_memcpy(vtransi[1], pij);
                    gsl_matrix_memcpy(vtransi[3], pij);
                } else if (i == 4) {
                    gsl_matrix_memcpy(vtransi[7], pij);
                } else if (i == 7) {
                    gsl_matrix_memcpy(vtransi[2], pij);
                    gsl_matrix_memcpy(vtransi[6], pij);
                } else if (i == 10) {
                    gsl_matrix_memcpy(vtransi[5], pij);
                    gsl_matrix_memcpy(vtransi[17], pij);
                } else if (i == 11) {
                    gsl_matrix_memcpy(vtransi[4], pij);
                } else if (i == 14) {
                    gsl_matrix_memcpy(vtransi[20], pij);
                } else if (i == 15) {
                    gsl_matrix_memcpy(vtransi[21], pij);
                } else if (i == 35) {
                    gsl_matrix_memcpy(vtransi[12], pij);
                    gsl_matrix_memcpy(vtransi[14], pij);
                } else if (i == 42) {
                    gsl_matrix_memcpy(vtransi[19], pij);
                } else if (i == 46) {
                    gsl_matrix_memcpy(vtransi[16], pij);
                } else if (i == 49) {
                    gsl_matrix_memcpy(vtransi[15], pij);
                } else if (i == 58) {
                    gsl_matrix_memcpy(vtransi[9], pij);
                } else if (i == 68) {
                    gsl_matrix_memcpy(vtransi[13], pij);
                } else if (i == 82) {
                    gsl_matrix_memcpy(vtransi[8], pij);
                } else if (i == 116) {
                    gsl_matrix_memcpy(vtransi[18], pij);
                }
            }
        } else if (species == "eutherian") {
            for (unsigned int i = 1; i < 28; i++) {
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, instrates, pij, 0.0, pijp);
                pmattemp = pij;
                pij = pijp;
                pijp = pmattemp;
                //      gsl_matrix_memcpy(pij,pijp);
                if (i == 1) {
                    gsl_matrix_memcpy(vtransi[2], pij);
                    gsl_matrix_memcpy(vtransi[3], pij);
                    gsl_matrix_memcpy(vtransi[4], pij);
                    gsl_matrix_memcpy(vtransi[13], pij);
                    gsl_matrix_memcpy(vtransi[14], pij);
                    gsl_matrix_memcpy(vtransi[15], pij);
                    gsl_matrix_memcpy(vtransi[19], pij);
                } else if (i == 2) {
                    gsl_matrix_memcpy(vtransi[5], pij);
                    gsl_matrix_memcpy(vtransi[16], pij);
                    gsl_matrix_memcpy(vtransi[18], pij);
                    gsl_matrix_memcpy(vtransi[20], pij);
                } else if (i == 3) {
                    gsl_matrix_memcpy(vtransi[21], pij);
                } else if (i == 4) {
                    gsl_matrix_memcpy(vtransi[6], pij);
                } else if (i == 7) {
                    gsl_matrix_memcpy(vtransi[7], pij);
                } else if (i == 8) {
                    gsl_matrix_memcpy(vtransi[0], pij);
                    gsl_matrix_memcpy(vtransi[10], pij);
                } else if (i == 9) {
                    gsl_matrix_memcpy(vtransi[17], pij);
                    gsl_matrix_memcpy(vtransi[1], pij);
                } else if (i == 11) {
                    gsl_matrix_memcpy(vtransi[8], pij);
                } else if (i == 15) {
                    gsl_matrix_memcpy(vtransi[9], pij);
                } else if (i == 17) {
                    gsl_matrix_memcpy(vtransi[11], pij);
                } else if (i == 27) {
                    gsl_matrix_memcpy(vtransi[12], pij);
                }
            }
        }
    }
    //   for (unsigned int i=0;i<4;i++){
    //      for (unsigned int j=0;j<4;j++){
    //         cout << gsl_matrix_get(vtransi[8],i,j) << " ";
    //      }
    //      cout << endl;
    //   }
    //   cout << endl;
    double logli = 0;
    for (ivma ima = mot.seqs.begin(); ima != mot.seqs.end(); ima++) {
        logli += loglikely_column(pos, *ima, vtransi, w);
    }
    double w0 = gsl_vector_get(w, 0);
    double w1 = gsl_vector_get(w, 1);
    double w2 = gsl_vector_get(w, 2);
    double w3 = 1 - w0 - w1 - w2;
    logli += alpha * (log(w0) + log(w3)) + beta * (log(w1) + log(w2));
    // logli+=(alpha-1.)*(log(w0)+log(w1))+(beta-1.)*(log(w2)+log(w3));
    return -logli;
}

double
likelyhood(vd x, void * params)
{
    gsl_vector * w;
    w = gsl_vector_alloc(3);
    gsl_vector_set(w, 0, x[0]);
    gsl_vector_set(w, 1, x[1]);
    gsl_vector_set(w, 2, x[2]);
    void ** par = (void **) params;
    Motif & mot = *((Motif *)(par[0]));
    const unsigned int pos = *((const unsigned int *)(par[1]));
    gsl_matrix * pmattemp;
    if (evolutionary_model == 2) {
        if (instant_rates(w, instrates)) return 1e10;
        gsl_matrix_memcpy(pij, id);
        if (species == "droso") {
            for (unsigned int i = 1; i < 117; i++) {
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, instrates, pij, 0.0, pijp);
                pmattemp = pij;
                pij = pijp;
                pijp = pmattemp;
                //      gsl_matrix_memcpy(pij,pijp);
                if (i == 1) {
                    gsl_matrix_memcpy(vtransi[10], pij);
                } else if (i == 2) {
                    gsl_matrix_memcpy(vtransi[0], pij);
                    gsl_matrix_memcpy(vtransi[11], pij);
                } else if (i == 3) {
                    gsl_matrix_memcpy(vtransi[1], pij);
                    gsl_matrix_memcpy(vtransi[3], pij);
                } else if (i == 4) {
                    gsl_matrix_memcpy(vtransi[7], pij);
                } else if (i == 7) {
                    gsl_matrix_memcpy(vtransi[2], pij);
                    gsl_matrix_memcpy(vtransi[6], pij);
                } else if (i == 10) {
                    gsl_matrix_memcpy(vtransi[5], pij);
                    gsl_matrix_memcpy(vtransi[17], pij);
                } else if (i == 11) {
                    gsl_matrix_memcpy(vtransi[4], pij);
                } else if (i == 14) {
                    gsl_matrix_memcpy(vtransi[20], pij);
                } else if (i == 15) {
                    gsl_matrix_memcpy(vtransi[21], pij);
                } else if (i == 35) {
                    gsl_matrix_memcpy(vtransi[12], pij);
                    gsl_matrix_memcpy(vtransi[14], pij);
                } else if (i == 42) {
                    gsl_matrix_memcpy(vtransi[19], pij);
                } else if (i == 46) {
                    gsl_matrix_memcpy(vtransi[16], pij);
                } else if (i == 49) {
                    gsl_matrix_memcpy(vtransi[15], pij);
                } else if (i == 58) {
                    gsl_matrix_memcpy(vtransi[9], pij);
                } else if (i == 68) {
                    gsl_matrix_memcpy(vtransi[13], pij);
                } else if (i == 82) {
                    gsl_matrix_memcpy(vtransi[8], pij);
                } else if (i == 116) {
                    gsl_matrix_memcpy(vtransi[18], pij);
                }
            }
        } else if (species == "eutherian") {
            for (unsigned int i = 1; i < 26; i++) {
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, instrates, pij, 0.0, pijp);
                pmattemp = pij;
                pij = pijp;
                pijp = pmattemp;
                //      gsl_matrix_memcpy(pij,pijp);
                if (i == 1) {
                    gsl_matrix_memcpy(vtransi[15], pij);
                    gsl_matrix_memcpy(vtransi[2], pij);
                    gsl_matrix_memcpy(vtransi[3], pij);
                    gsl_matrix_memcpy(vtransi[5], pij);
                    gsl_matrix_memcpy(vtransi[7], pij);
                } else if (i == 2) {
                    gsl_matrix_memcpy(vtransi[4], pij);
                    gsl_matrix_memcpy(vtransi[16], pij);
                } else if (i == 3) {
                    gsl_matrix_memcpy(vtransi[17], pij);
                } else if (i == 6) {
                    gsl_matrix_memcpy(vtransi[6], pij);
                } else if (i == 8) {
                    gsl_matrix_memcpy(vtransi[0], pij);
                    gsl_matrix_memcpy(vtransi[10], pij);
                    gsl_matrix_memcpy(vtransi[11], pij);
                    gsl_matrix_memcpy(vtransi[13], pij);
                    gsl_matrix_memcpy(vtransi[1], pij);
                } else if (i == 11) {
                    gsl_matrix_memcpy(vtransi[9], pij);
                    gsl_matrix_memcpy(vtransi[14], pij);
                } else if (i == 15) {
                    gsl_matrix_memcpy(vtransi[12], pij);
                } else if (i == 25) {
                    gsl_matrix_memcpy(vtransi[8], pij);
                }
            }
        }
    }
    //   for (unsigned int i=0;i<4;i++){
    //      for (unsigned int j=0;j<4;j++){
    //         cout << gsl_matrix_get(vtransi[8],i,j) << " ";
    //      }
    //      cout << endl;
    //   }
    //   cout << endl;
    double logli = 0;
    for (ivma ima = mot.seqs.begin(); ima != mot.seqs.end(); ima++) {
        logli += loglikely_column(pos, *ima, vtransi, w);
    }
    double w0 = gsl_vector_get(w, 0);
    double w1 = gsl_vector_get(w, 1);
    double w2 = gsl_vector_get(w, 2);
    double w3 = 1 - w0 - w1 - w2;
    logli += (alpha - 1.) * (log(w0) + log(w3)) + (beta - 1.) * (log(w1) + log(w2));
    return exp(logli);
}

double
loglikelyhood(vd x, void * params)
{
    gsl_vector * w;
    w = gsl_vector_alloc(3);
    gsl_vector_set(w, 0, x[0]);
    gsl_vector_set(w, 1, x[1]);
    gsl_vector_set(w, 2, x[2]);
    void ** par = (void **) params;
    Motif & mot = *((Motif *)(par[0]));
    const unsigned int pos = *((const unsigned int *)(par[1]));
    gsl_matrix * pmattemp;
    if (evolutionary_model == 2) {
        if (instant_rates(w, instrates)) return 1e10;
        gsl_matrix_memcpy(pij, id);
        if (species == "droso") {
            for (unsigned int i = 1; i < 117; i++) {
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, instrates, pij, 0.0, pijp);
                pmattemp = pij;
                pij = pijp;
                pijp = pmattemp;
                //      gsl_matrix_memcpy(pij,pijp);
                if (i == 1) {
                    gsl_matrix_memcpy(vtransi[10], pij);
                } else if (i == 2) {
                    gsl_matrix_memcpy(vtransi[0], pij);
                    gsl_matrix_memcpy(vtransi[11], pij);
                } else if (i == 3) {
                    gsl_matrix_memcpy(vtransi[1], pij);
                    gsl_matrix_memcpy(vtransi[3], pij);
                } else if (i == 4) {
                    gsl_matrix_memcpy(vtransi[7], pij);
                } else if (i == 7) {
                    gsl_matrix_memcpy(vtransi[2], pij);
                    gsl_matrix_memcpy(vtransi[6], pij);
                } else if (i == 10) {
                    gsl_matrix_memcpy(vtransi[5], pij);
                    gsl_matrix_memcpy(vtransi[17], pij);
                } else if (i == 11) {
                    gsl_matrix_memcpy(vtransi[4], pij);
                } else if (i == 14) {
                    gsl_matrix_memcpy(vtransi[20], pij);
                } else if (i == 15) {
                    gsl_matrix_memcpy(vtransi[21], pij);
                } else if (i == 35) {
                    gsl_matrix_memcpy(vtransi[12], pij);
                    gsl_matrix_memcpy(vtransi[14], pij);
                } else if (i == 42) {
                    gsl_matrix_memcpy(vtransi[19], pij);
                } else if (i == 46) {
                    gsl_matrix_memcpy(vtransi[16], pij);
                } else if (i == 49) {
                    gsl_matrix_memcpy(vtransi[15], pij);
                } else if (i == 58) {
                    gsl_matrix_memcpy(vtransi[9], pij);
                } else if (i == 68) {
                    gsl_matrix_memcpy(vtransi[13], pij);
                } else if (i == 82) {
                    gsl_matrix_memcpy(vtransi[8], pij);
                } else if (i == 116) {
                    gsl_matrix_memcpy(vtransi[18], pij);
                }
            }
        } else if (species == "eutherian") {
            //         for (unsigned int i=1;i<254;i++){
            //            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,instrates,pij,0.0,pijp);
            //            pmattemp=pij;
            //            pij=pijp;
            //            pijp=pmattemp;
            //            //      gsl_matrix_memcpy(pij,pijp);
            //            //      good approx, for integr_step=0.001
            //            if (i==5){
            //               gsl_matrix_memcpy(vtransi[15],pij);
            //            } else if (i==7){
            //               gsl_matrix_memcpy(vtransi[2],pij);
            //            } else if (i==8){
            //               gsl_matrix_memcpy(vtransi[3],pij);
            //            } else if (i==10){
            //               gsl_matrix_memcpy(vtransi[5],pij);
            //            } else if (i==12){
            //               gsl_matrix_memcpy(vtransi[7],pij);
            //            } else if (i==22){
            //               gsl_matrix_memcpy(vtransi[4],pij);
            //            } else if (i==23){
            //               gsl_matrix_memcpy(vtransi[16],pij);
            //            } else if (i==35){
            //               gsl_matrix_memcpy(vtransi[17],pij);
            //            } else if (i==59){
            //               gsl_matrix_memcpy(vtransi[6],pij);
            //            } else if (i==77){
            //               gsl_matrix_memcpy(vtransi[0],pij);
            //            } else if (i==80){
            //               gsl_matrix_memcpy(vtransi[10],pij);
            //               gsl_matrix_memcpy(vtransi[11],pij);
            //               gsl_matrix_memcpy(vtransi[13],pij);
            //            } else if (i==82){
            //               gsl_matrix_memcpy(vtransi[1],pij);
            //            } else if (i==107){
            //               gsl_matrix_memcpy(vtransi[9],pij);
            //            } else if (i==110){
            //               gsl_matrix_memcpy(vtransi[14],pij);
            //            } else if (i==148){
            //               gsl_matrix_memcpy(vtransi[12],pij);
            //            } else if (i==253){
            //               gsl_matrix_memcpy(vtransi[8],pij);
            //            }
            //            }
            // approx, integr_step=0.01
            for (unsigned int i = 1; i < 26; i++) {
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, instrates, pij, 0.0, pijp);
                pmattemp = pij;
                pij = pijp;
                pijp = pmattemp;
                //      gsl_matrix_memcpy(pij,pijp);
                if (i == 1) {
                    gsl_matrix_memcpy(vtransi[15], pij);
                    gsl_matrix_memcpy(vtransi[2], pij);
                    gsl_matrix_memcpy(vtransi[3], pij);
                    gsl_matrix_memcpy(vtransi[5], pij);
                    gsl_matrix_memcpy(vtransi[7], pij);
                } else if (i == 2) {
                    gsl_matrix_memcpy(vtransi[4], pij);
                    gsl_matrix_memcpy(vtransi[16], pij);
                } else if (i == 3) {
                    gsl_matrix_memcpy(vtransi[17], pij);
                } else if (i == 6) {
                    gsl_matrix_memcpy(vtransi[6], pij);
                } else if (i == 8) {
                    gsl_matrix_memcpy(vtransi[0], pij);
                    gsl_matrix_memcpy(vtransi[10], pij);
                    gsl_matrix_memcpy(vtransi[11], pij);
                    gsl_matrix_memcpy(vtransi[13], pij);
                    gsl_matrix_memcpy(vtransi[1], pij);
                } else if (i == 11) {
                    gsl_matrix_memcpy(vtransi[9], pij);
                    gsl_matrix_memcpy(vtransi[14], pij);
                } else if (i == 15) {
                    gsl_matrix_memcpy(vtransi[12], pij);
                } else if (i == 25) {
                    gsl_matrix_memcpy(vtransi[8], pij);
                }
            }
        }
    }
    //   for (unsigned int i=0;i<4;i++){
    //      for (unsigned int j=0;j<4;j++){
    //         cout << gsl_matrix_get(vtransi[8],i,j) << " ";
    //      }
    //      cout << endl;
    //   }
    //   cout << endl;
    double logli = 0;
    for (ivma ima = mot.seqs.begin(); ima != mot.seqs.end(); ima++) {
        logli += loglikely_column(pos, *ima, vtransi, w);
    }
    double w0 = gsl_vector_get(w, 0);
    double w1 = gsl_vector_get(w, 1);
    double w2 = gsl_vector_get(w, 2);
    double w3 = 1 - w0 - w1 - w2;
    logli += (alpha - 1.) * (log(w0) + log(w3)) + (beta - 1.) * (log(w1) + log(w2));
    return logli;
}
