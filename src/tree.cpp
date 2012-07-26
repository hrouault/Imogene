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

extern "C" {
    /* Tridiagonal reduction of a packed symmetric matrix, real single
     * precision
     */
    int dsptrd_(char *uplo, int *n, double *ap, double *d__, double *e,
                double *tau, int *info);

    /* Multiplication of matrix after reduction, real single precision
     */
    int dopmtr_(char *side, char *uplo, char *trans, int *m, int *n, double *ap,
                double *tau, double *c__, int *ldc, double *work, int *info);

    /* Eigenvalue decomposition of a tridiagonal matrix
     */
    int dstegr_(char *jobz, char *range, int *n, double *d__, double *e,
                double *vl, double *vu, int *il, int *iu, double *abstol, 
	            int *m, double *eigs, double *z__, int *ldz, int *isuppz,
                double * work, int *lwork, int *iwork, int *liwork, int *info);

    /* transform vector into another vector
     */
    int daxpy_(int *n, double *da, double *dx, int *incx, double *dy, int *incy);

    /* dot product
     */
    double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);
 
    /* Multiply a matrix by a vector
     */
    int dgemv_(char *trans, int *m, int *n, double * alpha, double *a,
               int *lda, double *x, int *incx, double *beta, double *y,
               int *incy);

    /* Multiply two matrices
     */
    int dgemm_(const char *transa, const char *transb, int *l, int *n, int *m,
               double *alpha, const void *a, int *lda, void *b, int *ldb,
               double *beta, void *c, int *ldc);
}

using namespace std;

vnoe treedist;

double pa;
double pc;

const double kappa = 2.0;

double w[4];
double rates[4 * 4];

double fat, fac, fag, fta, ftc, ftg, fca, fct, fcg, fga, fgt, fgc;

noeud::noeud(int e1, int e2, int n, double p1, double p2)
{
    esp1 = e1;
    esp2 = e2;
    noe = n;
    if (evolutionary_model == 2) {
        dist1 = p1 / (4 * kappa * pa * pc + 0.5);;
        dist2 = p2 / (4 * kappa * pa * pc + 0.5);;
        prox1 = 0;
        prox2 = 0;
        transi1 = new double[4 * 4];
        transi2 = new double[4 * 4];
    } else if (evolutionary_model == 1) {
        double correction = 0.5 + 4.0 * conca * concc;
        prox1 = exp(-p1 / correction);
        prox2 = exp(-p2 / correction);
    }
}

gsl_matrix * instrates;
gsl_vector * pnoe1out;
gsl_vector * pnoe2out;
gsl_matrix * id;

gsl_vector * proba2;
unsigned int noemax;

vvd
transition_felsen(double dist)
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
    vvd M = transition_felsen(dist);
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
evolvedist(vd probs, double dist)
{
    vd pf(4, 0.);
    vvd M;
    if (evolutionary_model == 1){
        M = transition_felsen(dist);
    } else {
        // !!! To be changed to use the new exponentiation
        //M = transition_halpern(freqs,dist);
    }
    for (unsigned int row = 0; row < 4; row++) {
        for (unsigned int k = 0; k < 4; k++) {
            pf[row] += M[k][row] * probs[k];
        }
    }
    return pf;
}

void   // Currently used phylogenetic tree for drosophilae :  Heger and Pontig, 2007
inittreedist()
{
    pa = conca;
    pc = concc;
    treedist.clear();
    if (species == "droso") {
        treedist.push_back(noeud(1,2,12,0.02,0.03)); // 12                           
        treedist.push_back(noeud(0,12,13,0.07,0.03)); // 13                        
        treedist.push_back(noeud(3,4,14,0.11,0.10)); // 14                         
        treedist.push_back(noeud(13,14,15,0.07,0.04)); // 15                       
        treedist.push_back(noeud(5,15,16,0.82,0.58)); // 16                        
        treedist.push_back(noeud(6,7,17,0.01,0.02)); // 17                         
        treedist.push_back(noeud(16,17,18,0.35,0.68)); // 18                       
        treedist.push_back(noeud(9,10,19,0.35,0.49)); // 19                        
        treedist.push_back(noeud(11,19,20,0.46,0.10)); // 20                       
        treedist.push_back(noeud(8,20,21,1.16,0.42)); // 21                        
        treedist.push_back(noeud(18,21,22,0.14,0.15)); // 22 
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
    if (evolutionary_model == 1) {
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
check_freq_matrix(const double * w_full)
{
    if (w_full[0] < 0 || w_full[1] < 0 || w_full[2] < 0 || w_full[3] < 0) {
        return -1;
    }
    return 0;
}

int
transition_rates_halpern(const double * w_full, double * rates)
{
    double w0 = w_full[0];
    double w1 = w_full[1];
    double w2 = w_full[2];
    double w3 = w_full[3];

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
    
    rates[0] = -(pa * fat + pc * fac + pc * kappa * fag);
    rates[1] =  pa * fca;
    rates[2] =  pa * kappa * fga;
    rates[3] =  pa * fta;

    rates[4 + 0] =  pc * fac;
    rates[4 + 1] =  -(pa * fca + pa * kappa * fct + pc * fcg);
    rates[4 + 2] =  pc * fgc;
    rates[4 + 3] =  pc * kappa * ftc;

    rates[8 + 0] =  pc * kappa * fag;
    rates[8 + 1] =  pc * fcg;
    rates[8 + 2] =  -(pa * kappa * fga + pa * fgt + pc * fgc);
    rates[8 + 3] =  pc * ftg;

    rates[12 + 0] =  pa * fat;
    rates[12 + 1] =  pa * kappa * fct;
    rates[12 + 2] =  pa * fgt;
    rates[12 + 3] =  -(pa * fta + pc * kappa * ftc + pc * ftg);
    
    return EXIT_SUCCESS;
}

int
printmat(double * m)
{
    for (unsigned int i = 0 ; i < 4 ; i++){
        for (unsigned int j = 0 ; j < 4 ; j++){
            cout << m[i * 4 + j] << " ";
        }
        cout << endl;
    }
    return EXIT_SUCCESS;
}

int
update_transi_halpern(double * w_full)
{
    // The linalg calculus consists in:
    // construct V = L^-1 T L with L_ij = sqrt(pi_i) delta_ij
    // (V is symmetric)
    //
    // Then diagonalise V:
    // V = P A P^t
    //
    // A = L^-1 P A P^t L
    //

    int n = 4;

    transition_rates_halpern(w_full, rates);

    double * l = new double[n];
    double * lm1 = new double[n];
    double * v = new double[n * (n + 1) / 2];

    // Build l
    for (unsigned int i = 0 ; i < n ; i++){
        l[i] = sqrt(w_full[i]);
        lm1[i] = 1 / l[i];
    }

    // Build upper part of V, in packed form, linewise
    v[0] = rates[0];
    v[1] = l[1] * lm1[0] * rates[1];
    v[2] = l[2] * lm1[0] * rates[2];
    v[3] = l[3] * lm1[0] * rates[3];

    v[4] = rates[4 + 1];
    v[5] = l[2] * lm1[1] * rates[4 + 2];
    v[6] = l[3] * lm1[1] * rates[4 + 3];

    v[7] = rates[8 + 2];
    v[8] = l[3] * lm1[2] * rates[8 + 3];

    v[9] = rates[12 + 3];

    // Diagonalize V
    char uplo = 'L'; // the upper part is stored but fortran conventions are columnwise
    double d[n];
    double e[n];
    double tau[n - 1];
    int info;

    dsptrd_(&uplo, &n, v, d, e, tau, &info);

    char job = 'V'; // Compute eigenvalues AND eigenvectors
    char range = 'A'; // Compute all the eigenvalues

    double abstol = 1e-4;

    int m;
    double eigs[n];
    int ldz = n;
    double z[ldz * n];

    int isuppz[2 * n];

    int lwork = 128; // Hard coded but sufficiently large 
    double work[lwork];
    int liwork = 128;
    int iwork[liwork];

    dstegr_(&job, &range, &n, d, e, NULL, NULL, NULL, NULL, &abstol, &m, eigs,
            z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);

    // Compute Lq and qm1lm1
    double lq[n * n];
    double qm1lm1[n * n];


    double zcop[n * n];
    double zt[n * n];
    for (unsigned int i = 0 ; i < n * n ; i++){
        zcop[i] = z[i];
    }
    for (unsigned int i = 0 ; i < n ; i++){
        for (unsigned int j = 0 ; j < n ; j++){
            zt[i * n + j] = z[j * n + i];
        }
    }
    char side = 'L';
    char trans = 'N';
    double work2[n];
    dopmtr_(&side, &uplo, &trans, &n, &n, v, tau, zcop, &n, work2, &info);

    side = 'R';
    trans = 'T';
    dopmtr_(&side, &uplo, &trans, &n, &n, v, tau, zt, &n, work2, &info);

    // We now have to compute LQ and Q^-1L^-1
    for (unsigned int i = 0 ; i < n ; i++){
        for (unsigned int j = 0 ; j < n ; j++){
            lq[i * n + j] = l[j] * zcop[i * n + j]; // fortran convention
            qm1lm1[i * n + j] = lm1[i] * zt[i * n + j]; // fortran convention
        }
    }

    double expd1[4];
    double expd2[4];

    double dum1[4 * 4];
    double dum2[4 * 4];
    for (ivnoe iv = treedist.begin(); iv != treedist.end(); iv++) {
        for (unsigned int i = 0 ; i < 3 ; i++){
            double dist = iv->dist1;
            expd1[i] = exp(dist * eigs[i]);

            dist = iv->dist2;
            expd2[i] = exp(dist * eigs[i]);
        }
        expd1[3] = 1.0;
        expd2[3] = 1.0;

        for (unsigned int i = 0 ; i < n ; i++){
            for (unsigned int j = 0 ; j < n ; j++){
                dum1[i * n + j] = lq[i * n + j] * expd1[i];
                dum2[i * n + j] = lq[i * n + j] * expd2[i];
            }
        }
        char trans = 'T';
        double alpha = 1.0;
        double beta = 0.0;
        dgemm_(&trans, &trans, &n, &n, &n, &alpha, qm1lm1 , &n, dum1, &n,
               &beta, iv -> transi1, &n );
        dgemm_(&trans, &trans, &n, &n, &n, &alpha, qm1lm1 , &n, dum2, &n,
               &beta, iv -> transi2, &n );
    }

    return EXIT_SUCCESS;
}

double
proba_fixation_rel(double ratio)
{
    if (ratio == 1.0) ratio = 1.0 + 1e-6;
    return ratio * log(ratio) / (ratio - 1.0);
}

void
initprobatree(const unsigned int pos, Motalign & ma, double * probastree)
{
    for (unsigned int i = 0; i < nbspecies; i++) {
        double * node = &probastree[i * 4];
        if (ma.matches[i]) {
            for (unsigned int j = 0 ; j < 4 ; j++){
                node[j] = 0;
            }
            node[ma.alignseq[i][pos]] = 1;
        } else {
            node[0] = -1;
        }
    }
}

double
loglikely_column(const unsigned int pos, Motalign & ma)
{
    // Probablitity vector for each node of the tree
    double probastree[4 * (noemax + 1)];

    initprobatree(pos, ma, probastree);

    int inc = 1;
    int n = 4;

    double dum[4];

    for (ivnoe iv = treedist.begin(); iv != treedist.end(); iv++) {
        int n1 = iv->esp1;
        int n2 = iv->esp2;
        double * sourc1 = &probastree[n1 * 4];
        double * sourc2 = &probastree[n2 * 4];
        double * targ = &probastree[iv -> noe * 4];

        if (evolutionary_model == 2) {
            char trans = 'N'; // backward should impose the transpose but
                              // fortran convention add another transpose
            double alpha = 1.0;
            double beta = 0.0;

            if (*sourc1 < -0.5 && *sourc2 < -0.5){
                *targ = -1.0;
            } else if (*sourc2 < -0.5){
                dgemv_(&trans, &n, &n, &alpha, iv -> transi1, &n, sourc1, &inc,
                       &beta, targ, &inc);
            } else if (*sourc1 < -0.5){
                dgemv_(&trans, &n, &n, &alpha, iv -> transi2, &n, sourc2, &inc,
                       &beta, targ, &inc);
            } else {
                dgemv_(&trans, &n, &n, &alpha, iv -> transi1, &n, sourc1, &inc,
                       &beta, targ, &inc);
                dgemv_(&trans, &n, &n, &alpha, iv -> transi2, &n, sourc2, &inc,
                       &beta, dum, &inc);
                for (unsigned int i = 0 ; i < 4 ; i++){
                    targ[i] *= dum[i];
                }
            }

        } else if (evolutionary_model == 1) {
            double prox1 = iv->prox1;
            double prox2 = iv->prox2;
            double sum;
            if (*sourc1 < -0.5 && *sourc2 < -0.5) {
                *targ = -1.0;
            } else if (*sourc2 < -0.5) {
                sum = (1 - prox1) * ddot_(&n, w, &inc, sourc1, &inc);
                for (unsigned int i = 0 ; i < 4 ; i++){
                    targ[i] = sum;
                }
                daxpy_(&n, &prox1, sourc1, &inc, targ, &inc);
            } else if (*sourc1 < -0.5) {
                sum = (1 - prox2) * ddot_(&n, w, &inc, sourc2, &inc);
                for (unsigned int i = 0 ; i < 4 ; i++){
                    targ[i] = sum;
                }
                daxpy_(&n, &prox2, sourc2, &inc, targ, &inc);
            } else {
                sum = (1 - prox1) * ddot_(&n, w, &inc, sourc1, &inc);
                for (unsigned int i = 0 ; i < 4 ; i++){
                    targ[i] = sum;
                }
                daxpy_(&n, &prox1, sourc1, &inc, targ, &inc);
                sum = (1 - prox2) * ddot_(&n, w, &inc, sourc2, &inc);
                for (unsigned int i = 0 ; i < 4 ; i++){
                    dum[i] = sum;
                }
                daxpy_(&n, &prox2, sourc2, &inc, dum, &inc);
                for (unsigned int i = 0 ; i < 4 ; i++){
                    targ[i] *= dum[i];
                }
            }
        }
    }
    if (probastree[4 * noemax] < -0.5) cout << "error!!!\n";

    double * lnode = &probastree[4 * noemax];
    double res = log(w[0] * lnode[0] + w[1] * lnode[1] + w[2] * lnode[2] +
                     w[3] * lnode[3]);

    return res;
}

double
loglikelyhood(void * params)
{
    void ** par = (void **) params;
    Motif & mot = *((Motif *)(par[0]));
    const unsigned int pos = *((const unsigned int *)(par[1]));

    double logli = 0;
    for (ivma ima = mot.seqs.begin(); ima != mot.seqs.end(); ima++) {
        logli += loglikely_column(pos, *ima);
    }
    return logli;
}

int
initw(const double * w_pack){
    w[0] = w_pack[0];
    w[1] = w_pack[1];
    w[2] = w_pack[2];
    w[3] = 1.0 - w[0] - w[1] - w[2];
}

double
posterior_priormax(const gsl_vector * w_pack, void * params)
{
    double w_pack2[3]; 
    for (unsigned int i = 0 ; i < 3 ; i++){
        w_pack2[i] = gsl_vector_get(w_pack, i);
    }
    initw(w_pack2);
    double logli = loglikelyhood(params);
    logli += alpha * (log(w[0]) + log(w[3])) + beta * (log(w[1]) + log(w[2]));
    return -logli;
}

double
posterior_priormean(vd w_pack, void * params)
{
    double w_pack2[3]; 
    for (unsigned int i = 0 ; i < 3 ; i++){
        w_pack2[i] = w_pack[i];
    }
    initw(w_pack2);
    double logli = loglikelyhood(params);
    logli += alpha * (log(w[0]) + log(w[3])) + beta * (log(w[1]) + log(w[2]));
    return -logli;
}

