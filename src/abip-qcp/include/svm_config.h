#ifndef SVM_CONFIG_H_GUARD
#define SVM_CONFIG_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "abip.h"
#include "amatrix.h"
#include "linsys.h"
typedef struct Svm svm;

struct Svm{

    /*----------common parts-------------*/
    enum problem_type pro_type;
    abip_int m;  //rows of input data A
    abip_int n;  //cols of input data A
    abip_int p;  //rows of ABIP constraint matrix A
    abip_int q;  //cols of ABIP constraint matrix A
    ABIPLinSysWork *L;
    ABIPSettings *stgs;
    ABIPData *data; //original data
    abip_float sparsity;

    abip_float *rho_dr; // non-identity DR scaling

    /* scaled data */
    ABIPMatrix *A;
    ABIPMatrix *Q;
    abip_float *b;
    abip_float *c;
    /*-------------*/
   
    void (*scaling_data)(svm *self, ABIPCone *k);
    void (*un_scaling_sol)(svm *self, ABIPSolution *sol);
    void (*calc_residuals)(svm *self, ABIPWork *w, ABIPResiduals *r, abip_int ipm_iter, abip_int admm_iter);
    abip_int (*init_spe_linsys_work)(svm *self);
    abip_int (*solve_spe_linsys)(svm *self, abip_float *b, abip_float *pcg_warm_start, abip_int iter, abip_float error_ratio);
    void (*free_spe_linsys_work)(svm *self);
    void (*spe_A_times)(svm *self, const abip_float *x, abip_float *y);//y += Ax, where A is the reformulated constraint matrix of ABIP
    void (*spe_AT_times)(svm *self, const abip_float *x, abip_float *y);//y += A'x, where A is the reformulated constraint matrix of ABIP
    abip_float (*inner_conv_check)(svm *self, ABIPWork *w);
    /*--------------------------------*/

    abip_float lambda; // 'C' in svm
    /*-----scaling-------*/
    abip_float *sc_D;
    abip_float *sc_E;
    abip_float *sc_F;
    abip_float sc_b;
    abip_float sc_c;
    abip_float sc;
    abip_float sc_cone1;
    abip_float sc_cone2;
    /*-------------*/
    ABIPMatrix *wA;
    abip_float *wy;
    abip_float *wB;
    abip_float *wC;
    abip_float *wD;
    abip_float *wE;
    abip_float *wF;
    abip_float *wG;
    abip_float *wH;
    ABIPMatrix *wX;
};

abip_int init_svm(svm **gen_svm, ABIPData *d, ABIPSettings *stgs);

void scaling_svm_data(svm *self, ABIPCone *k);
void un_scaling_svm_sol(svm *self, ABIPSolution *sol);
void calc_svm_residuals(svm *self, ABIPWork *w, ABIPResiduals *r, abip_int ipm_iter, abip_int admm_iter);
abip_int init_svm_linsys_work(svm *self);
abip_int solve_svm_linsys(svm *self, abip_float *b, abip_float *pcg_warm_start, abip_int iter, abip_float error_ratio);
void free_svm_linsys_work(svm *self);
void svm_A_times(svm *self, const abip_float *x, abip_float *y);//y += Ax, where A is the reformulated constraint matrix of ABIP
void svm_AT_times(svm *self, const abip_float *x, abip_float *y);//y += A'x, where A is the reformulated constraint matrix of ABIP
abip_float svm_inner_conv_check(svm *self, ABIPWork *w);

#ifdef __cplusplus
}
#endif
#endif