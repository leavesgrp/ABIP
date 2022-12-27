#ifndef LASSO_CONFIG_H_GUARD
#define LASSO_CONFIG_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "abip.h"
#include "amatrix.h"
#include "linsys.h"
typedef struct Lasso lasso;

struct Lasso{

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
   
    void (*scaling_data)(lasso *self, ABIPCone *k);
    void (*un_scaling_sol)(lasso *self, ABIPSolution *sol);
    void (*calc_residuals)(lasso *self, ABIPWork *w, ABIPResiduals *r, abip_int ipm_iter, abip_int admm_iter);
    abip_int (*init_spe_linsys_work)(lasso *self);
    abip_int (*solve_spe_linsys)(lasso *self, abip_float *b, abip_float *pcg_warm_start, abip_int iter, abip_float error_ratio);
    void (*free_spe_linsys_work)(lasso *self);
    void (*spe_A_times)(lasso *self, const abip_float *x, abip_float *y);//y += Ax, where A is the reformulated constraint matrix of ABIP
    void (*spe_AT_times)(lasso *self, const abip_float *x, abip_float *y);//y += A'x, where A is the reformulated constraint matrix of ABIP
    abip_float (*inner_conv_check)(lasso *self, ABIPWork *w);
    /*--------------------------------*/

    abip_float lambda;
    /*-----scaling-------*/
    abip_float *D_hat;
    abip_float *D;
    abip_float *E;
    abip_float sc_b;
    abip_float sc_c;
    abip_float sc;
    abip_float sc_cone1;
    abip_float sc_cone2;
    /*-------------*/

};

abip_int init_lasso(lasso **gen_lasso, ABIPData *d, ABIPSettings *stgs);

void scaling_lasso_data(lasso *self, ABIPCone *k);
void un_scaling_lasso_sol(lasso *self, ABIPSolution *sol);
void calc_lasso_residuals(lasso *self, ABIPWork *w, ABIPResiduals *r, abip_int ipm_iter, abip_int admm_iter);
abip_int init_lasso_linsys_work(lasso *self);
abip_int solve_lasso_linsys(lasso *self, abip_float *b, abip_float *pcg_warm_start, abip_int iter, abip_float error_ratio);
void free_lasso_linsys_work(lasso *self);
void lasso_A_times(lasso *self, const abip_float *x, abip_float *y);//y += Ax, where A is the reformulated constraint matrix of ABIP
void lasso_AT_times(lasso *self, const abip_float *x, abip_float *y);//y += A'x, where A is the reformulated constraint matrix of ABIP
abip_float lasso_inner_conv_check(lasso *self, ABIPWork *w);

#ifdef __cplusplus
}
#endif
#endif