#ifndef SVM_QP_CONFIG_H_GUARD
#define SVM_QP_CONFIG_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "abip.h"
#include "amatrix.h"
#include "linsys.h"
typedef struct SVMqp svmqp;

struct SVMqp{

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
   
    void (*scaling_data)(svmqp *self, ABIPCone *k);
    void (*un_scaling_sol)(svmqp *self, ABIPSolution *sol);
    void (*calc_residuals)(svmqp *self, ABIPWork *w, ABIPResiduals *r, abip_int ipm_iter, abip_int admm_iter);
    abip_int (*init_spe_linsys_work)(svmqp *self);
    abip_int (*solve_spe_linsys)(svmqp *self, abip_float *b, abip_float *pcg_warm_start, abip_int iter, abip_float error_ratio);
    void (*free_spe_linsys_work)(svmqp *self);
    void (*spe_A_times)(svmqp *self, const abip_float *x, abip_float *y);//y += Ax, where A is the reformulated constraint matrix of ABIP
    void (*spe_AT_times)(svmqp *self, const abip_float *x, abip_float *y);//y += A'x, where A is the reformulated constraint matrix of ABIP
    abip_float (*inner_conv_check)(svmqp *self, ABIPWork *w);
    /*--------------------------------*/

    abip_float lambda; 
    /*-----scaling-------*/
    abip_float *D;
    abip_float *E;
    abip_float *F;
    abip_float *H;
    abip_float sc_b;
    abip_float sc_c;
};

abip_int init_svmqp(svmqp **gen_svm, ABIPData *d, ABIPSettings *stgs);

void scaling_svmqp_data(svmqp *self, ABIPCone *k);
void un_scaling_svmqp_sol(svmqp *self, ABIPSolution *sol);
void calc_svmqp_residuals(svmqp *self, ABIPWork *w, ABIPResiduals *r, abip_int ipm_iter, abip_int admm_iter);
abip_int init_svmqp_linsys_work(svmqp *self);
abip_int solve_svmqp_linsys(svmqp *self, abip_float *b, abip_float *pcg_warm_start, abip_int iter, abip_float error_ratio);
void free_svmqp_linsys_work(svmqp *self);
void svmqp_A_times(svmqp *self, const abip_float *x, abip_float *y);//y += Ax, where A is the reformulated constraint matrix of ABIP
void svmqp_AT_times(svmqp *self, const abip_float *x, abip_float *y);//y += A'x, where A is the reformulated constraint matrix of ABIP
abip_float svmqp_inner_conv_check(svmqp *self, ABIPWork *w);

#ifdef __cplusplus
}
#endif
#endif