#ifndef QCP_CONFIG_H_GUARD
#define QCP_CONFIG_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "abip.h"
#include "amatrix.h"
#include "linsys.h"
typedef struct qcp qcp;

struct qcp{

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
   
    void (*scaling_data)(qcp *self, ABIPCone *k);
    void (*un_scaling_sol)(qcp *self, ABIPSolution *sol);
    void (*calc_residuals)(qcp *self, ABIPWork *w, ABIPResiduals *r, abip_int ipm_iter, abip_int admm_iter);
    abip_int (*init_spe_linsys_work)(qcp *self);
    abip_int (*solve_spe_linsys)(qcp *self, abip_float *b, abip_float *pcg_warm_start, abip_int iter, abip_float error_ratio);
    void (*free_spe_linsys_work)(qcp *self);
    void (*spe_A_times)(qcp *self, const abip_float *x, abip_float *y);//y += Ax, where A is the reformulated constraint matrix of ABIP
    void (*spe_AT_times)(qcp *self, const abip_float *x, abip_float *y);//y += A'x, where A is the reformulated constraint matrix of ABIP
    abip_float (*inner_conv_check)(qcp *self, ABIPWork *w);
    /*--------------------------------*/

    /*-----scaling-------*/
    abip_float *D;
    abip_float *E;
    abip_float sc_b;
    abip_float sc_c;
    /*-------------*/

};

abip_int init_qcp(qcp **QCP, ABIPData *d, ABIPSettings *stgs);

void scaling_qcp_data(qcp *self, ABIPCone *k);
void un_scaling_qcp_sol(qcp *self, ABIPSolution *sol);
void calc_qcp_residuals(qcp *self, ABIPWork *w, ABIPResiduals *r, abip_int ipm_iter, abip_int admm_iter);
abip_int init_qcp_linsys_work(qcp *self);
abip_int solve_qcp_linsys(qcp *self, abip_float *b, abip_float *pcg_warm_start, abip_int iter, abip_float error_ratio);
void free_qcp_linsys_work(qcp *self);
void qcp_A_times(qcp *self, const abip_float *x, abip_float *y);//y += Ax, where A is the reformulated constraint matrix of ABIP
void qcp_AT_times(qcp *self, const abip_float *x, abip_float *y);//y += A'x, where A is the reformulated constraint matrix of ABIP
abip_float qcp_inner_conv_check(qcp *self, ABIPWork *w);
#ifdef __cplusplus
}
#endif
#endif