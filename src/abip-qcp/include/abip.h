#ifndef ABIP_H_GUARD
#define ABIP_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include <string.h>
#include "amatrix.h"

/*-----special problem to solve-----------*/
enum problem_type{LASSO, SVM, QCP, SVMQP};

typedef struct ABIP_A_DATA_MATRIX ABIPMatrix;
typedef struct ABIP_LIN_SYS_WORK ABIPLinSysWork;

typedef struct ABIP_PROBLEM_DATA ABIPData;
typedef struct ABIP_SETTINGS ABIPSettings;
typedef struct ABIP_SOL_VARS ABIPSolution;
typedef struct ABIP_INFO ABIPInfo;
typedef struct ABIP_WORK ABIPWork;
typedef struct ABIP_ADAPTIVE_WORK ABIPAdaptWork;
typedef struct ABIP_RESIDUALS ABIPResiduals;
typedef struct ABIP_CONE ABIPCone;
typedef struct mkl_lin_sys MKLlinsys;
typedef struct solve_specific_problem spe_problem;

struct solve_specific_problem{

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
      ABIPMatrix *A;  //scaled original constraint matrix
      ABIPMatrix *Q;  //scaled original quadratic matrix
      abip_float *b;  //scaled reformulated vector for abip
      abip_float *c;  //scaled reformulated vector for abip
      /*-------------*/
      
      void (*scaling_data)(spe_problem *self, ABIPCone *k);
      void (*un_scaling_sol)(spe_problem *self, ABIPSolution *sol);
      void (*calc_residuals)(spe_problem *self, ABIPWork *w, ABIPResiduals *r, abip_int ipm_iter, abip_int admm_iter);
      abip_int (*init_spe_linsys_work)(spe_problem *self);
      abip_int (*solve_spe_linsys)(spe_problem *self, abip_float *b, abip_float *pcg_warm_start, abip_int iter, abip_float error_ratio);
      void (*free_spe_linsys_work)(spe_problem *self);
      void (*spe_A_times)(spe_problem *self, const abip_float *x, abip_float *y);//y += Ax, where A is the reformulated constraint matrix of ABIP
      void (*spe_AT_times)(spe_problem *self, const abip_float *x, abip_float *y);//y += A'x, where A is the reformulated constraint matrix of ABIP
      abip_float (*inner_conv_check)(spe_problem *self, ABIPWork *w);

};


/*  cols of data matrix A must be specified in this exact order 
    if change the definition of this order, remember to change the order in
    'solve_sub_problem' too
*/
struct ABIP_CONE {
    abip_int* q;    /* array of second-order cone constraints */
    abip_int qsize; /* length of SOC array */
    abip_int* rq;    /* array of rotated second-order cone constraints */
    abip_int rqsize; /* length of RSOC array */
    abip_int f;     /* length of free cone */
    abip_int z;     /* length of zero cone */
    abip_int l;     /* length of LP cone */

};


struct ABIP_PROBLEM_DATA 
{
      abip_int m;                                             
      abip_int n;   
      ABIPMatrix *A;
      ABIPMatrix *Q;                                       

      abip_float *b;                                         
      abip_float *c;  

      abip_float lambda;
      ABIPSettings *stgs;                             
};

struct ABIP_SETTINGS 
{
      abip_int normalize;      
      abip_int scale_E;
      abip_int scale_bc;                       
      abip_float scale;
      abip_float rho_x;
      abip_float rho_y;
      abip_float rho_tau; 

      abip_int max_ipm_iters;                     
      abip_int max_admm_iters;                 
      abip_float eps;
      abip_float eps_p;
      abip_float eps_d;
      abip_float eps_g;
      abip_float eps_inf;
      abip_float eps_unb;

      abip_float err_dif; //tol between max(dres,pres,dgap) of two consecutive inters                                  
      abip_float alpha; 
      abip_float cg_rate;  /* for indirect, tolerance goes down like (1/iter)^cg_rate: 2 */                                
      
      abip_int use_indirect;
      abip_int inner_check_period;
      abip_int outer_check_period;      

      abip_int verbose;    /* boolean, write out progress: 1 */      
      abip_int linsys_solver;   // 0:mkl_dss, 1:qdldl, 2:sparse cholesky, 3:pcg,  4:pardiso, 5:dense cholesky
      abip_int prob_type;       // 0:general_qp, 1:lasso, 2:svm, 3:qcp
      abip_float time_limit;  // in s
      abip_float psi;

      abip_int origin_scaling;
      abip_int ruiz_scaling;
      abip_int pc_scaling;

};

struct ABIP_SOL_VARS 
{
      abip_float *x; 
      abip_float *y; 
      abip_float *s;
};

struct ABIP_INFO 
{
      char status[32];                 
      abip_int status_val;             
      abip_int ipm_iter;            
      abip_int admm_iter;             
      
      abip_float pobj;
      abip_float dobj;
      abip_float res_pri;
      abip_float res_dual;
      abip_float rel_gap;
      abip_float res_infeas;
      abip_float res_unbdd;
  
      abip_float setup_time;         
      abip_float solve_time;      
      abip_float avg_linsys_time;   
      abip_float avg_cg_iters;
};


struct ABIP_WORK 
{
      abip_float sigma;                    
      abip_float gamma; 
      abip_float mu;                        
      abip_float beta;                      
      abip_float *u; 
      abip_float *v;
      abip_float *v_origin; 
      abip_float *u_t;  
      abip_float *rel_ut;
      abip_float nm_inf_b; 
      abip_float nm_inf_c;
      abip_int m; 
      abip_int n;
      ABIPMatrix *A;
      abip_float *r;
      abip_float a;

};

struct ABIP_RESIDUALS 
{
      abip_int last_ipm_iter;
      abip_int last_admm_iter;
      abip_float last_mu;
      
      abip_float res_pri;
      abip_float res_dual;
      abip_float rel_gap;
      abip_float res_infeas;
      abip_float res_unbdd;
      
      abip_float ct_x_by_tau;              
      abip_float bt_y_by_tau;        

      abip_float pobj;
      abip_float dobj;

      abip_float tau;
      abip_float kap;

      abip_float res_dif;
      abip_float error_ratio;

      abip_float Ax_b_norm;
      abip_float Qx_ATy_c_s_norm;
};


ABIPWork *ABIP(init)
(
      const ABIPData *d, 
      ABIPInfo *info,
      spe_problem *s,
      ABIPCone *c
) ;
abip_int ABIP(solve)
(
    ABIPWork *w, 
    const ABIPData *d,
    ABIPSolution *sol, 
    ABIPInfo *info,
    ABIPCone *c,
    spe_problem *s
);
void ABIP(finish)
(
    ABIPWork *w,
    spe_problem *spe
);

abip_int ABIP(main)(const ABIPData *d, ABIPSolution *sol, ABIPInfo *info);
const char *ABIP(version)(void);
abip_int abip
(
    const ABIPData* d,
    ABIPSolution* sol,
    ABIPInfo* info,
    ABIPCone *K
);


#ifdef __cplusplus
}
#endif
#endif
