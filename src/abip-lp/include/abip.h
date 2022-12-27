#ifndef ABIP_H_GUARD
#define ABIP_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include <string.h>

typedef struct ABIP_A_DATA_MATRIX ABIPMatrix;
typedef struct ABIP_LIN_SYS_WORK ABIPLinSysWork;

typedef struct ABIP_PROBLEM_DATA ABIPData;
typedef struct ABIP_SETTINGS ABIPSettings;
typedef struct ABIP_SOL_VARS ABIPSolution;
typedef struct ABIP_INFO ABIPInfo;
typedef struct ABIP_SCALING ABIPScaling;
typedef struct ABIP_WORK ABIPWork;
typedef struct ABIP_ADAPTIVE_WORK ABIPAdaptWork;
typedef struct ABIP_RESIDUALS ABIPResiduals;

struct ABIP_PROBLEM_DATA 
{
      abip_int m;                                             
      abip_int n;                                              
      ABIPMatrix *A;                                       

      abip_float *b;                                         
      abip_float *c;  
      abip_float sp; 

      ABIPSettings *stgs;                             
};

struct ABIP_SETTINGS 
{
      abip_int normalize;    
      abip_int pfeasopt;                         
      abip_float scale;                                  
      abip_float rho_y; 
      abip_float sparsity_ratio;                                  

      abip_int max_ipm_iters;                     
      abip_int max_admm_iters;
      abip_float max_time;

      abip_float eps;                                    
      abip_float alpha; 
      abip_float cg_rate;  /* for indirect, tolerance goes down like (1/iter)^cg_rate: 2 */                                
      
      abip_int adaptive;                              
      abip_float eps_cor;                             
      abip_float eps_pen;     

      abip_float dynamic_sigma; /* Dynamic strategy A to update the barrier parameter */

      abip_float dynamic_x; /* Dynamic strategy B to update the barrier parameter */
      abip_float dynamic_eta;

      abip_int   restart_fre; /* Frequency of performing restart*/
      abip_int   restart_thresh; /* Threshold of restart */
  
      abip_int verbose;    /* boolean, write out progress: 1 */                         
      abip_int warm_start; /* boolean, warm start (put initial guess in ABIPSolution struct): 0 */                         

      abip_int adaptive_lookback; 


      abip_int origin_rescale;            /* boolean, use the origin_rescale or not */
      abip_int pc_ruiz_rescale;           /* boolean, use the pc_rescale or not */
      abip_int qp_rescale;                /* boolean, use the qp rescale or not*/
      abip_int ruiz_iter;                 /* int, number of ruiz rescaling */
      abip_int hybrid_mu;                 /* boolean, use the hybrid_mu strategy or not */
      abip_float hybrid_thresh;             /* float, control when to use the LOQO mu strategy, only used when hybrid_mu = 1 */
      abip_float dynamic_sigma_second;    /* float, control the dynamic sigma when the mu strategy becomes LOQO */
      abip_int half_update;              /* boolean, use the half update technique or not*/
      abip_int avg_criterion;           /*  boolean, use the avg criterion technique or not*/
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
};

struct ABIP_SCALING 
{
      abip_float *D; 
      abip_float *E;
      
      abip_float mean_norm_row_A; 
      abip_float mean_norm_col_A;
};

/* main functions: ABIP(init): allocates memory of matrix, e.g., [I A; A^T -I]. 
                               ABIP(solve): can be called many times with different b,c data per init call. 
                               ABIP(finish): cleans up the memory (one per init call) */
ABIPWork *ABIP(init)(const ABIPData *d, ABIPInfo *info);
abip_int ABIP(solve)(ABIPWork *w, const ABIPData *d, ABIPSolution *sol, ABIPInfo *info);
void ABIP(finish)(ABIPWork *w);

abip_int ABIP(main)(const ABIPData *d, ABIPSolution *sol, ABIPInfo *info);
const char *ABIP(version)(void);

struct ABIP_WORK 
{
      abip_float sigma;                    
      abip_float gamma; 
      abip_int final_check; 
      abip_int double_check; 

      abip_float mu;                        
      abip_float beta;
      
      abip_float *u; 
      abip_float *v; 
      abip_float *u_t;  
      abip_float *u_prev; 
      abip_float *v_prev;

      abip_float* u_avg;
      abip_float* v_avg;

      abip_float* u_avgcon;
      abip_float* v_avgcon;

      abip_float* u_sumcon;
      abip_float* v_sumcon;


      abip_int    fre_old;

      abip_float *h; 
      abip_float *g; 
      abip_float *pr;  
      abip_float *dr;
  
      abip_float g_th; 
      abip_float sc_b; 
      abip_float sc_c; 
      abip_float nm_b; 
      abip_float nm_c;

      abip_float *b; 
      abip_float *c;
      abip_int m; 
      abip_int n;
      ABIPMatrix *A;
      abip_float sp; 

      ABIPLinSysWork *p;                   
      ABIPAdaptWork *adapt;            
      ABIPSettings *stgs;                    
      ABIPScaling *scal;                      
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
      
      abip_float tau;
      abip_float kap;
};

#ifdef __cplusplus
}
#endif
#endif
