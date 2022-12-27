#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "abip.h"
#include "glbopts.h"
#include "ctrlc.h"
#include "linalg.h"
#include "linsys.h"
#include "util.h"
#include "mkl.h"
#include "mkl_lapacke.h"
#include "cs.h"
#include "cones.h"
#include "lasso_config.h"
#include "svm_config.h"
#include "qcp_config.h"
#include "svm_qp_config.h"


ABIP(timer) global_timer;

/* printing header */
static const char* HEADER[] = {
      " ipm iter ", " admm iter ", "     mu ",
      " pri res ", " dua res ", " rel gap ",
      " pri obj ", " dua obj ", " kap/tau ", " time (s)",
};

static const abip_int HSPACE = 9;
static const abip_int HEADER_LEN = 10;
static const abip_int LINE_LEN = 150;

static abip_int abip_isnan
(
    abip_float x
) 
{
    DEBUG_FUNC
    RETURN(x == NAN || x != x);
}

static void free_work
(
    ABIPWork *w
) 
{
    DEBUG_FUNC
    
    if (!w) 
    {
        RETURN;
    }
    
    if (w->u) 
    {
        abip_free(w->u);
    }
    
    if (w->v) 
    {
        abip_free(w->v);
    }
    
    if (w->u_t) 
    {
        abip_free(w->u_t);
    }

    if (w->rel_ut) 
    {
        abip_free(w->rel_ut);
    }

    if (w->A)
    {
        ABIP(free_A_matrix)(w->A);
    }
    
    abip_free(w);
    
    RETURN;
}

static void print_init_header
(
    spe_problem *spe
) 
{
      DEBUG_FUNC
      
      abip_int i;
      ABIPSettings *stgs = spe->stgs;
      char *lin_sys_method = ABIP(get_lin_sys_method)(spe);
      
      
      for (i = 0; i < LINE_LEN; ++i) 
      {
            abip_printf("-");
      }
      
      abip_printf("\n\tABIP v%s - First-Order Interior-Point Solver for Conic Programming\n\t(c)  Jinsong Liu\n", 
            ABIP(version)());
      
      for (i = 0; i < LINE_LEN; ++i) 
      {
            abip_printf("-");
      }
      
      abip_printf("\n");
      
     if (lin_sys_method) 
      {
            abip_printf("Lin-sys: %s\n", lin_sys_method);
            abip_free(lin_sys_method);
      }

      if (stgs->normalize) 
      {
            abip_printf("eps_p = %.2e, eps_d = %.2e, eps_g = %.2e, alpha = %.2f, max_ipm_iters = %i, max_admm_iters = %i, normalize = %i\n"
                  "rho_y = %.2e\n", 
                  stgs->eps_p, stgs->eps_d, stgs->eps_g, stgs->alpha, (int)stgs->max_ipm_iters, (int)stgs->max_admm_iters, 
                  (int)stgs->normalize, stgs->rho_y);
      } 
      else 
      {
            abip_printf("eps_p = %.2e, eps_d = %.2e, eps_g = %.2e, alpha = %.2f, max_ipm_iters = %i, max_admm_iters = %i, normalize = %i\n"
                  "rho_y = %.2e\n", 
                  stgs->eps_p, stgs->eps_d, stgs->eps_g, stgs->alpha, (int)stgs->max_ipm_iters, (int)stgs->max_admm_iters, 
                  (int)stgs->normalize, stgs->rho_y);
      }

      abip_printf("constraints m = %i, Variables n = %i\n", (int)spe->p, (int)spe->q);
      

      #ifdef MATLAB_MEX_FILE

      mexEvalString("drawnow;");

      #endif

      RETURN;
}

static void populate_on_failure
(
    abip_int m,
    abip_int n,
    ABIPSolution* sol,
    ABIPInfo* info,
    abip_int status_val,
    const char* msg
)
{
    DEBUG_FUNC

        if (info)
        {
            info->res_pri = NAN;
            info->res_dual = NAN;
            info->rel_gap = NAN;
            info->res_infeas = NAN;
            info->res_unbdd = NAN;

            info->pobj = NAN;
            info->dobj = NAN;

            info->ipm_iter = -1;
            info->admm_iter = -1;
            info->status_val = status_val;
            info->solve_time = NAN;
            strcpy(info->status, msg);
        }

    if (sol)
    {
        if (n > 0)
        {
            if (!sol->x)
            {
                sol->x = (abip_float*)abip_malloc(sizeof(abip_float) * n);
            }
            ABIP(scale_array)(sol->x, NAN, n);

            if (!sol->s)
            {
                sol->s = (abip_float*)abip_malloc(sizeof(abip_float) * n);
            }
            ABIP(scale_array)(sol->s, NAN, m);
        }

        if (m > 0)
        {
            if (!sol->y)
            {
                sol->y = (abip_float*)abip_malloc(sizeof(abip_float) * m);
            }
            ABIP(scale_array)(sol->y, NAN, m);
        }
    }

    RETURN;
}

static abip_int failure
(
    ABIPWork* w,
    abip_int m,
    abip_int n,
    ABIPSolution* sol,
    ABIPInfo* info,
    abip_int stint,
    const char* msg,
    const char* ststr
)
{
    DEBUG_FUNC

        abip_int status = stint;
    populate_on_failure(m, n, sol, info, status, ststr);

    abip_printf("Failure:%s\n", msg);
    abip_end_interrupt_listener();

    RETURN status;
}



static abip_int projection
(
    ABIPWork *w, 
    abip_int iter,
    spe_problem *spe,
    ABIPResiduals *r
) 
{
    DEBUG_FUNC
    
    abip_int status;

    if(spe->Q != ABIP_NULL || spe->stgs->linsys_solver != 3){
        abip_int n = w->n; 
        abip_int m = w->m; 
        abip_int l = n + m + 1; 
        abip_float a = w->a;

        abip_float *mu = (abip_float*)abip_malloc((m+n)*sizeof(abip_float));
        memcpy(mu, w->u, (m+n)*sizeof(abip_float));
        ABIP(add_scaled_array)(mu, w->v, m+n, 1);
        ABIP(c_dot)(mu, spe->rho_dr, m+n);
        abip_float eta = spe->rho_dr[m+n] * (w->u[m+n] + w->v[m+n]);

        abip_float *warm_start = (abip_float*)abip_malloc((m+n)*sizeof(abip_float));
        memcpy(warm_start, w->u, (m+n)*sizeof(abip_float));

        abip_float *tem = (abip_float*)abip_malloc((m+n)*sizeof(abip_float));


        abip_float pcg_tol = 0;
        if(spe->stgs->linsys_solver == 3){// create warm start for pcg
            ABIP(add_scaled_array)(warm_start, w->r, m+n, w->u[m+n]);


            pcg_tol = MIN(r->Ax_b_norm, r->Qx_ATy_c_s_norm);
            pcg_tol = 0.2 * MIN(pcg_tol, ABIP(norm_inf)(warm_start, n) / POWF((abip_float)iter + 1, 1.5));

            pcg_tol = MAX(pcg_tol, 1e-12);
        }

        abip_float *p = (abip_float*)abip_malloc((m+n)*sizeof(abip_float));
        memcpy(p, mu, (m+n)*sizeof(abip_float));

        status = spe->solve_spe_linsys(spe, p, warm_start, iter, pcg_tol);

        memcpy(tem, p, (m+n)*sizeof(abip_float));
        ABIP(c_dot)(tem, spe->rho_dr, m+n);

        abip_float b = ABIP(dot)(w->r, mu, m+n) - 2*ABIP(dot)(w->r, tem, m+n) - eta;

        abip_float *Qp = (abip_float*)abip_malloc(n*sizeof(abip_float));
        memset(Qp, 0, n*sizeof(abip_float));

        if(spe->Q != ABIP_NULL){

            ABIP(accum_by_A)(spe->Q, &p[m], Qp);

        }

        abip_float c = -ABIP(dot)(&p[m], Qp, n);

        if(iter > 0){
            w->u_t[m+n] = (-b + SQRTF(MAX(0, b*b - 4*a*c))) / (2*a);
        }
        else{
            w->u_t[m+n] = 1;
        }

        memcpy(w->u_t, p, (m+n)*sizeof(abip_float));
        ABIP(add_scaled_array)(w->u_t, w->r, m+n, -w->u_t[m+n]);

        abip_free(mu);
        abip_free(warm_start);
        abip_free(p);
        abip_free(tem);
        abip_free(Qp);
    }
    else{
        abip_int n = w->n; 
        abip_int m = w->m; 
        abip_int l = n + m + 1; 
        abip_float a = w->a;

        abip_float *mu = (abip_float*)abip_malloc((m+n)*sizeof(abip_float));
        memcpy(mu, w->u, (m+n)*sizeof(abip_float));
        ABIP(add_scaled_array)(mu, w->v, m+n, 1);
        ABIP(c_dot)(mu, spe->rho_dr, m+n);

        abip_float eta = spe->rho_dr[m+n] * (w->u[m+n] + w->v[m+n]);

        abip_float *w3 = (abip_float*)abip_malloc((m+n)*sizeof(abip_float));
        memcpy(w3, mu, (m+n)*sizeof(abip_float));
        ABIP(add_scaled_array)(w3, spe->b, m, eta);
        ABIP(add_scaled_array)(&w3[m], spe->c, n, -eta);

        abip_float *warm_start = (abip_float*)abip_malloc((m+n)*sizeof(abip_float));
        memcpy(warm_start, w->u, (m+n)*sizeof(abip_float));

        abip_float pcg_tol = 0;
        if(spe->stgs->linsys_solver == 3){// create warm start for pcg

            ABIP(add_scaled_array)(warm_start, w->r, m+n, w->u[m+n]);

            pcg_tol = MIN(r->Ax_b_norm, r->Qx_ATy_c_s_norm);
            pcg_tol = 0.2 * MIN(pcg_tol, ABIP(norm_inf)(warm_start, n) / POWF((abip_float)iter + 1, 1.5));

            pcg_tol = MAX(pcg_tol, 1e-12);
        }

        memset(warm_start, 0, (m+n)*sizeof(abip_float));

        status = spe->solve_spe_linsys(spe, w3, warm_start, iter, r->error_ratio);

        abip_float coef = -(-ABIP(dot)(spe->b, w3, m) + ABIP(dot)(spe->c, &w3[m], n)) / (1 -ABIP(dot)(spe->b, w->r, m) + ABIP(dot)(spe->c, &w->r[m], n));

        memcpy(w->u_t, w3, (m+n)*sizeof(abip_float));
        ABIP(add_scaled_array)(w->u_t, w->r, m+n, coef);
        if(iter > 0){
            w->u_t[m+n] = eta -ABIP(dot)(spe->b, w->u_t, m) + ABIP(dot)(spe->c, &w->u_t[m], n);
        }
        else{
            w->u_t[m+n] = 1;
        }

        abip_free(mu);
        abip_free(w3);
        abip_free(warm_start);
    }
    

    RETURN status;
}




static void update_dual_vars
(
      ABIPWork *w
) 
{
    DEBUG_FUNC
      
    abip_int l = w->m + w->n + 1;

    memcpy(w->v, w->u, l*sizeof(abip_float));

    ABIP(add_scaled_array)(w->v,w->rel_ut,l,-1);

    RETURN;
}




static void solve_barrier_subproblem
(
    ABIPWork *w,
    const ABIPCone *c,
    spe_problem *spe
)
{
    DEBUG_FUNC

    abip_int n = w->n;
    abip_int m = w->m;
    abip_int l = m + n + 1;
    abip_float lambda = w->mu / w->beta;

    // over relaxation: rel_ut = alpha * ut + (1 - alpha) * u
    memcpy(w->rel_ut, w->u_t, l * sizeof(abip_float));
    ABIP(scale_array)(w->rel_ut, spe->stgs->alpha, l);
    ABIP(add_scaled_array)(w->rel_ut, w->u, l, 1 - spe->stgs->alpha);

    abip_float* tmp = (abip_float*)abip_malloc(l * sizeof(abip_float));

    ABIP(add_scaled_array)(w->rel_ut, w->v, l, -1);
    memcpy(tmp, w->rel_ut, l * sizeof(abip_float));

    abip_float* y = (abip_float*)abip_malloc(m * sizeof(abip_float));
    memcpy(y, tmp, m * sizeof(abip_float));

    abip_float tau = (tmp[l - 1] + SQRTF(tmp[l - 1]*tmp[l - 1] + 4 * lambda/spe->rho_dr[l - 1])) / 2;

    memcpy(w->u, y, m * sizeof(abip_float));
    w->u[l - 1] = tau;

    abip_int count = 0;
    abip_int i;

    /*soc*/
    if (c->qsize && c->q) {
        for (i = 0; i < c->qsize; ++i) {
            if (c->q[i] == 0) {
                continue;
            }
            if (c->q[i] == 1) {
                ABIP(positive_orthant_barrier_subproblem)(&w->u[m + count], &tmp[m + count], lambda/spe->rho_dr[m + count], 1);
            }
            else {
                ABIP(soc_barrier_subproblem)(&w->u[m + count], &tmp[m + count], lambda/spe->rho_dr[m + count], c->q[i]);
            }
            count += c->q[i];
        }
    }

    /*rsoc*/
    if (c->rqsize && c->rq) {
        for (i = 0; i < c->rqsize; ++i) {
            if (c->rq[i] < 3) {
                continue;
            }
            else {
                ABIP(rsoc_barrier_subproblem)(&w->u[m + count], &tmp[m + count], lambda/spe->rho_dr[m + count], c->rq[i]);
            }
            count += c->rq[i];
        }
    }

    /*free cone*/
    if (c->f) {
        ABIP(free_barrier_subproblem)(&w->u[m + count], &tmp[m + count], lambda/spe->rho_dr[m + count], c->f);
        count += c->f;
    }

    /*zero cone*/
    if (c->z) {
        ABIP(zero_barrier_subproblem)(&w->u[m + count], &tmp[m + count], lambda/spe->rho_dr[m + count], c->z);
        count += c->z;
    }

    /*positive orthant*/
    if (c->l) {
        ABIP(positive_orthant_barrier_subproblem)(&w->u[m + count], &tmp[m + count], lambda/spe->rho_dr[m + count], c->l);
        count += c->l;
    }
    
    abip_free(tmp);
    abip_free(y);
}


static abip_int indeterminate
(
    ABIPWork* w,
    ABIPSolution* sol,
    ABIPInfo* info
)
{
    DEBUG_FUNC

        strcpy(info->status, "Indeterminate");

    ABIP(scale_array)(sol->x, NAN, w->n);
    ABIP(scale_array)(sol->y, NAN, w->m);
    ABIP(scale_array)(sol->s, NAN, w->n);

    RETURN ABIP_INDETERMINATE;
}

static abip_int solved
(
    ABIPWork* w,
    ABIPSolution* sol,
    ABIPInfo* info,
    abip_float tau
)
{
    DEBUG_FUNC

    ABIP(scale_array)(sol->x, SAFEDIV_POS(1.0, tau), w->n);
    ABIP(scale_array)(sol->y, SAFEDIV_POS(1.0, tau), w->m);
    ABIP(scale_array)(sol->s, SAFEDIV_POS(1.0, tau), w->n);

    if ((info->status_val == 0)||(info->status_val == 2))
    {
        strcpy(info->status, "Solved/Inaccurate");
        RETURN ABIP_SOLVED_INACCURATE;
    }

    strcpy(info->status, "Solved");

    RETURN ABIP_SOLVED;
}
static void sety
(
      ABIPWork *w, 
      ABIPSolution *sol
) 
{
      DEBUG_FUNC
      
      if (!sol->y) 
      {
            sol->y = (abip_float *) abip_malloc(sizeof(abip_float) * w->m);
      }
      
      memcpy(sol->y, w->u, w->m * sizeof(abip_float));
      
      RETURN;
}

static void setx
(
      ABIPWork *w, 
      ABIPSolution *sol
) 
{
      DEBUG_FUNC
      
      if (!sol->x) 
      {
            sol->x = (abip_float *) abip_malloc(sizeof(abip_float) * w->n);
      }
      
      memcpy(sol->x, &(w->u[w->m]), w->n * sizeof(abip_float));
      
      RETURN;
}

static void sets
(
      ABIPWork *w, 
      ABIPSolution *sol
) 
{
      DEBUG_FUNC
      
      if (!sol->s) 
      {
            sol->s = (abip_float *) abip_malloc(sizeof(abip_float) * w->n);
      }
      
      memcpy(sol->s, &(w->v[w->m]), w->n * sizeof(abip_float));
      
      RETURN;
}

static abip_int infeasible
(
    ABIPWork* w,
    ABIPSolution* sol,
    ABIPInfo* info,
    abip_float bt_y
)
{
    DEBUG_FUNC

    ABIP(scale_array)(sol->y, 1 / bt_y, w->m);
    ABIP(scale_array)(sol->s, 1 / bt_y, w->n);
    ABIP(scale_array)(sol->x, NAN, w->n);

    if (info->status_val == 0)
    {
        strcpy(info->status, "Infeasible/Inaccurate");
        RETURN ABIP_INFEASIBLE_INACCURATE;
    }

    strcpy(info->status, "Infeasible");
    RETURN ABIP_INFEASIBLE;
}

static abip_int unbounded
(
    ABIPWork* w,
    ABIPSolution* sol,
    ABIPInfo* info,
    abip_float ct_x
)
{
    DEBUG_FUNC

    ABIP(scale_array)(sol->x, -1 / ct_x, w->n);
    ABIP(scale_array)(sol->y, NAN, w->m);
    ABIP(scale_array)(sol->s, NAN, w->n);

    if (info->status_val == 0)
    {
        strcpy(info->status, "Unbounded/Inaccurate");
        RETURN ABIP_UNBOUNDED_INACCURATE;
    }

    strcpy(info->status, "Unbounded");
    RETURN ABIP_UNBOUNDED;
}

static abip_int is_solved_status
(
    abip_int status
)
{
    RETURN status == ABIP_SOLVED || status == ABIP_SOLVED_INACCURATE;
}

static abip_int is_infeasible_status
(
    abip_int status
)
{
    RETURN status == ABIP_INFEASIBLE || status == ABIP_INFEASIBLE_INACCURATE;
}

static abip_int is_unbounded_status
(
    abip_int status
)
{
    RETURN status == ABIP_UNBOUNDED || status == ABIP_UNBOUNDED_INACCURATE;
}

static void get_info
(
    ABIPWork* w,
    ABIPSolution* sol,
    ABIPInfo* info,
    ABIPResiduals* r,
    abip_int ipm_iter,
    abip_int admm_iter
)
{
    DEBUG_FUNC

    info->ipm_iter = ipm_iter + 1;
    info->admm_iter = admm_iter;

    info->res_infeas = r->res_infeas;
    info->res_unbdd = r->res_unbdd;

    if (is_solved_status(info->status_val))
    {
        info->rel_gap = r->rel_gap;
        info->res_pri = r->res_pri;
        info->res_dual = r->res_dual;
        info->pobj = r->pobj;
        info->dobj = r->dobj;
    }
    else if (is_unbounded_status(info->status_val))
    {
        info->rel_gap = NAN;
        info->res_pri = NAN;
        info->res_dual = NAN;
        info->pobj = -INFINITY;
        info->dobj = -INFINITY;
    }
    else if (is_infeasible_status(info->status_val))
    {
        info->rel_gap = NAN;
        info->res_pri = NAN;
        info->res_dual = NAN;
        info->pobj = INFINITY;
        info->dobj = INFINITY;
    }

    RETURN;
}

static void get_solution
(
    ABIPWork* w,
    spe_problem *spe,
    ABIPSolution* sol,
    ABIPInfo* info,
    ABIPResiduals* r,
    abip_int ipm_iter,
    abip_int admm_iter
)
{
    DEBUG_FUNC

    abip_int l = w->m + w->n + 1;

    setx(w, sol);
    sety(w, sol);
    sets(w, sol);

    if (info->status_val == ABIP_UNFINISHED)
    {
        info->status_val = solved(w, sol, info, r->tau);
    }
    else if (is_solved_status(info->status_val))
    {
        info->status_val = solved(w, sol, info, r->tau);
    }
    else if (is_infeasible_status(info->status_val))
    {
        info->status_val = infeasible(w, sol, info, r->dobj * r->tau);
    }
    else
    {
        info->status_val = unbounded(w, sol, info, r->pobj * r->tau);
    }

    if (spe->stgs->normalize)
    {
        spe->un_scaling_sol(spe, sol);
    }

    get_info(w, sol, info, r, ipm_iter, admm_iter);

    RETURN;
}


static void print_summary
(
    ABIPWork *w, 
    spe_problem *spe,
    abip_int i, 
    abip_int j,
    abip_int k,
    ABIPResiduals *r,
    ABIP(timer) *solve_timer
) 
{
    DEBUG_FUNC

    if(!spe->stgs->verbose){
        RETURN;
    }
    
    abip_printf("%*i|", (int)strlen(HEADER[0]), (int)i+1);
    abip_printf("%*i|", (int)strlen(HEADER[1]), (int)j+1);
    abip_printf("%*.2e|", (int)strlen(HEADER[2]), w->mu);

    abip_printf("%*.2e|", (int)HSPACE, r->res_pri/spe->stgs->eps_p);
    abip_printf("%*.2e|", (int)HSPACE, r->res_dual/spe->stgs->eps_d);
    abip_printf("%*.2e|", (int)HSPACE, r->rel_gap/spe->stgs->eps_g);
    abip_printf("%*.2e|", (int)HSPACE, r->pobj);
    abip_printf("%*.2e|", (int)HSPACE, r->dobj);
    abip_printf("%*.2e|", (int)HSPACE, r->tau);

    abip_printf("%*.2e ", (int)HSPACE, ABIP(tocq)(solve_timer) / 1e3);
    abip_printf("\n");

    #if EXTRA_VERBOSE > 0
    
    abip_printf("Norm u = %4f, ", ABIP(norm)(w->u, w->n + w->m + 2));
    abip_printf("Norm u_t = %4f, ", ABIP(norm)(w->u_t, w->n + w->m + 2));
    abip_printf("Norm v = %4f, ", ABIP(norm)(w->v, w->n + w->m + 2));
    abip_printf("tau = %4f, ", r->tau);
    abip_printf("kappa = %4f, ", r->kap);
    abip_printf("|u - u_t| = %1.2e, ", ABIP(norm_diff)(w->u, w->u_t, w->n + w->m + 2));
    abip_printf("res_infeas = %1.2e, ", r->res_infeas);
    abip_printf("res_unbdd = %1.2e\n", r->res_unbdd);

    #endif

    #ifdef MATLAB_MEX_FILE

    mexEvalString("drawnow;");

    #endif
    
    RETURN;
}

static void print_header
(
      ABIPWork *w
) 
{
      DEBUG_FUNC
      
      abip_int i;
      
      for (i = 0; i < LINE_LEN; ++i) 
      {
                  abip_printf("-");
      }
      abip_printf("\n");
      
      for (i = 0; i < HEADER_LEN - 1; ++i) 
      {
                  abip_printf("%s|", HEADER[i]);
      }
      abip_printf("%s\n", HEADER[HEADER_LEN - 1]);
      
      for (i = 0; i < LINE_LEN; ++i) 
      {
                  abip_printf("-");
      }
      abip_printf("\n");
      
      #ifdef MATLAB_MEX_FILE

      mexEvalString("drawnow;");

      #endif
      
      RETURN;
}

static void print_footer
(
    const ABIPData* d,
    spe_problem *spe,
    ABIPSolution* sol,
    ABIPWork* w,
    ABIPInfo* info,
    abip_int k
)
{
    DEBUG_FUNC

        abip_int i;

    char* lin_sys_str = ABIP(get_lin_sys_summary)(spe, info);

    for (i = 0; i < LINE_LEN; ++i)
    {
        abip_printf("-");
    }

    abip_printf("\n");

    abip_printf("Status: %s\n", info->status);

    if (info->ipm_iter + 1 == spe->stgs->max_ipm_iters)
    {
        abip_printf("Hit max_ipm_iters, solution may be inaccurate\n");
    }

    if (info->admm_iter + 1 >= spe->stgs->max_admm_iters)
    {
        abip_printf("Hit max_admm_iters, solution may be inaccurate\n");
    }

    abip_printf("Timing: Solve time: %1.2es\n", info->solve_time / 1e3);
    abip_printf("          per iter: %1.2es\n", info->solve_time / (k*1e3));
    abip_printf("        Total time: %1.2es\n", (info->solve_time + info->setup_time)/ 1e3);

    if (lin_sys_str)
    {
        abip_printf("%s", lin_sys_str);
        abip_free(lin_sys_str);
    }

    for (i = 0; i < LINE_LEN; ++i)
    {
        abip_printf("-");
    }

    abip_printf("\n");

    if (is_infeasible_status(info->status_val))
    {
        abip_printf("Certificate of primal infeasibility:\n");
        abip_printf("|A'y + s|_2 * |b|_2 = %.4e\n", info->res_infeas);
        abip_printf("b'y = %.4f\n", ABIP(dot)(d->b, sol->y, d->m));
    }
    else if (is_unbounded_status(info->status_val))
    {
        abip_printf("Certificate of dual infeasibility:\n");
        abip_printf("|Ax|_2 * |c|_2 = %.4e\n", info->res_unbdd);
        abip_printf("c'x = %.4f\n", ABIP(dot)(d->c, sol->x, d->n));
    }
    else
    {
        abip_printf("Error metrics:\n");
        abip_printf("primal res: |Ax - b|_inf / (1 + max(|Ax|_inf, |b|_inf)) = %.4e\n", info->res_pri);
        abip_printf("dual res: |Qx - A'y + c - s|_inf / (1 + max(|Qx|_inf + |c|_inf)) = %.4e\n", info->res_dual);
        abip_printf("rel gap: |x'Qx + c'x - b'y| / (1 + max(|x'Qx| + |c'x| + |b'y|)) = %.4e\n", info->rel_gap);

        for (i = 0; i < LINE_LEN; ++i)
        {
            abip_printf("-");
        }

        abip_printf("\n");
        abip_printf("1/2x'Qx + c'x = %.4e, -1/2x'Qx + b'y = %.4e\n", info->pobj, info->dobj);
    }
    
    for (i = 0; i < LINE_LEN; ++i)
    {
        abip_printf("=");
    }

    abip_printf("\n");

#ifdef MATLAB_MEX_FILE

    mexEvalString("drawnow;");

#endif

    RETURN;
}


static abip_int has_converged
(
    ABIPWork* w,
    spe_problem *spe,
    ABIPResiduals* r,
    abip_int ipm_iter,
    abip_int admm_iter
)
{
    DEBUG_FUNC

    abip_float eps_p = spe->stgs->eps_p;
    abip_float eps_d = spe->stgs->eps_d;
    abip_float eps_g = spe->stgs->eps_g;
    abip_float eps_inf = spe->stgs->eps_inf;
    abip_float eps_unb = spe->stgs->eps_unb;


    if (r->res_pri < eps_p && r->res_dual < eps_d && r->rel_gap < eps_g)
    {
        RETURN ABIP_SOLVED;
    }

    if (r->res_dif <spe->stgs->err_dif * MAX(MAX(eps_p,eps_d), eps_g))
    {
        RETURN ABIP_SOLVED_INACCURATE;
    }

    if (r->res_unbdd < eps_unb && ipm_iter > 0 && admm_iter > 0)
    {
        RETURN ABIP_UNBOUNDED;
    }

    if (r->res_infeas < eps_inf && ipm_iter > 0 && admm_iter > 0)
    {
        RETURN ABIP_INFEASIBLE;
    }

    RETURN 0;
}

static abip_int validate
(
      const ABIPData *d,
      const ABIPCone *k,
      spe_problem *spe
) 
{
      DEBUG_FUNC
      
      ABIPSettings *stgs = d->stgs;
      
      if (d->n <= 0) 
      {
                  abip_printf("n must be greater than 0; n = %li\n", (long) d->n);
                  RETURN - 1;
      }
      
      if (spe->p > spe->q) 
      {
                  abip_printf("WARN: m larger than n, problem likely degenerate\n");
                  RETURN - 1;
      }
      
      if (ABIP(validate_lin_sys)(d->A) < 0) 
      {
                  abip_printf("invalid linear system input data\n");
                  RETURN - 1;
      }

      if (ABIP(validate_cones)(spe, k) < 0) {
          abip_printf("cone validation error\n");
          return -1;
      }
      
      if (stgs->max_ipm_iters <= 0) 
      {
                  abip_printf("max_ipm_iters must be positive\n");
                  RETURN - 1;
      }

      if (stgs->max_admm_iters <= 0) 
      {
                  abip_printf("max_admm_iters must be positive\n");
                  RETURN - 1;
      }

      if (stgs->eps_p <= 0 || stgs->eps_d <= 0 || stgs->eps_g <= 0 || stgs->eps_inf <= 0 || stgs->eps_unb <= 0) 
      {
                  abip_printf("eps tolerance must be positive\n");
                  RETURN - 1;
      }
      
      if (stgs->alpha <= 0 || stgs->alpha >= 2) 
      {
                  abip_printf("alpha must be in (0,2)\n");
                  RETURN - 1;
      }
      
      if (stgs->rho_y <= 0) 
      {
                  abip_printf("rho_y must be positive (1e-3 works well).\n");
                  RETURN - 1;
      }
      
      RETURN 0;
}

static ABIPWork *init_work
(     
    spe_problem *s, 
    ABIPCone *k
) 
{
    DEBUG_FUNC
    
    ABIPWork *w = (ABIPWork *) abip_calloc(1, sizeof(ABIPWork));
    abip_int l = s->p + s->q + 1;

    if (s->stgs->verbose) 
    {
        print_init_header(s);
    }

    if (!w) 
    {
        abip_printf("ERROR: allocating work failure\n");
        RETURN ABIP_NULL;
    }
    
    w->sigma = SIGMA;
    w->gamma = GAMMA;

    w->mu = 1.0;
    w->beta = 1.0;

    w->m = s->p;
    w->n = s->q;

    w->u = (abip_float *) abip_malloc(l * sizeof(abip_float));
    w->v = (abip_float *) abip_malloc(l * sizeof(abip_float));
    memset(w->u, 0, l * sizeof(abip_float));
    memset(w->v, 0, l * sizeof(abip_float));

    w->v_origin = (abip_float *) abip_malloc(l * sizeof(abip_float));
    w->u_t = (abip_float *) abip_malloc(l * sizeof(abip_float));
    w->rel_ut = (abip_float *) abip_malloc(l * sizeof(abip_float));
    w->r = (abip_float *)abip_malloc((w->n + w->m) * sizeof(abip_float));


    if (!w->u || !w->v || !w->u_t || !w->rel_ut || !w->v_origin || !w->r) 
    {
                abip_printf("ERROR: work memory allocation failure\n");
                RETURN ABIP_NULL;
    }
    
    w->nm_inf_b = ABIP(norm_inf)(s->data->b, s->p);
    w->nm_inf_c = ABIP(norm_inf)(s->data->c, s->q);
                  
    s->scaling_data(s, k);       
     
    if (s->init_spe_linsys_work(s))
    {
                abip_printf("ERROR: init_lin_sys_work failure\n");
                RETURN ABIP_NULL;
    }

    RETURN w;
}



static abip_int pre_calculate
(
    ABIPWork *w,
    spe_problem *spe   
)
{
    DEBUG_FUNC

    abip_int m = w->m;
    abip_int n = w->n;

    abip_float *zeros = (abip_float*)abip_malloc((m+n)*sizeof(abip_float));
    memset(zeros,0,(m+n)*sizeof(abip_float));

    memcpy(w->r, spe->b, m*sizeof(abip_float));
    ABIP(scale_array)(w->r, -1, m);
    memcpy(&w->r[m], spe->c, n*sizeof(abip_float));

    spe->solve_spe_linsys(spe, w->r, ABIP_NULL, -1, 1e-12);

    abip_float *tem = (abip_float*)abip_malloc((m+n)*sizeof(abip_float));
    memcpy(tem, w->r, (m+n)*sizeof(abip_float));
    ABIP(c_dot)(tem, spe->rho_dr, m+n);
    w->a = spe->rho_dr[m+n] + ABIP(dot)(tem, w->r, m+n);

    abip_free(tem);
    abip_free(zeros);

    RETURN 0;
}


abip_int update_work
(
    const ABIPData *d, 
    ABIPWork *w,
    const ABIPSolution *sol,
    const ABIPCone *c,
    spe_problem *spe
     
) 
{
    DEBUG_FUNC
    
    abip_int n = w->n;
    abip_int m = w->m; 
        
    //initialize x,y
    abip_float *x = (abip_float*)abip_malloc(n*sizeof(abip_float));
    abip_float *y = (abip_float*)abip_malloc(m*sizeof(abip_float));
    memset(y,0,m*sizeof(abip_float));
    /*-----initialize x----------*/
    abip_int i;
    abip_int count = 0;

    //soc
    if (c->qsize && c->q) {
        for (i = 0; i < c->qsize; ++i) {
            if (c->q[i] == 0) {
                continue;
            }
            memset(&x[count], 0, c->q[i] * sizeof(abip_float));
            x[count] = 1;
            count += c->q[i];
        }
    }

    //rsoc
    if (c->rqsize && c->rq) {
        for (i = 0; i < c->rqsize; ++i) {
            if (c->rq[i] < 3) {
                continue;
            }
            memset(&x[count], 0, c->rq[i] * sizeof(abip_float));
            x[count] = 1;
            x[count + 1] = 1;
            count += c->rq[i];
        }
    }

    //free cone
    if (c->f) {
        for (i = count; i < count+c->f; i++) {
            x[i] = 0;
        }
        count += c->f;
    }

    //zero cone
    if (c->z) {
        for (i = count; i < count+c->z; i++) {
            x[i] = 0;
        }
        count += c->z;
    }

    //positive orthant
    if (c->l) {
        for (i = count; i < count+c->l; i++) {
            x[i] = 1;
        }
        count += c->l;
    }

    //initialize u,v
    memcpy(w->u,y,m*sizeof(abip_float));
    memcpy(&w->u[m],x,n*sizeof(abip_float));
    w->u[m+n] = 1.0;
    // w->u[m+n+1] = 1;
    memcpy(w->v,y,m*sizeof(abip_float));
    memcpy(&w->v[m],x,n*sizeof(abip_float));
    w->v[m+n] = 1.0;

    pre_calculate(w, spe);
    
    abip_free(x);
    abip_free(y);
    return 0;
}


abip_float adjust_barrier(ABIPWork *w, ABIPResiduals *r, spe_problem *spe){
    
    abip_float error_ratio = r->error_ratio;



    abip_float sigma = 0.8;
    abip_float ratio = w->mu / MIN(MIN(spe->stgs->eps_p, spe->stgs->eps_d), spe->stgs->eps_g);
    abip_float gamma;

    if(ratio>50 && ratio<=100){
        gamma=1.5;
    }
    else if(ratio>10 && ratio<=50){
        gamma=1.3;
    }
    else if(ratio>5 && ratio<=10){
        gamma=1.2;
    }
    else if(ratio>1 && ratio<=5){
        gamma=1.1;
    }
    else if(ratio>0.5 && ratio<=1){
        gamma=1;
    }
    else if(ratio>0.1 && ratio<=0.5){
        gamma=0.9;
    }
    else if(ratio>0.05 && ratio<=0.1){
        gamma=0.9;
    }
    else if(ratio>0.01 && ratio<=0.05){
        gamma=0.8;
    }
    else if(ratio>0.005 && ratio<=0.01){
        gamma=0.8;
    }
    else if(ratio>0.001 && ratio<=0.005){
        gamma=0.7;
    }
    else if(ratio>0.0005 && ratio<=0.001){
        gamma=0.7;
    }
    else if(ratio>0.0001 && ratio<=0.0005){
        gamma=0.6;
    }
    else if(ratio>0.00005 && ratio<=0.0001){
        gamma=0.6;
    }
    else{
        gamma=0.5;        
    } 
    
    abip_float mix_ratio = error_ratio;

    if(mix_ratio>22){
        gamma=gamma*4.4;
    }
    else if(mix_ratio>18 && mix_ratio<=22){
        gamma=gamma*4.2;
    }
    else if(mix_ratio>15 && mix_ratio<=18){
        gamma=gamma*4;
    }
    else if(mix_ratio>12 && mix_ratio<=15){
        gamma=gamma*3.8;
    }
    else if(mix_ratio>8 && mix_ratio<=12){
        gamma=gamma*3.6;
    }
    else if(mix_ratio>6 && mix_ratio<=8){
        sigma=0.81;
        gamma=gamma*3.4;
    }
    else if(mix_ratio>4 && mix_ratio<=6){
        sigma=0.82;
        gamma=gamma*3.4;
    }
    else if(mix_ratio>3 && mix_ratio<=4){
        sigma=0.83;
        gamma=gamma*3.2;
    }
    else if(mix_ratio>3 && mix_ratio<=4){
        sigma=0.84;
        gamma=gamma*3;
    }
    else if(mix_ratio>2 && mix_ratio<=3){
        sigma=0.85;
        gamma=gamma*2.8;
    }
    else if(mix_ratio>1.5 && mix_ratio<=2){
        sigma=0.85;
        gamma=gamma*2.6;
    }
    else if(mix_ratio<1.5){
        sigma=0.85;
        gamma=gamma*2.4;
    }

    sigma = sigma * 0.2;
    
    w->mu = sigma * w->mu;
    return gamma * POWF(w->mu, spe->stgs->psi);
}


abip_int ABIP(solve)
(
    ABIPWork *w, 
    const ABIPData *d,
    ABIPSolution *sol, 
    ABIPInfo *info,
    ABIPCone *c,
    spe_problem *s
) 
{
    DEBUG_FUNC
    
    abip_int i;
    abip_int j; 
    abip_int k; 
    ABIP(timer) solve_timer;
    ABIP(timer) lin_timer;
    ABIP(timer) barrier_timer;
    ABIP(timer) res_timer;
    ABIP(timer) P_timer;
    ABIP(timer) uw_timer;
    abip_float lin_time=0;
    abip_float barrier_time=0;
    abip_float res_time=0;
    abip_float P_time=0;
    abip_float uw_time=0;

    abip_float time_limit_left = 1e3*s->stgs->time_limit - info->setup_time;
    
    ABIPResiduals *r = (ABIPResiduals*)abip_malloc(sizeof(ABIPResiduals));

    abip_int l = w->m + w->n + 1;
    
    if (!d || !sol || !info || !w ) 
    {
                abip_printf("ERROR: ABIP_NULL input\n");
                RETURN ABIP_FAILED;
    }

    abip_start_interrupt_listener();
    ABIP(tic)(&solve_timer);
    
    info->status_val = ABIP_UNFINISHED; 
    r->last_ipm_iter = -1;
    r->last_admm_iter = -1; 
    r->res_pri = 1e8;
    r->res_dual = 1e8;
    r->rel_gap = 1e8;
    r->error_ratio = 1e8;

    abip_float tol_inner = 4 * POWF(w->mu, s->stgs->psi);

    ABIP(tic)(&uw_timer);
    update_work(d, w, sol, c, s);
    uw_time += ABIP(tocq)(&uw_timer) / 1e3;

    if (s->stgs->verbose) 
    {
                print_header(w);
    }

    k = 0;
    
    for (i = 0; i < s->stgs->max_ipm_iters; ++i)
    {
                for (j = 0; j < s->stgs->max_admm_iters; ++j) 
                {
                    ABIP(tic)(&lin_timer);
                    if (projection(w, k, s, r) < 0) 
                    {
                        RETURN failure(w, w->m, w->n, sol, info, ABIP_FAILED, "error in project_lin_sys", "Failure");
                    }
                    lin_time += ABIP(tocq)(&lin_timer) / 1e3;
                    ABIP(tic)(&barrier_timer);
                    solve_barrier_subproblem(w, c, s);
                    barrier_time += ABIP(tocq)(&barrier_timer) / 1e3;

                    update_dual_vars(w);

                    memcpy(w->v_origin, w->v, l*sizeof(abip_float));
                    ABIP(c_dot)(w->v_origin, s->rho_dr, l);

                    k += 1; 
                        
                    ABIP(tic)(&P_timer);

                    abip_float err_inner= s->inner_conv_check(s,w);

                    if (err_inner < tol_inner || ABIP(tocq)(&solve_timer) > time_limit_left)
                    {
                        P_time += ABIP(tocq)(&P_timer) / 1e3;
                                
                        break; 
                    }
                    P_time += ABIP(tocq)(&P_timer) / 1e3;

                    if (abip_is_interrupted()) 
                    {
                                RETURN failure(w, w->m, w->n, sol, info, ABIP_SIGINT, "Interrupted", "Interrupted");
                    }
                    
                    #if EXTRA_VERBOSE > 0
                            abip_printf("primal error: %.4f, dual error: %.4f, gap: %.4f\n",r.err_pri/w->eps_p,r.err_dual/w->eps_d,r.gap/w->eps_g);
                    #endif

                    if((j+1) % s->stgs->inner_check_period == 0 || r->error_ratio <= 8){
                        ABIP(tic)(&res_timer);
                        s->calc_residuals(s,w,r,i,k);
                        res_time += ABIP(tocq)(&res_timer) / 1e3;

                        if((j+1) % s->stgs->inner_check_period == 0){
                            print_summary(w, s, i, j, k, r, &solve_timer);
                        }

                        if ((info->status_val = has_converged(w, s, r, i ,k)) != 0 || k+1 >= s->stgs->max_admm_iters*s->stgs->max_ipm_iters || i+1 >= s->stgs->max_ipm_iters||ABIP(tocq)(&solve_timer) > time_limit_left) //max running time is time_limit s
                        {
                            if (s->stgs->verbose && k>0) 
                            {
                            printf("\nin last admm iter:\n");
                            print_summary(w, s, i, j, k, r, &solve_timer);
                            abip_printf("total admm iter is %i\n", k);
                            }
                            
                            get_solution(w, s, sol, info, r, i, k);
                            info->solve_time = ABIP(tocq)(&solve_timer);

                            if (s->stgs->verbose) 
                            {
                                print_footer(d, s, sol, w, info, k);
                                printf("\ntotal time of project_lin_sys: %.2es\ntotal time of solve_barrier_subproblem: %.2es\ntotal time of calculate res: %.2es\ntotal time of calculate err_inner: %.2es\ntotal time of updating work: %.2es\n", lin_time, barrier_time, res_time, P_time, uw_time);

                            }

                            abip_end_interrupt_listener();

                            RETURN info->status_val;
                        }
                    }

                                                                
                }//inner for

                if(s->sparsity || (i+1) % s->stgs->outer_check_period == 0){
                    ABIP(tic)(&res_timer);
                    s->calc_residuals(s,w,r,i,k);
                    res_time += ABIP(tocq)(&res_timer) / 1e3;
                    print_summary(w, s, i, j, k, r, &solve_timer);
                    if ((info->status_val = has_converged(w, s, r, i ,k)) != 0 || k+1 >= s->stgs->max_admm_iters*s->stgs->max_ipm_iters || i+1 >= s->stgs->max_ipm_iters||ABIP(tocq)(&solve_timer) > time_limit_left) 
                    {
                        if (s->stgs->verbose && k>0) 
                        {
                        printf("\nin last admm iter:\n");
                        print_summary(w, s, i, j, k, r, &solve_timer);
                        abip_printf("total admm iter is %i\n", k);
                        }
                        
                        get_solution(w, s, sol, info, r, i, k);
                        info->solve_time = ABIP(tocq)(&solve_timer);

                        if (s->stgs->verbose) 
                        {
                            print_footer(d, s, sol, w, info, k);
                            printf("\ntotal time of project_lin_sys: %.2es\ntotal time of solve_barrier_subproblem: %.2es\ntotal time of calculate res: %.2es\ntotal time of calculate err_inner: %.2es\n", lin_time, barrier_time, res_time, P_time);

                        }

                        abip_end_interrupt_listener();

                        RETURN info->status_val;
                    }
                }

            tol_inner = adjust_barrier(w,r,s);

    }
    
    RETURN info->status_val;
}


void ABIP(finish)
(
    ABIPWork *w,
    spe_problem *spe
) 
{
    DEBUG_FUNC
    
    if (w) 
    {          
    free_work(w);
    }
    
    if(spe->L)
    {
        spe->free_spe_linsys_work(spe);
    }

      RETURN;
}

ABIPWork *ABIP(init)
(
    const ABIPData *d, 
    ABIPInfo *info,
    spe_problem *s,
    ABIPCone *c
) 
{
    DEBUG_FUNC
    
    #if EXTRA_VERBOSE > 1
    ABIP(tic)(&global_timer);
    #endif
    
    ABIPWork *w;
    ABIP(timer) init_timer;
    abip_start_interrupt_listener();
    
    if (!d || !info) 
    {
                abip_printf("ERROR: Missing ABIPData or ABIPInfo input\n");
                RETURN ABIP_NULL;
    }
    
    #if EXTRA_VERBOSE > 0
    ABIP(print_data)(d);
    #endif
    
    #ifndef NOVALIDATE
    if (validate(d,c,s) < 0) 
    {
                abip_printf("ERROR: Validation returned failure\n");
                RETURN ABIP_NULL;
    }
    #endif
    
    ABIP(tic)(&init_timer);
    
    w = init_work(s,c);
    info->setup_time = ABIP(tocq)(&init_timer);
    
    if (d->stgs->verbose) 
    {
                abip_printf("Setup time: %1.2es\n", info->setup_time / 1e3);
    }
    
    abip_end_interrupt_listener();
    
    RETURN w;
}

abip_int ABIP(init_problem)(
    spe_problem **s,
    ABIPData *d, 
    ABIPSettings *stgs,
    enum problem_type special_problem
)
{
    switch(special_problem){

        case LASSO:
            return init_lasso((lasso **)s,d,stgs);
        case SVM:
            return init_svm((svm **)s,d,stgs);
        case QCP:
            return init_qcp((qcp **)s,d,stgs);
        case SVMQP:
            return init_svmqp((svmqp **)s,d,stgs);
        default:
            return init_qcp((qcp **)s,d,stgs);
    }
}

abip_int abip
(
    const ABIPData* d,
    ABIPSolution* sol,
    ABIPInfo* info,
    ABIPCone *K
)
{
    DEBUG_FUNC

    spe_problem *s = (spe_problem *)abip_malloc(sizeof(spe_problem));
    enum problem_type prob_type;
    if(d->stgs->prob_type == 0)prob_type = LASSO;
    else if(d->stgs->prob_type == 1)prob_type = SVM;
    else if(d->stgs->prob_type == 2)prob_type = QCP;
    else if(d->stgs->prob_type == 3)prob_type = SVMQP;

    ABIP(init_problem)(&s,d,d->stgs,prob_type);
    abip_int status;
    ABIPWork* w = ABIP(init)(d, info, s, K);//call init_work() to call init_lin_sys_work() to perform LDL' factorization

#if EXTRA_VERBOSE > 0
    abip_printf("size of abip_int = %lu, size of abip_float = %lu\n", sizeof(abip_int), sizeof(abip_float));
#endif

    if (w)
    {
        ABIP(solve)(w, d, sol, info, K, s);
        status = info->status_val;
        ABIP(finish)(w,s);
    }
    else
    {
        status = failure(ABIP_NULL, d ? d->m : -1, d ? d->n : -1, sol, info, ABIP_FAILED, "could not initialize work", "Failure");
    }
    
    RETURN status;
}
