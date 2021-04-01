#include "abip.h"
#include "glbopts.h"
#include "adaptive.h"
#include "ctrlc.h"
#include "linalg.h"
#include "linsys.h"
#include "normalize.h"
#include "util.h"

ABIP(timer) global_timer;

/* printing header */
static const char *HEADER[] = {
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
      
      if (w->u_prev) 
      {
            abip_free(w->u_prev);
      }
      
      if (w->v_prev) 
      {
            abip_free(w->v_prev);
      }
      
      if (w->h) 
      {
            abip_free(w->h);
      }
      
      if (w->g) 
      {
            abip_free(w->g);
      }

      if (w->pr) 
      {
            abip_free(w->pr);
      }
      
      if (w->dr) 
      {
            abip_free(w->dr);
      }

      if (w->b) 
      {
            abip_free(w->b);
      }
      
      if (w->c) 
      {
            abip_free(w->c);
      }
  
      if (w->scal) 
      {
            if (w->scal->D) 
            {
                  abip_free(w->scal->D);
            }
            
            if (w->scal->E) 
            {
                  abip_free(w->scal->E);
            }
            
            abip_free(w->scal);
      }

      abip_free(w);
      
      RETURN;
}

static void print_init_header
(
      const ABIPData *d
) 
{
      DEBUG_FUNC
      
      abip_int i;
      ABIPSettings *stgs = d->stgs;
      char *lin_sys_method = ABIP(get_lin_sys_method)(d->A, d->stgs);
      
      abip_int adaptive_lookback = stgs->adaptive_lookback;

      for (i = 0; i < LINE_LEN; ++i) 
      {
            abip_printf("-");
      }
      
      abip_printf("\n\tABIP v%s - First-Order Interior-Point Solver\n\t(c) Tianyi Lin, UC Berkeley, 2017-2018\n", 
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
            abip_printf("eps = %.2e, alpha = %.2f, max_ipm_iters = %i, max_admm_iters = %i, normalize = %i\n"
                  "scale = %2.2f, adaptive = %i, adaptive_lookback = %i, rho_y = %.2e\n", 
                  stgs->eps, stgs->alpha, (int)stgs->max_ipm_iters, (int)stgs->max_admm_iters, 
                  (int)stgs->normalize, stgs->scale,  (int)stgs->adaptive, stgs->adaptive_lookback, stgs->rho_y);
      } 
      else 
      {
            abip_printf("eps = %.2e, alpha = %.2f, max_ipm_iters = %i, max_admm_iters = %i, normalize = %i\n"
                  "adaptive = %i, adaptive_lookback = %i, rho_y = %.2e\n", 
                  stgs->eps, stgs->alpha, (int)stgs->max_ipm_iters, (int)stgs->max_admm_iters, 
                  (int)stgs->normalize, (int)stgs->adaptive, stgs->adaptive_lookback, stgs->rho_y);
      }

      abip_printf("Variables n = %i, constraints m = %i\n", (int)d->n, (int)d->m);

      #ifdef MATLAB_MEX_FILE

      mexEvalString("drawnow;");

      #endif

      RETURN;
}

static void populate_on_failure
( 
      abip_int m, 
      abip_int n, 
      ABIPSolution *sol,
      ABIPInfo *info, 
      abip_int status_val,
      const char *msg
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
                          sol->x = (abip_float *)abip_malloc(sizeof(abip_float) * n);
                  }
                  ABIP(scale_array)(sol->x, NAN, n);

                  if (!sol->s) 
                  {
                          sol->s = (abip_float *)abip_malloc(sizeof(abip_float) * n);
                  }
                  ABIP(scale_array)(sol->s, NAN, m);
            }
            
            if (m > 0) 
            {
                  if (!sol->y) 
                  {
                          sol->y = (abip_float *)abip_malloc(sizeof(abip_float) * m);
                  }
                  ABIP(scale_array)(sol->y, NAN, m);
            }
      }

      RETURN;
}

static abip_int failure
(
      ABIPWork *w, 
      abip_int m, 
      abip_int n, 
      ABIPSolution *sol,
      ABIPInfo *info, 
      abip_int stint, 
      const char *msg,
      const char *ststr
) 
{
      DEBUG_FUNC
      
      abip_int status = stint;
      populate_on_failure(m, n, sol, info, status, ststr);
      
      abip_printf("Failure:%s\n", msg);
      abip_end_interrupt_listener();
      
      RETURN status;
}

static void warm_start_vars
(
      ABIPWork *w, 
      const ABIPSolution *sol
) 
{
      DEBUG_FUNC
      
      abip_int i; 
      abip_int n = w->n; 
      abip_int m = w->m;
      
      memset(w->v, 0, m * sizeof(abip_float));
      memcpy(w->u, sol->y, m * sizeof(abip_float));
      memcpy(&(w->u[m]), sol->x, n * sizeof(abip_float));
      memcpy(&(w->v[m]), sol->s, n * sizeof(abip_float));
      w->u[n + m] = 1.0;
      w->v[n + m] = 0.0;
      
      #ifndef NOVALIDATE
      
      for (i = 0; i < n + m + 1; ++i) 
      {
            if (abip_isnan(w->u[i]) && i<m) 
            {
                  w->u[i] = 0; 
            }
            else
            {
                  w->u[i] = SQRTF(w->mu/w->beta);
            }
                  
            if (abip_isnan(w->v[i])) 
            {
                  w->v[i] = 0;
            }
            else
            {
                  w->v[i] = SQRTF(w->mu/w->beta);
            }
      }
      
      #endif
      
      if (w->stgs->normalize) 
      {
            ABIP(normalize_warm_start)(w);
      }
      
      RETURN;
}
 
static void cold_start_vars
( 
      ABIPWork *w
) 
{
      DEBUG_FUNC
      
      abip_int l = w->m + w->n + 1;
      abip_int i;
      
      memset(w->u, 0, w->m * sizeof(abip_float));
      memset(w->v, 0, w->m * sizeof(abip_float));

      for (i = w->m; i < l; ++i) 
      {
            w->u[i] = SQRTF(w->mu/w->beta);
            w->v[i] = SQRTF(w->mu/w->beta);
      }
      
      RETURN;
}

static abip_float calc_primal_resid
(
      ABIPWork *w, 
      const abip_float *x, 
      const abip_float tau,
      abip_float *nm_A_x
) 
{
      DEBUG_FUNC
      
      abip_int i;
      
      abip_float pres = 0; 
      abip_float scale; 
      abip_float *pr = w->pr;
      
      *nm_A_x = 0;
      
      memset(pr, 0, w->m * sizeof(abip_float));
      ABIP(accum_by_A)(w->A, w->p, x, pr);             
      
      for (i = 0; i < w->m; ++i) 
      {
            scale = w->stgs->normalize ? w->scal->D[i] / (w->sc_b * w->stgs->scale) : 1;
            scale = scale * scale;
            *nm_A_x += (pr[i] * pr[i]) * scale;
            pres += (pr[i] - w->b[i] * tau) * (pr[i] - w->b[i] * tau) * scale;
      }
      
      *nm_A_x = SQRTF(*nm_A_x);                           
      RETURN SQRTF(pres);                                      
}

static abip_float calc_dual_resid
(
      ABIPWork *w, 
      const abip_float *y,
      const abip_float *s,
      const abip_float tau, 
      abip_float *nm_At_ys
) 
{
      DEBUG_FUNC
      
      abip_int i;
      abip_float dres = 0; 
      abip_float scale; 
      abip_float *dr = w->dr;
      
      *nm_At_ys = 0;
      
      memset(dr, 0, w->n * sizeof(abip_float));
      ABIP(accum_by_Atrans)(w->A, w->p, y, dr);      
      ABIP(add_scaled_array)(dr, s, w->n, 1.0);          

      for (i = 0; i < w->n; ++i) 
      {
            scale = w->stgs->normalize ? w->scal->E[i] / (w->sc_c * w->stgs->scale) : 1;
            scale = scale * scale;
            *nm_At_ys += (dr[i] * dr[i]) * scale;
            dres += (dr[i] - w->c[i] * tau) * (dr[i] - w->c[i] * tau) * scale;
      }
      
      *nm_At_ys = SQRTF(*nm_At_ys);                    
      RETURN SQRTF(dres);                                      
}

static void calc_residuals
(
      ABIPWork *w, 
      ABIPResiduals *r, 
      abip_int ipm_iter, 
      abip_int admm_iter
) 
{
      DEBUG_FUNC
      
      abip_float *y = w->u; 
      abip_float *x = &(w->u[w->m]); 
      abip_float *s = &(w->v[w->m]);

      abip_float nmpr_tau; 
      abip_float nmdr_tau; 
      abip_float nm_A_x_tau; 
      abip_float nm_At_ys_tau; 
      abip_float ct_x; 
      abip_float bt_y;
      
      abip_int n = w->n; 
      abip_int m = w->m;

      if (admm_iter && r->last_admm_iter == admm_iter) 
      {
            RETURN;
      }
      
      r->last_ipm_iter = ipm_iter;
      r->last_admm_iter = admm_iter; 

      r->tau = ABS(w->u[n + m]);
      r->kap = ABS(w->v[n + m]) / (w->stgs->normalize ? (w->stgs->scale * w->sc_c * w->sc_b) : 1);

      nmpr_tau = calc_primal_resid(w, x, r->tau, &nm_A_x_tau);
      nmdr_tau = calc_dual_resid(w, y, s, r->tau, &nm_At_ys_tau);

      r->bt_y_by_tau = ABIP(dot)(y, w->b, m) / (w->stgs->normalize ? (w->stgs->scale * w->sc_c * w->sc_b) : 1);
      r->ct_x_by_tau = ABIP(dot)(x, w->c, n) / (w->stgs->normalize ? (w->stgs->scale * w->sc_c * w->sc_b) : 1);

      r->res_infeas = r->bt_y_by_tau > 0 ? w->nm_b * nm_At_ys_tau / r->bt_y_by_tau : NAN;
      r->res_unbdd = r->ct_x_by_tau < 0 ? w->nm_c * nm_A_x_tau / -r->ct_x_by_tau : NAN;

      bt_y = SAFEDIV_POS(r->bt_y_by_tau, r->tau);
      ct_x = SAFEDIV_POS(r->ct_x_by_tau, r->tau);

      r->res_pri = SAFEDIV_POS(nmpr_tau / (1 + w->nm_b), r->tau);
      r->res_dual = SAFEDIV_POS(nmdr_tau / (1 + w->nm_c), r->tau);
      r->rel_gap = ABS(ct_x - bt_y) / (1 + ABS(ct_x) + ABS(bt_y));

      RETURN;
}

static abip_int project_lin_sys
(
      ABIPWork *w, 
      abip_int iter
) 
{
      DEBUG_FUNC
      
      abip_int n = w->n; 
      abip_int m = w->m; 
      abip_int l = n + m + 1; 

      abip_int status;
      memcpy(w->u_t, w->u, l * sizeof(abip_float));
      ABIP(add_scaled_array)(w->u_t, w->v, l, 1.0);

      ABIP(scale_array)(w->u_t, w->stgs->rho_y, m);
      ABIP(add_scaled_array)(w->u_t, w->h, l - 1, -w->u_t[l - 1]);
      ABIP(add_scaled_array)(w->u_t, w->h, l - 1, -ABIP(dot)(w->u_t, w->g, l - 1) / (w->g_th + 1));
      ABIP(scale_array)(&(w->u_t[m]), -1, n);
      status = ABIP(solve_lin_sys)(w->A, w->stgs, w->p, w->u_t, w->u, iter);
      w->u_t[l - 1] += ABIP(dot)(w->u_t, w->h, l - 1);
      RETURN status;
}

static void update_dual_vars
(
      ABIPWork *w
) 
{
      DEBUG_FUNC
      
      abip_int i; 
      abip_int m = w->m; 
      abip_int l = m + w->n + 1;
      
      for (i = m; i < l; ++i) 
      {
             w->v[i] += (w->u[i] - w->stgs->alpha * w->u_t[i] - (1.0 - w->stgs->alpha) * w->u_prev[i]);
      }
      
      RETURN;
}

static void project_barrier
(
      ABIPWork *w
) 
{
      DEBUG_FUNC
      
      abip_int i; 
      abip_int m = w->m; 
      abip_int l = m + w->n + 1; 
      abip_int status;

      abip_float tmp; 
      
      for (i = 0; i < m; ++i) 
      {
            w->u[i] = w->u_t[i] - w->v[i];
      }
      
      for (i = m; i < l; ++i) 
      {
            w->u[i] = w->stgs->alpha * w->u_t[i] + (1 - w->stgs->alpha) * w->u_prev[i] - w->v[i];
      }
      
      for(i = m; i < l; ++i)
      {
            tmp = w->u[i] / 2; 
            w->u[i] = tmp + SQRTF(tmp * tmp + w->mu / w->beta); 
      }

      RETURN;
}

static void update_barrier
(
      ABIPWork *w,
      ABIPResiduals *r
) 
{
      abip_float sigma = w->sigma; 
      abip_float mu = w->mu; 
      abip_float gamma = w->gamma; 

      abip_float ratio = w->mu / w->stgs->eps; 
      abip_float err_ratio = MAX(MAX(r->res_pri, r->res_dual), r->rel_gap) / w->stgs->eps; 
      
      if (MAX(w->sp,w->stgs->sparsity_ratio) > 0.4 || MIN(w->sp,w->stgs->sparsity_ratio) > 0.1)
      {
            if (ratio > 10.0)
            {
                  gamma = 2.0; 
            } 
            else if (ratio > 1.0 && ratio <= 10.0)
            {
                  gamma = 1.0; 
            }
            else if (ratio > 0.5 && ratio <= 1.0)
            {
                  gamma = 0.9; 
            }
            else if (ratio > 0.1 && ratio <= 0.5)
            {
                  gamma = 0.8; 
            }
            else if (ratio > 0.05 && ratio <= 0.1)
            {
                  gamma = 0.7;
            }
            else if (ratio > 0.01 && ratio <= 0.05)
            {
                  gamma = 0.6;
            }
            else if (ratio > 0.005 && ratio <= 0.01)
            {
                  gamma = 0.5; 
            }
            else if (ratio > 0.001 && ratio <= 0.005)
            {
                  gamma = 0.4; 
            }
            else
            {
                  gamma = 0.3; 
            }

            if (err_ratio > 6 && err_ratio <= 10)
            {
                  sigma = 0.5;
            }
            else if (err_ratio > 3 && err_ratio <= 6)
            {
                  sigma = 0.6; 
                  gamma = gamma * 0.8;
            }
            else if (err_ratio > 1 && err_ratio <= 3)
            {
                  w->final_check = 1;
                  gamma = gamma * 0.4; 
                  if (ratio < 0.1) 
                  {
                        sigma = 0.8;
                  }
                  else
                  {
                        sigma = 0.7; 
                  }
                       
            }
            else
            {
                  sigma = w->sigma; 
            }
      }
      else
      {
            if (ratio > 10.0)
            {
                  gamma = 3.0; 
            } 
            else if (ratio > 1.0 && ratio <= 10.0)
            {
                  gamma = 1.0; 
            }
            else if (ratio > 0.5 && ratio <= 1.0)
            {
                  gamma = 0.9; 
            }
            else if (ratio > 0.1 && ratio <= 0.5)
            {
                  gamma = 0.8; 
            }
            else if (ratio > 0.05 && ratio <= 0.1)
            {
                  gamma = 0.7;
            }
            else if (ratio > 0.01 && ratio <= 0.05)
            {
                  gamma = 0.6;
            }
            else if (ratio > 0.005 && ratio <= 0.01)
            {
                  gamma = 0.5; 
            }
            else if (ratio > 0.001 && ratio <= 0.005)
            {
                  gamma = 0.4; 
            }
            else
            {
                  gamma = 0.3; 
            }

            if (err_ratio > 6 && err_ratio <= 10)
            {
                  sigma = 0.82; 
                  gamma = gamma * 0.8; 
            }
            else if (err_ratio > 4 && err_ratio <= 6)
            {
                  sigma = 0.84; 
                  gamma = gamma * 0.6; 
            }
            else if (err_ratio > 3 && err_ratio <= 4)
            {
                  sigma = 0.85; 
                  gamma = gamma * 0.5; 
                  w->final_check = 1;      
            }
            else if (err_ratio > 1 && err_ratio <= 3)
            {
                  w->final_check = 1;
                  if (ratio < 0.1) 
                  {
                        if (w->double_check)
                        {
                              sigma = 0.9; 
                              gamma = gamma * 0.4; 
                              w->double_check = 0; 
                        }
                        else
                        {
                              sigma = 1.0; 
                              gamma = gamma * 0.1;
                              w->double_check = 1;
                        }
                  }
                  else
                  {
                        sigma = 0.88; 
                        gamma = gamma * 0.4; 
                  }
            }
            else
            {
                  sigma = w->sigma; 
            }
      }

      mu = mu * sigma; 

      w->mu = mu; 
      w->sigma = sigma; 
      w->gamma = gamma; 
}

static void reinitialize_vars
(
      ABIPWork *w
) 
{
      DEBUG_FUNC
      
      abip_int i; 
      abip_int m = w->m; 
      abip_int l = m + w->n + 1;
      
      for (i = m; i < l; ++i) 
      {
            if (w-> u[i] > w->v[i])
            {
                  w->v[i] = w->sigma * w->v[i];
            }
            else if (w-> u[i] < w->v[i]) 
            {
                  w->u[i] = w->sigma * w->u[i];
            }
            else
            {
                  w->u[i] = SQRTF(w->sigma) * w->u[i];
                  w->v[i] = SQRTF(w->sigma) * w->v[i];
            }      
      }      
      RETURN;
}

static abip_int indeterminate
(
      ABIPWork *w, 
      ABIPSolution *sol, 
      ABIPInfo *info
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
      ABIPWork *w, 
      ABIPSolution *sol, 
      ABIPInfo *info,
      abip_float tau
) 
{
      DEBUG_FUNC
      
      ABIP(scale_array)(sol->x, SAFEDIV_POS(1.0, tau), w->n);
      ABIP(scale_array)(sol->y, SAFEDIV_POS(1.0, tau), w->m);
      ABIP(scale_array)(sol->s, SAFEDIV_POS(1.0, tau), w->n);
      
      if (info->status_val == 0) 
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
      ABIPWork *w, 
      ABIPSolution *sol, 
      ABIPInfo *info,
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
      ABIPWork *w, 
      ABIPSolution *sol, 
      ABIPInfo *info,
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
      ABIPWork *w, 
      ABIPSolution *sol, 
      ABIPInfo *info,
      ABIPResiduals *r, 
      abip_int ipm_iter, 
      abip_int admm_iter
) 
{
      DEBUG_FUNC
      
      info->ipm_iter = ipm_iter + 1;
      info->admm_iter = admm_iter + 1; 

      info->res_infeas = r->res_infeas;
      info->res_unbdd = r->res_unbdd;
  
      if (is_solved_status(info->status_val)) 
      {
                  info->rel_gap = r->rel_gap;
                  info->res_pri = r->res_pri;
                  info->res_dual = r->res_dual;
                  info->pobj = r->ct_x_by_tau / r->tau;
                  info->dobj = r->bt_y_by_tau / r->tau;
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
      ABIPWork *w, 
      ABIPSolution *sol, 
      ABIPInfo *info,
      ABIPResiduals *r, 
      abip_int ipm_iter, 
      abip_int admm_iter
) 
{
      DEBUG_FUNC
      
      abip_int l = w->m + w->n + 1;
      
      calc_residuals(w, r, ipm_iter, admm_iter);
      setx(w, sol);
      sety(w, sol);
      sets(w, sol);
      
      if (info->status_val == ABIP_UNFINISHED) 
      {
                  if (r->tau > INDETERMINATE_TOL && r->tau > r->kap) 
                  {
                              info->status_val = solved(w, sol, info, r->tau);
                  } 
                  else if (ABIP(norm)(w->u, l) < INDETERMINATE_TOL * SQRTF((abip_float)l)) 
                  {
                              info->status_val = indeterminate(w, sol, info);
                  } 
                  else if (-r->bt_y_by_tau < r->ct_x_by_tau) 
                  {
                              info->status_val = infeasible(w, sol, info, r->bt_y_by_tau);
                  } 
                  else
                  {
                              info->status_val = unbounded(w, sol, info, r->ct_x_by_tau);
                  }
      } 
      else if (is_solved_status(info->status_val)) 
      {
                  info->status_val = solved(w, sol, info, r->tau);
      } 
      else if (is_infeasible_status(info->status_val)) 
      {
                  info->status_val = infeasible(w, sol, info, r->bt_y_by_tau);
      } 
      else
      {
                  info->status_val = unbounded(w, sol, info, r->ct_x_by_tau);
      }
      
      if (w->stgs->normalize) 
      {
                  ABIP(un_normalize_sol)(w, sol);
      }
      
      get_info(w, sol, info, r, ipm_iter, admm_iter);
      
      RETURN;
}

static void print_summary
(
      ABIPWork *w, 
      abip_int i, 
      abip_int j,
      ABIPResiduals *r,
      ABIP(timer) *solve_timer
) 
{
      DEBUG_FUNC
      
      abip_printf("%*i|", (int) strlen(HEADER[0]), (int) i);
      abip_printf("%*i|", (int) strlen(HEADER[1]), (int) j);
      abip_printf("%*.2e|", (int) strlen(HEADER[2]), w->mu);

      abip_printf("%*.2e|", (int) HSPACE, r->res_pri);
      abip_printf("%*.2e|", (int) HSPACE, r->res_dual);
      abip_printf("%*.2e|", (int) HSPACE, r->rel_gap);
      abip_printf("%*.2e|", (int) HSPACE, SAFEDIV_POS(r->ct_x_by_tau, r->tau));
      abip_printf("%*.2e|", (int) HSPACE, SAFEDIV_POS(r->bt_y_by_tau, r->tau));
      abip_printf("%*.2e|", (int) HSPACE, SAFEDIV_POS(r->kap, r->tau));
      abip_printf("%*.2e ", (int) HSPACE, ABIP(tocq)(solve_timer) / 1e3);
      abip_printf("\n");

      #if EXTRA_VERBOSE > 0
      
      abip_printf("Norm u = %4f, ", ABIP(norm)(w->u, w->n + w->m + 1));
      abip_printf("Norm u_t = %4f, ", ABIP(norm)(w->u_t, w->n + w->m + 1));
      abip_printf("Norm v = %4f, ", ABIP(norm)(w->v, w->n + w->m + 1));
      abip_printf("tau = %4f, ", r->tau);
      abip_printf("kappa = %4f, ", r->kap);
      abip_printf("|u - u_prev| = %1.2e, ", ABIP(norm_diff)(w->u, w->u_prev, w->n + w->m + 1));
      abip_printf("|u - u_t| = %1.2e, ", ABIP(norm_diff)(w->u, w->u_t, w->n + w->m + 1));
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
      
      if (w->stgs->warm_start) 
      {
                  abip_printf("ABIP using variable warm-starting\n");
      }
      
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
      const ABIPData *d, 
      ABIPSolution *sol,
      ABIPWork *w, 
      ABIPInfo *info
) 
{
      DEBUG_FUNC
      
      abip_int i;
      
      char *lin_sys_str = ABIP(get_lin_sys_summary)(w->p, info);
      char *adapt_str = ABIP(get_adapt_summary)(info, w->adapt);
      
      for (i = 0; i < LINE_LEN; ++i) 
      {
                  abip_printf("-");
      }

      abip_printf("\n");

      abip_printf("Status: %s\n", info->status);
      
      if (info->ipm_iter+1 == w->stgs->max_ipm_iters) 
      {
                  abip_printf("Hit max_ipm_iters, solution may be inaccurate\n");
      }

      if (info->admm_iter+1 >= w->stgs->max_admm_iters) 
      {
                  abip_printf("Hit max_admm_iters, solution may be inaccurate\n");
      }

      abip_printf("Timing: Solve time: %1.2es\n", info->solve_time / 1e3);

      if (lin_sys_str) 
      {
                  abip_printf("%s", lin_sys_str);
                  abip_free(lin_sys_str);
      }

      if (adapt_str) 
      {
                  abip_printf("%s", adapt_str);
                  abip_free(adapt_str);
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
                  abip_printf("primal res: |Ax - b|_2 / (1 + |b|_2) = %.4e\n", info->res_pri);
                  abip_printf("dual res: |A'y + s - c|_2 / (1 + |c|_2) = %.4e\n", info->res_dual);
                  abip_printf("rel gap: |c'x - b'y| / (1 + |c'x| + |b'y|) = %.4e\n", info->rel_gap);
                  
                  for (i = 0; i < LINE_LEN; ++i) 
                  {
                                    abip_printf("-");
                  }
                  
                  abip_printf("\n");
                  abip_printf("c'x = %.4e, b'y = %.4e\n", info->pobj, info->dobj);
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
      ABIPWork *w, 
      ABIPResiduals *r, 
      abip_int ipm_iter, 
      abip_int admm_iter
) 
{
      DEBUG_FUNC
      
      abip_float eps = w->stgs->eps;
      
      if (r->res_pri < eps && r->res_dual < eps && r->rel_gap < eps) 
      {
                  RETURN ABIP_SOLVED;
      }
      
      if (r->res_unbdd < eps && ipm_iter > 0 && admm_iter > 0) 
      {
                  RETURN ABIP_UNBOUNDED;
      }
      
      if (r->res_infeas < eps && ipm_iter > 0 && admm_iter > 0) 
      {
                  RETURN ABIP_INFEASIBLE;
      }
      
      RETURN 0;
}

static abip_int validate
(
      const ABIPData *d
) 
{
      DEBUG_FUNC
      
      ABIPSettings *stgs = d->stgs;
      
      if (d->m <= 0 || d->n <= 0) 
      {
                  abip_printf("m and n must both be greater than 0; m = %li, n = %li\n", (long) d->m, (long) d->n);
                  RETURN - 1;
      }
      
      if (d->m > d->n) 
      {
                  abip_printf("WARN: m larger than n, problem likely degenerate\n");
                  RETURN - 1;
      }
      
      if (ABIP(validate_lin_sys)(d->A) < 0) 
      {
                  abip_printf("invalid linear system input data\n");
                  RETURN - 1;
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

      if (stgs->eps <= 0) 
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
      
      if (stgs->scale <= 0) 
      {
                  abip_printf("scale must be positive (1 works well).\n");
                  RETURN - 1;
      }

      if (stgs->eps_cor <= 0) 
      {
                  abip_printf("eps_cor tolerance must be positive\n");
                  RETURN - 1;
      }

      if (stgs->eps_pen <= 0) 
      {
                  abip_printf("eps_pen tolerance must be positive\n");
                  RETURN - 1;
      }

      if (stgs->adaptive_lookback <= 0) 
      {
                  abip_printf("adaptive_lookback must be positive\n");
                  RETURN - 1;
      }
      
      RETURN 0;
}

static ABIPWork *init_work
(     
      const ABIPData *d
) 
{
      DEBUG_FUNC
      
      ABIPWork *w = (ABIPWork *) abip_calloc(1, sizeof(ABIPWork));
      abip_int l = d->n + d->m + 1;

      if (d->stgs->verbose) 
      {
            print_init_header(d);
      }

      if (!w) 
      {
            abip_printf("ERROR: allocating work failure\n");
            RETURN ABIP_NULL;
      }
      
      w->stgs = d->stgs;
      w->m = d->m;
      w->n = d->n;
      
      w->u = (abip_float *) abip_malloc(l * sizeof(abip_float));
      w->v = (abip_float *) abip_malloc(l * sizeof(abip_float));
      w->u_t = (abip_float *) abip_malloc(l * sizeof(abip_float));
      w->u_prev = (abip_float *) abip_malloc(l * sizeof(abip_float));
      w->v_prev = (abip_float *) abip_malloc(l * sizeof(abip_float));
      w->h = (abip_float *) abip_malloc((l - 1) * sizeof(abip_float));
      w->g = (abip_float *) abip_malloc((l - 1) * sizeof(abip_float));
      w->pr = (abip_float *) abip_malloc(d->m * sizeof(abip_float));
      w->dr = (abip_float *) abip_malloc(d->n * sizeof(abip_float));
      w->b = (abip_float *) abip_malloc(d->m * sizeof(abip_float));
      w->c = (abip_float *) abip_malloc(d->n * sizeof(abip_float));

      if (!w->u || !w->v || !w->u_t || !w->u_prev || !w->h || !w->g || !w->pr || !w->dr || !w->b || !w->c) 
      {
                  abip_printf("ERROR: work memory allocation failure\n");
                  RETURN ABIP_NULL;
      }
      
      w->A = d->A;
      w->sp = d->sp;
      
      if (w->stgs->normalize) 
      {
                  #ifdef COPYAMATRIX
                  
                  if (!ABIP(copy_A_matrix)(&(w->A), d->A)) 
                  {
                              abip_printf("ERROR: copy A matrix failed\n");
                              RETURN ABIP_NULL;
                  }
                  
                  #endif

                  w->scal = (ABIPScaling *)abip_malloc(sizeof(ABIPScaling));
                  ABIP(normalize_A)(w->A, w->stgs, w->scal);

                  #if EXTRA_VERBOSE > 0
                  
                  ABIP(print_array)(w->scal->D, d->m, "D");
                  abip_printf("ABIP(norm) D = %4f\n", ABIP(norm)(w->scal->D, d->m));
                  ABIP(print_array)(w->scal->E, d->n, "E");
                  abip_printf("ABIP(norm) E = %4f\n", ABIP(norm)(w->scal->E, d->n));

                  #endif
      } 
      else 
      {
                  w->scal = ABIP_NULL;
      }

      if (!(w->p = ABIP(init_lin_sys_work)(w->A, w->stgs))) 
      {
                  abip_printf("ERROR: init_lin_sys_work failure\n");
                  RETURN ABIP_NULL;
      }
      
      if (!(w->adapt = ABIP(init_adapt)(w))) 
      {
                  abip_printf("ERROR: init_adapt failure\n");
                  RETURN ABIP_NULL;
      }

      RETURN w;
}

static abip_int update_work
(
      const ABIPData *d, 
      ABIPWork *w,
      const ABIPSolution *sol
) 
{
      DEBUG_FUNC
      
      abip_int n = d->n;
      abip_int m = d->m; 

      w->nm_b = ABIP(norm)(d->b, m);
      w->nm_c = ABIP(norm)(d->c, n);
      memcpy(w->b, d->b, d->m * sizeof(abip_float));
      memcpy(w->c, d->c, d->n * sizeof(abip_float));

      #if EXTRA_VERBOSE > 0
      
      ABIP(print_array)(w->b, m, "b");
      abip_printf("pre-normalized norm b = %4f\n", ABIP(norm)(w->b, m));
      ABIP(print_array)(w->c, n, "c");
      abip_printf("pre-normalized norm c = %4f\n", ABIP(norm)(w->c, n));
      
      #endif
      
      if (w->stgs->normalize) 
      {
                  ABIP(normalize_b_c)(w);
                  
                  #if EXTRA_VERBOSE > 0
                  
                  ABIP(print_array)(w->b, m, "bn");
                  abip_printf("sc_b = %4f\n", w->sc_b);
                  abip_printf("post-normalized norm b = %4f\n", ABIP(norm)(w->b, m));
                  
                  ABIP(print_array)(w->c, n, "cn");
                  abip_printf("sc_c = %4f\n", w->sc_c);
                  abip_printf("post-normalized norm c = %4f\n", ABIP(norm)(w->c, n));
                  
                  #endif
      }
      
      if (MAX(w->sp,w->stgs->sparsity_ratio) > 0.4 || (MIN(w->sp,w->stgs->sparsity_ratio) > 0.1 && MIN(w->sp,w->stgs->sparsity_ratio) < 0.2))
      {
            w->sigma = 0.3; 
            w->gamma = 2.0; 
      }
      else if (MIN(w->sp,w->stgs->sparsity_ratio) > 0.2)
      {
            w->sigma = 0.5; 
            w->gamma = 3.0;
      }
      else
      {
            w->sigma = 0.8; 
            w->gamma = 3.0;
      }

      w->final_check = 0; 
      w->double_check = 0; 

      w->mu = 1.0; 
      w->beta = 1.0; 
      
      if (w->stgs->warm_start) 
      {
                  warm_start_vars(w, sol);
      } 
      else 
      {
                  cold_start_vars(w);
      }
      
      memcpy(w->h, w->b, m * sizeof(abip_float));
      memcpy(&(w->h[m]), w->c, n * sizeof(abip_float));
      ABIP(scale_array)(w->h, -1, m);
      memcpy(w->g, w->h, (n + m) * sizeof(abip_float));
      
      ABIP(solve_lin_sys)(w->A, w->stgs, w->p, w->g, ABIP_NULL, -1);
      ABIP(scale_array)(&(w->g[m]), -1, n);
      w->g_th = ABIP(dot)(w->h, w->g, n + m);
      
      RETURN 0;
}

static abip_float iterate_norm_diff
(     
      ABIPWork *w
) 
{
      DEBUG_FUNC
      
      abip_int l = w->m + w->n + 1;
      
      abip_float u_norm_difference = ABIP(norm_diff)(w->u, w->u_prev, l);
      abip_float v_norm_difference = ABIP(norm_diff)(w->v, w->v_prev, l);
      abip_float norm = 1 + SQRTF(ABIP(norm_sq)(w->u, l) + ABIP(norm_sq)(w->v, l)) + SQRTF(ABIP(norm_sq)(w->u_prev, l) + ABIP(norm_sq)(w->v_prev, l));
      abip_float norm_diff = SQRTF(u_norm_difference * u_norm_difference + v_norm_difference * v_norm_difference);
      
      RETURN norm_diff / norm;
}

static abip_float iterate_Q_norm_resd
(     
      ABIPWork *w
) 
{
      DEBUG_FUNC

      abip_int i;
      abip_int l = w->m + w->n + 1;

      abip_float *y = w->u; 
      abip_float *x = &(w->u[w->m]); 
      abip_float *s = &(w->v[w->m]);
      abip_float tau = w->u[w->m + w->n];
      abip_float kap = w->v[w->m + w->n]; 

      abip_float Qres = 0; 
      abip_float *pr = w->pr;;
      abip_float *dr = w->dr;;
           
      memset(pr, 0, w->m * sizeof(abip_float));
      memset(dr, 0, w->n * sizeof(abip_float));
      ABIP(accum_by_A)(w->A, w->p, x, pr);
      ABIP(accum_by_Atrans)(w->A, w->p, y, dr);      
      ABIP(add_scaled_array)(dr, s, w->n, 1.0);
      
      for (i = 0; i < w->m; ++i) 
      {
            Qres += (pr[i] - w->b[i] * tau) * (pr[i] - w->b[i] * tau);
      }
      
      for (i = 0; i < w->n; ++i) 
      {
            Qres += (dr[i] - w->c[i] * tau) * (dr[i] - w->c[i] * tau);
      }

      abip_float cTx = ABIP(dot)(x, w->c, w->n); 
      abip_float bTy = ABIP(dot)(y, w->b, w->m); 
      Qres += (bTy - cTx - kap) * (bTy - cTx - kap); 
      abip_float norm = 1 + SQRTF(ABIP(norm_sq)(w->u, l) + ABIP(norm_sq)(w->v, l));
      
      RETURN SQRTF(Qres) / norm; 
}

abip_int ABIP(solve)
(
      ABIPWork *w, 
      const ABIPData *d,
      ABIPSolution *sol, 
      ABIPInfo *info
) 
{
      DEBUG_FUNC
      
      abip_int i;
      abip_int j; 
      abip_int k; 
      abip_int inner_stopper; 
      ABIP(timer) solve_timer;
      
      ABIPResiduals r;
      abip_int l = w->m + w->n + 1;
      
      if (!d || !sol || !info || !w || !d->b || !d->c) 
      {
                  abip_printf("ERROR: ABIP_NULL input\n");
                  RETURN ABIP_FAILED;
      }

      abip_start_interrupt_listener();
      ABIP(tic)(&solve_timer);
      
      info->status_val = ABIP_UNFINISHED; 
      r.last_ipm_iter = -1;
      r.last_admm_iter = -1; 
      update_work(d, w, sol);

      if (w->stgs->verbose) 
      {
                  print_header(w);
      }

      k = 0;
      
      for (i = 0; i < w->stgs->max_ipm_iters; ++i)
      {
                  if (MIN(w->sp,w->stgs->sparsity_ratio) > 0.5)
                  {
                        inner_stopper = (int)round(POWF(w->mu, -0.35));
                  }
                  else if (MIN(w->sp,w->stgs->sparsity_ratio) > 0.2)
                  {
                        inner_stopper = (int)round(POWF(w->mu, -1));
                  }
                  else
                  {
              
                        inner_stopper = w->stgs->max_admm_iters; 
                  }  
              
                  for (j = 0; j < inner_stopper; ++j) 
                  {
                                      memcpy(w->u_prev, w->u, l * sizeof(abip_float));
                                      memcpy(w->v_prev, w->v, l * sizeof(abip_float));

                                      if (project_lin_sys(w, k) < 0) 
                                      {
                                                      RETURN failure(w, w->m, w->n, sol, info, ABIP_FAILED, "error in project_lin_sys", "Failure");
                                      }
                              
                                      project_barrier(w);
                               
                                      update_dual_vars(w);

                                      if (abip_is_interrupted()) 
                                      {
                                                      RETURN failure(w, w->m, w->n, sol, info, ABIP_SIGINT, "Interrupted", "Interrupted");
                                      }
                                      
                                      /* abip_printf("||Qu-v|| = %3.6f, gamma = %3.6f, mu = %3.6f\n", iterate_Q_norm_resd(w), w->gamma, w->mu);  */
                                      
                                      k += 1; 

                                      if (iterate_Q_norm_resd(w) < w->gamma*w->mu)
                                      {
                                                      break; 
                                      }

                                      if (w-> final_check && (j+1) % CONVERGED_INTERVAL == 0)
                                      {
                                                        calc_residuals(w, &r, i, k);
                                                      
                                                        if ((info->status_val = has_converged(w, &r, i, k)) != 0 || k+1 >= w->stgs->max_admm_iters || i+1 >= w->stgs->max_ipm_iters) 
                                                        {
                                                                            if (w->stgs->verbose && k>0) 
                                                                            {
                                                                                    print_summary(w, i, k, &r, &solve_timer);
                                                                            }
                                                                            
                                                                            get_solution(w, sol, info, &r, i, k);
                                                                            info->solve_time = ABIP(tocq)(&solve_timer);

                                                                            if (w->stgs->verbose) 
                                                                            {
                                                                                    print_footer(d, sol, w, info);
                                                                            }
                  
                                                                            abip_end_interrupt_listener();

                                                                            RETURN info->status_val;
                                                        }
                                      }
                  }
                  

                  if (k > (int) (w->stgs->max_admm_iters*0.8))
                  {
                        w-> final_check = 1; 
                  }
                  
                  calc_residuals(w, &r, i, k);
                  if (w->stgs->verbose)
                  {
                              print_summary(w, i, k, &r, &solve_timer);
                  }

                  if ((info->status_val = has_converged(w, &r, i, k)) != 0 || k+1 >= w->stgs->max_admm_iters) 
                  {
                                      get_solution(w, sol, info, &r, i, k);
                                      info->solve_time = ABIP(tocq)(&solve_timer);

                                      if (w->stgs->verbose) 
                                      {
                                                         print_footer(d, sol, w, info);
                                      }
                  
                                      abip_end_interrupt_listener();

                                      RETURN info->status_val;
                  }

                  update_barrier(w, &r); 

                  reinitialize_vars(w); 

                  if (w->stgs->adaptive) 
                  {
                              for (i = w->m; i < l; ++i) 
                              {
                                    w->u[i] = SQRTF(w->beta) * w->u[i]; 
                                    w->v[i] = SQRTF(w->beta) * w->v[i]; 
                              }
                              w->beta = 1;                               
                              if (ABIP(adaptive)(w, k) < 0) 
                              {
                                                RETURN failure(w, w->m, w->n, sol, info, ABIP_FAILED, "error in adaptive", "Failure");
                              }
                              for (i = w->m; i < l; ++i) 
                              {
                                    w->u[i] = SQRTF(1.0/w->beta) * w->u[i]; 
                                    w->v[i] = SQRTF(1.0/w->beta) * w->v[i]; 
                              }
                  }
        }
        RETURN info->status_val;
}

void ABIP(finish)
(
      ABIPWork *w
) 
{
      DEBUG_FUNC
      
      if (w) 
      {
                  if (w->stgs && w->stgs->normalize) 
                  {
                              #ifndef COPYAMATRIX
                              ABIP(un_normalize_A)(w->A, w->stgs, w->scal);
                              #else
                              ABIP(free_A_matrix)(w->A);
                              #endif
                  }
                  
                  if (w->p) 
                  {
                              ABIP(free_lin_sys_work)(w->p);
                  }
                  
                  if (w->adapt) 
                  {
                              ABIP(free_adapt)(w->adapt);
                  }
                  
                  free_work(w);
      }
      
      RETURN;
}

ABIPWork *ABIP(init)
(
      const ABIPData *d, 
      ABIPInfo *info
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
      if (validate(d) < 0) 
      {
                  abip_printf("ERROR: Validation returned failure\n");
                  RETURN ABIP_NULL;
      }
      #endif
      
      ABIP(tic)(&init_timer);
      
      w = init_work(d);
      info->setup_time = ABIP(tocq)(&init_timer);
      
      if (d->stgs->verbose) 
      {
                  abip_printf("Setup time: %1.2es\n", info->setup_time / 1e3);
      }
      
      abip_end_interrupt_listener();
      
      RETURN w;
}

abip_int ABIP(main)
(
      const ABIPData *d, 
      ABIPSolution *sol,
      ABIPInfo *info
) 
{
      DEBUG_FUNC
      
      abip_int status;
      ABIPWork *w = ABIP(init)(d, info);
      
      #if EXTRA_VERBOSE > 0
      abip_printf("size of abip_int = %lu, size of abip_float = %lu\n", sizeof(abip_int), sizeof(abip_float));
      #endif
      
      if (w) 
      {
                  ABIP(solve)(w, d, sol, info);
                  status = info->status_val;
      } 
      else 
      {
                  status = failure(ABIP_NULL, d ? d->m : -1, d ? d->n : -1, sol, info, ABIP_FAILED, "could not initialize work", "Failure");
      }
      
      ABIP(finish)(w);
      
      RETURN status;
}
