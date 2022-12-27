#include "adaptive.h"
#include "linalg.h"
#include "linsys.h"
#include "abip.h"
#include "abip_blas.h"
#include "util.h"

/* This file uses adaption to improve the convergence rate of the ADMM in each inner loop. 
 * At each iteration we need to select a nearly-optimal penalty parameter beta, we do this using Barzilai-Borwein spectral method.
 * Adaptive_lookback is the number of lookback iterations.
 */

struct ABIP_ADAPTIVE_WORK 
{
    abip_float *u_prev;
    abip_float *v_prev;
    abip_float *ut;
    abip_float *u;
    abip_float *v;
    abip_float *ut_next;
    abip_float *u_next;
    abip_float *v_next;
    
    abip_float *delta_ut;
    abip_float *delta_u;
    abip_float *delta_v;
    
    abip_int l;
    abip_int k;
    
    abip_float total_adapt_time;
};

static abip_int update_adapt_params
 (
  ABIPWork *w,
  abip_int iter
  )
{
    DEBUG_FUNC
    
    abip_float *u_prev = w->adapt->u_prev;
    abip_float *v_prev = w->adapt->v_prev;
    abip_float *ut = w->adapt->ut;
    abip_float *u = w->adapt->u;
    abip_float *v = w->adapt->v;
    abip_float *ut_next = w->adapt->ut_next;
    abip_float *u_next = w->adapt->u_next;
    abip_float *v_next = w->adapt->v_next;
    
    abip_float *delta_ut = w->adapt->delta_ut;
    abip_float *delta_u = w->adapt->delta_u;
    abip_float *delta_v = w->adapt->delta_v;
    
    abip_int l = w->adapt->l;
    abip_int n = w->n;
    abip_int m = w->m;
    abip_int k = w->adapt->k;
    abip_int i;
    abip_int j;
    
    abip_int status_1 = 0;
    abip_int status_2 = 0;
    
    abip_float tmp;
    abip_float beta_prev = 1.0;
    abip_float beta = 0.0;
    
    abip_float uu;
    abip_float uv;
    abip_float vv;
    abip_float utut;
    abip_float utv;
    abip_float norm_ut;
    abip_float norm_u;
    abip_float norm_v;
    
    abip_float alpha_SD;
    abip_float alpha_MG;
    abip_float gamma_SD;
    abip_float gamma_MG;
    abip_float alpha_cor;
    abip_float gamma_cor;
    abip_float alpha_ss;
    abip_float gamma_ss;
    
    memcpy(u_prev, w->u, sizeof(abip_float) * l);
    memcpy(v_prev, w->v, sizeof(abip_float) * l);
    
    for (i = 0; i < k; ++i)
    {
        memcpy(ut, u_prev, l * sizeof(abip_float));
        ABIP(add_scaled_array)(ut, v_prev, l, 1.0);
        ABIP(scale_array)(ut, w->stgs->rho_y, m);
        ABIP(add_scaled_array)(ut, w->h, l - 1, -ut[l - 1]);
        ABIP(add_scaled_array)(ut, w->h, l - 1, -ABIP(dot)(ut, w->g, l - 1) / (w->g_th + 1));
        ABIP(scale_array)(&(ut[m]), -1, n);
        status_1 = ABIP(solve_lin_sys)(w->A, w->stgs, w->p, ut, u_prev, iter);
        ut[l - 1] += ABIP(dot)(ut, w->h, l - 1);
        
        for (j = 0; j < m; ++j)
        {
            u[j] = ut[j] - v_prev[j];
        }
        
        for (j = m; j < l; ++j)
        {
            u[j] = w->stgs->alpha * ut[j] + (1 - w->stgs->alpha) * u_prev[j] - v_prev[j];
        }
        
        for(j = m; j < l; ++j)
        {
            tmp = u[j] / 2;
            u[j] = tmp + SQRTF(tmp * tmp + w->mu / beta_prev);
        }
        
        
        for (j = m; j < l; ++j)
        {
            v[j] =  v_prev[j] + (u[j] - w->stgs->alpha * ut[j] - (1 - w->stgs->alpha) * u_prev[j]);
        }
        
        
        memcpy(ut_next, u, l * sizeof(abip_float));
        ABIP(add_scaled_array)(ut_next, v, l, 1.0);
        ABIP(scale_array)(ut_next, w->stgs->rho_y, m);
        ABIP(add_scaled_array)(ut_next, w->h, l - 1, -ut_next[l - 1]);
        ABIP(add_scaled_array)(ut_next, w->h, l - 1, -ABIP(dot)(ut_next, w->g, l - 1) / (w->g_th + 1));
        ABIP(scale_array)(&(ut_next[m]), -1, n);
        status_2 = ABIP(solve_lin_sys)(w->A, w->stgs, w->p, ut_next, u, iter);
        ut_next[l - 1] += ABIP(dot)(ut_next, w->h, l - 1);
        
        for (j = 0; j < m; ++j)
        {
            u_next[j] = ut_next[j] - v[j];
        }
        
        for (j = m; j < l; ++j)
        {
            u_next[j] = w->stgs->alpha * ut_next[j] + (1 - w->stgs->alpha) * u[j] - v[j];
        }
        
        for(j = m; j < l; ++j)
        {
            tmp = u_next[j] / 2;
            u_next[j] = tmp + SQRTF(tmp * tmp + w->mu / beta_prev);
        }
        
        for (j = m; j < l; ++j)
        {
            v_next[j] =  v[j] + (u_next[j] - w->stgs->alpha * ut_next[j] - (1 - w->stgs->alpha) * u[j]);
        }
        
        memcpy(delta_ut, v, l * sizeof(abip_float));
        ABIP(scale_array)(delta_ut, 2.0, l);
        ABIP(add_scaled_array)(delta_ut, u_next, l, 1.0);
        ABIP(add_scaled_array)(delta_ut, u, l, -1.0);
        ABIP(add_scaled_array)(delta_ut, v_next, l, -1.0);
        ABIP(add_scaled_array)(delta_ut, v_prev, l, -1.0);
        
        memcpy(delta_u, u, l * sizeof(abip_float));
        ABIP(add_scaled_array)(delta_u, u_next, l, -1.0);
        
        memcpy(delta_v, u_next, l * sizeof(abip_float));
        ABIP(add_scaled_array)(delta_v, u, l, -1.0);
        ABIP(scale_array)(delta_v, w->stgs->alpha-1.0, l);
        ABIP(add_scaled_array)(delta_v, v_next, l, 1.0);
        ABIP(add_scaled_array)(delta_v, v, l, -1.0);
        
        utut = ABIP(dot)(delta_ut, delta_ut, l);
        utv = ABIP(dot)(delta_ut, delta_v, l);
        uu = ABIP(dot)(delta_u, delta_u, l);
        vv = ABIP(dot)(delta_v, delta_v, l);
        uv = ABIP(dot)(delta_u, delta_v, l);
        
        norm_ut = ABIP(norm)(delta_ut, l);
        norm_u = ABIP(norm)(delta_u, l);
        norm_v = ABIP(norm)(delta_v, l);
        
        alpha_SD = vv/utv;
        alpha_MG = utv/utut;
        gamma_SD = vv/uv;
        gamma_MG = uv/uu;
        /*
         abip_printf("alpha_SD = %3.6f, alpha_MG = %3.6f, gamma_SD = %3.6f, gamma_MG = %3.6f\n", alpha_SD, alpha_MG, gamma_SD, gamma_MG);
         */
        if (2*alpha_MG > alpha_SD)
        {
            alpha_ss = alpha_MG;
        }
        else
        {
            alpha_ss = alpha_SD - 0.5*alpha_MG;
        }
        
        if (2*gamma_MG > gamma_SD)
        {
            gamma_ss = gamma_MG;
        }
        else
        {
            gamma_ss = gamma_SD - 0.5*gamma_MG;
        }
        
        alpha_cor = utv / (norm_v*norm_ut);
        gamma_cor = uv / (norm_v*norm_u);
        
        if (alpha_cor > w->stgs->eps_cor && gamma_cor > w->stgs->eps_cor)
        {
            beta = SQRTF(alpha_ss*gamma_ss);
        }
        else if (alpha_cor > w->stgs->eps_cor && gamma_cor <= w->stgs->eps_cor)
        {
            beta = alpha_ss;
        }
        else if (alpha_cor <= w->stgs->eps_cor && gamma_cor > w->stgs->eps_cor)
        {
            beta = gamma_ss;
        }
        else
        {
            beta = beta_prev;
        }
        
        if (ABS(beta-beta_prev) > 0 && ABS(beta-beta_prev) <= w->stgs->eps_pen)
        {
            beta = (beta+beta_prev)/2;
            break;
        }
        else if (ABS(beta-beta_prev) > w->stgs->eps_pen)
        {
            beta_prev = beta;
            memcpy(u_prev, u, l * sizeof(abip_float));
            for (j = 0; j < m; ++j)
            {
                v_prev[j] = v[j];
            }
            for (j = m; j < l; ++j)
            {
                v_prev[j] = (w->mu / beta_prev) / u_prev[j];
            }
        }
        else
        {
            memcpy(u_prev, u, l * sizeof(abip_float));
            memcpy(v_prev, v, l * sizeof(abip_float));
        }
        /*
         abip_printf("beta = %3.7e, vv = %3.7e, uu = %3.7e, utv = %3.7e, uv = %3.7e, utut = %3.7e\n", beta, vv, uu, utv, uv, utut);
         */
    }
    
    w->beta = beta;
    
    RETURN MIN(status_1, status_2) ;
}

ABIPAdaptWork *ABIP(init_adapt)
(
 ABIPWork *w
 )
{
    DEBUG_FUNC
    
    ABIPAdaptWork *a = (ABIPAdaptWork *) abip_calloc(1, sizeof(ABIPAdaptWork));
    
    if (!a)
    {
        RETURN ABIP_NULL;
    }
    
    a->l = w->m + w->n + 1;
    
    a->k = w->stgs->adaptive_lookback;
    
    if (a->k <= 0)
    {
        RETURN a;
    }
    
    a->u_prev = (abip_float *) abip_calloc(a->l, sizeof(abip_float));
    a->v_prev = (abip_float *) abip_calloc(a->l, sizeof(abip_float));
    a->ut = (abip_float *) abip_calloc(a->l, sizeof(abip_float));
    a->u = (abip_float *) abip_calloc(a->l, sizeof(abip_float));
    a->v = (abip_float *) abip_calloc(a->l, sizeof(abip_float));
    a->ut_next = (abip_float *) abip_calloc(a->l, sizeof(abip_float));
    a->u_next = (abip_float *) abip_calloc(a->l, sizeof(abip_float));
    a->v_next = (abip_float *) abip_calloc(a->l, sizeof(abip_float));
    
    a->delta_ut = (abip_float *) abip_calloc(a->l, sizeof(abip_float));
    a->delta_u = (abip_float *) abip_calloc(a->l, sizeof(abip_float));
    a->delta_v = (abip_float *) abip_calloc(a->l, sizeof(abip_float));
    
    a->total_adapt_time = 0.0;
    
    if (!a->u_prev || !a->v_prev || !a->ut || !a->u || !a->v || !a->ut_next || !a->u_next || !a->v_next || !a->delta_ut || !a->delta_u || !a->delta_v)
    {
        ABIP(free_adapt)(a);
        a = ABIP_NULL;
    }
    
    RETURN a;
}

abip_int ABIP(adaptive)
(
 ABIPWork *w,
 abip_int iter
 )
{
    DEBUG_FUNC
    
    abip_int k = w->adapt->k;
    abip_int info;
    
    ABIP(timer) adapt_timer;
    if (k <= 0)
    {
        RETURN -1;
    }
    
    ABIP(tic)(&adapt_timer);
    
    info = update_adapt_params(w, iter);
    
    if (iter == -1)
    {
        RETURN -1;
    }
    
    w->adapt->total_adapt_time += ABIP(tocq)(&adapt_timer);
    
    RETURN info;
}

void ABIP(free_adapt)
(
 ABIPAdaptWork *a
 )
{
    DEBUG_FUNC
    
    if (a)
    {
        if (a->u_prev)
        {
            abip_free(a->u_prev);
        }
        
        if (a->v_prev)
        {
            abip_free(a->v_prev);
        }
        
        if (a->ut)
        {
            abip_free(a->ut);
        }
        
        if (a->u)
        {
            abip_free(a->u);
        }
        
        if (a->v)
        {
            abip_free(a->v);
        }
        
        if (a->ut_next)
        {
            abip_free(a->ut_next);
        }
        
        if (a->u_next)
        {
            abip_free(a->u_next);
        }
        
        if (a->v_next)
        {
            abip_free(a->v_next);
        }
        
        if (a->delta_ut)
        {
            abip_free(a->delta_ut);
        }
        
        if (a->delta_u)
        {
            abip_free(a->delta_u);
        }
        
        if (a->delta_v)
        {
            abip_free(a->delta_v);
        }
        
        abip_free(a);
    }
    
    RETURN;
}

char *ABIP(get_adapt_summary)
(
 const ABIPInfo *info,
 ABIPAdaptWork *a
 )
{
    DEBUG_FUNC
    
    char *str = (char *) abip_malloc(sizeof(char) * 64);
    sprintf(str, "\tBarzilai-Borwein spectral method: avg step time: %1.2es\n", a->total_adapt_time / (info->admm_iter + 1) / 1e3);
    
    a->total_adapt_time = 0.0;
    RETURN str;
}
