#include "abip.h"
#include "normalize.h"
#include "linalg.h"

#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)

/**
@brief normalize b and c
*/
void ABIP(normalize_b_c)
(
 ABIPWork *w
 )
{
    abip_int i;
    
    abip_float nm;
    abip_float *D = w->scal->D;
    abip_float *E = w->scal->E;
    abip_float *b = w->b;
    abip_float *c = w->c;
    
    for (i = 0; i < w->n; ++i)
    {
        c[i] /= E[i];
    }
    nm = ABIP(norm)(c, w->n);
    w->sc_c = w->scal->mean_norm_row_A / MAX(nm, MIN_SCALE);
    
    for (i = 0; i < w->m; ++i)
    {
        b[i] /= D[i];
    }
    nm = ABIP(norm)(b, w->m);
    w->sc_b = w->scal->mean_norm_col_A / MAX(nm, MIN_SCALE);
    
    ABIP(scale_array)(c, w->sc_c * w->stgs->scale, w->n);
    ABIP(scale_array)(b, w->sc_b * w->stgs->scale, w->m);
}
/**
@brief calculate the scaled residuals
*/
void ABIP(calc_scaled_resids)
(
 ABIPWork *w,
 ABIPResiduals *r
 )
{
    abip_float *D = w->scal->D;
    abip_float *E = w->scal->E;
    
    abip_float *u = w->u;
    abip_float *u_t = w->u_t;
    abip_float *u_prev = w->u_prev;
    abip_float tmp;
    
    abip_int i;
    abip_int n = w->n;
    abip_int m = w->m;
    
    r->res_pri = 0;
    for (i = 0; i < m; ++i)
    {
        tmp = (u[i] - u_t[i]) / (D[i] * w->sc_c);
        r->res_pri += tmp * tmp;
    }
    
    for (i = 0; i < n; ++i)
    {
        tmp = (u[i + m] - u_t[i + m]) / (E[i] * w->sc_b);
        r->res_pri += tmp * tmp;
    }
    
    tmp = u[n + m] - u_t[n + m];
    r->res_pri += tmp * tmp;
    r->res_pri = sqrt(r->res_pri);
    
    r->res_dual = 0;
    for (i = 0; i < m; ++i)
    {
        tmp = (u[i] - u_prev[i]) * D[i] / w->sc_c;
        r->res_dual += tmp * tmp;
    }
    
    for (i = 0; i < n; ++i)
    {
        tmp = (u[i + m] - u_prev[i + m]) * E[i] / w->sc_b;
        r->res_dual += tmp * tmp;
    }
    
    tmp = u[n + m] - u_prev[n + m];
    r->res_dual += tmp * tmp;
    r->res_dual = sqrt(r->res_dual);
}

/**
@brief normalize the warm start solution
*/
void ABIP(normalize_warm_start)
(
 ABIPWork *w
 )
{
    abip_int i;
    
    abip_float *D = w->scal->D;
    abip_float *E = w->scal->E;
    
    abip_float *y = w->u;
    abip_float *x = &(w->u[w->m]);
    abip_float *s = &(w->v[w->m]);
    
    for (i = 0; i < w->n; ++i)
    {
        x[i] *= (E[i] * w->sc_b);
    }
    
    for (i = 0; i < w->m; ++i)
    {
        y[i] *= (D[i] * w->sc_c);
    }
    
    for (i = 0; i < w->n; ++i)
    {
        s[i] /= (E[i] / (w->sc_c * w->stgs->scale));
    }
}

/**
@brief recover the optimal solution
*/
void ABIP(un_normalize_sol)
(
 ABIPWork *w,
 ABIPSolution *sol
 )
{
    abip_int i;
    
    abip_float *D = w->scal->D;
    abip_float *E = w->scal->E;
    
    for (i = 0; i < w->n; ++i)
    {
        sol->x[i] /= (E[i] * w->sc_b);
    }
    
    for (i = 0; i < w->m; ++i)
    {
        sol->y[i] /= (D[i] * w->sc_c);
    }
    
    for (i = 0; i < w->n; ++i)
    {
        sol->s[i] *= E[i] / (w->sc_c * w->stgs->scale);
    }
}
