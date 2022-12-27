#include "glbopts.h"
#include "util.h"
#include "linsys.h"

#if (defined NOTIMER)

void ABIP(tic)
(
 ABIP(timer) *t
 ) {}

abip_float ABIP(tocq)
(
 ABIP(timer) *t
 ) 
{
    return NAN;
}

#elif(defined _WIN32 || _WIN64 || defined _WINDLL)

void ABIP(tic)
(
 ABIP(timer) *t
 ) 
{
    QueryPerformanceFrequency(&t->freq);
    QueryPerformanceCounter(&t->tic);
}

abip_float ABIP(tocq)
(
 ABIP(timer) *t
 ) 
{
    QueryPerformanceCounter(&t->toc);
    return (1e3 * (t->toc.QuadPart - t->tic.QuadPart) / (abip_float)t->freq.QuadPart);
}

#elif(defined __APPLE__)

void ABIP(tic)
(
 ABIP(timer) *t
 ) 
{
    /* read current clock cycles */
    t->tic = mach_absolute_time();
}

abip_float ABIP(tocq)
(
 ABIP(timer) *t
 ) 
{
    uint64_t duration;
    
    t->toc = mach_absolute_time();
    duration = t->toc - t->tic;
    
    mach_timebase_info(&(t->tinfo));
    duration *= t->tinfo.numer;
    duration /= t->tinfo.denom;
    
    return (abip_float)duration / 1e6;
}

#else

void ABIP(tic)
(
 ABIP(timer) *t
 ) 
{
    clock_gettime(CLOCK_MONOTONIC, &t->tic);
}

abip_float ABIP(tocq)
(
 ABIP(timer) *t
 ) 
{
    struct timespec temp;
    
    clock_gettime(CLOCK_MONOTONIC, &t->toc);
    
    if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) 
    {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
        temp.tv_nsec = 1e9 + t->toc.tv_nsec - t->tic.tv_nsec;
    } 
    else 
    {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
        temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
    }
    
    return (abip_float) temp.tv_sec * 1e3 + (abip_float) temp.tv_nsec / 1e6;
}

#endif

abip_float ABIP(toc)
(
 ABIP(timer) *t
 ) 
{
    abip_float time = ABIP(tocq)(t);
    abip_printf("time: %8.4f milli-seconds.\n", time);
    return time;
}

abip_float ABIP(str_toc)
(
 char *str, 
 ABIP(timer) *t
 ) 
{
    abip_float time = ABIP(tocq)(t);
    abip_printf("%s - time: %8.4f milli-seconds.\n", str, time);
    return time;
}

void ABIP(print_work)
(
 const ABIPWork *w
 ) 
{
    abip_int i; 
    abip_int l = w->n + w->m;
    
    abip_printf("\n u_t is \n");
    for (i = 0; i < l; i++) 
    {
        abip_printf("%f\n", w->u_t[i]);
    }
    
    abip_printf("\n u is \n");
    for (i = 0; i < l; i++) 
    {
        abip_printf("%f\n", w->u[i]);
    }
    
    abip_printf("\n v is \n");
    for (i = 0; i < l; i++) 
    {
        abip_printf("%f\n", w->v[i]);
    }
}

void ABIP(print_data)
(
 const ABIPData *d
 ) 
{
    abip_printf("m = %i\n", (int)d->m);
    abip_printf("n = %i\n", (int)d->n);
    
    abip_printf("max_ipm_iters = %i\n", (int)d->stgs->max_ipm_iters);
    abip_printf("max_admm_iters = %i\n", (int)d->stgs->max_admm_iters);
    
    abip_printf("verbose = %i\n", (int)d->stgs->verbose);
    abip_printf("normalize = %i\n", (int)d->stgs->normalize);
    abip_printf("warm_start = %i\n", (int)d->stgs->warm_start); 
    abip_printf("adaptive = %i\n", (int)d->stgs->adaptive);
    abip_printf("adaptive_lookback = %i\n", (int)d->stgs->adaptive_lookback);
    
    abip_printf("eps = %4f\n", d->stgs->eps);
    abip_printf("alpha = %4f\n", d->stgs->alpha);
    abip_printf("rho_y = %4f\n", d->stgs->rho_y);
    abip_printf("scale = %4f\n", d->stgs->scale);
    
    abip_printf("eps_cor = %4f\n", d->stgs->eps_cor);
    abip_printf("eps_pen = %4f\n", d->stgs->eps_pen);
    
}

void ABIP(print_array)
(
 const abip_float *arr, 
 abip_int n, 
 const char *name
 ) 
{
    abip_int i; 
    abip_int j; 
    abip_int k = 0;
    
    abip_int num_on_one_line = 10;
    
    abip_printf("\n");
    for (i = 0; i < n / num_on_one_line; ++i) 
    {
        for (j = 0; j < num_on_one_line; ++j) 
        {
            abip_printf("%s[%li] = %4f, ", name, (long)k, arr[k]);
            k++;
        }
        abip_printf("\n");
    }
    
    for (j = k; j < n; ++j) 
    {
        abip_printf("%s[%li] = %4f, ", name, (long)j, arr[j]);
    }
    
    abip_printf("\n");
}

void ABIP(free_data)
(
 ABIPData *d
 ) 
{
    if (d) 
    {
        if (d->b) 
        {
            abip_free(d->b);
        }
        
        if (d->c) 
        {
            abip_free(d->c);
        }
        
        if (d->stgs) 
        {
            abip_free(d->stgs);
        }
        
        if (d->A) 
        {
            ABIP(free_A_matrix)(d->A);
        }
        
        abip_free(d);
    }
}

void ABIP(free_sol)
(
 ABIPSolution *sol
 ) 
{
    if (sol) 
    {
        if (sol->x) 
        {
            abip_free(sol->x);
        }
        
        if (sol->y) 
        {
            abip_free(sol->y);
        }
        
        if (sol->s) 
        {
            abip_free(sol->s);
        }
        
        abip_free(sol);
    }
}

void ABIP(set_default_settings)
(
 ABIPData *d
 ) 
{
    d->stgs->max_ipm_iters = MAX_IPM_ITERS;                         
    d->stgs->max_admm_iters = MAX_ADMM_ITERS;                
    d->stgs->eps = EPS;                                                              
    d->stgs->alpha = ALPHA;
    d->stgs->cg_rate = CG_RATE;                                                       
    
    d->stgs->normalize = NORMALIZE;                                       
    d->stgs->scale = SCALE;                                                        
    d->stgs->rho_y = RHO_Y;
    d->stgs->sparsity_ratio = SPARSITY_RATIO; 
    
    d->stgs->adaptive = ADAPTIVE;                                            
    d->stgs->eps_cor = EPS_COR;                                               
    d->stgs->eps_pen = EPS_PEN;                                               
    d->stgs->adaptive_lookback = ADAPTIVE_LOOKBACK;
    
    d->stgs->dynamic_x = 0.8;
    d->stgs->dynamic_eta = 1.1;
    
    d->stgs->restart_fre = 1000;
    d->stgs->restart_thresh = 100000;
    
    // add by Kurt. 22.05.03
    d->stgs->origin_rescale = 0;
    d->stgs->pc_ruiz_rescale = 1;
    d->stgs->qp_rescale = 0;
    d->stgs->ruiz_iter = 10;
    d->stgs->hybrid_mu = 1;
    d->stgs->dynamic_sigma = -1.0;
    d->stgs->hybrid_thresh = 1000;
    d->stgs->dynamic_sigma_second = 0.5;
    
    d->stgs->half_update = 0;
    d->stgs->avg_criterion = 0;
    
    d->stgs->verbose = VERBOSE;                                             
    d->stgs->warm_start = WARM_START;
}
