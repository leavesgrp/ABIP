#include "indirect.h"

#define CG_BEST_TOL 1e-9
#define CG_MIN_TOL 1e-1

// use cg to solve the linear system

char *ABIP(get_lin_sys_method)
(
 const ABIPMatrix *A,
 const ABIPSettings *stgs
 )
{
    char *str = (char *)abip_malloc(sizeof(char) * 128);
    sprintf(str, "sparse-indirect, nnz in A = %li, CG tol ~ 1/iter^(%2.2f)",
            (long)A->p[A->n], stgs->cg_rate);
    return str;
}

char *ABIP(get_lin_sys_summary)
(
 ABIPLinSysWork *p,
 const ABIPInfo *info
 )
{
    char *str = (char *)abip_malloc(sizeof(char) * 128);
    sprintf(str, "\tLin-sys: avg # CG iterations: %2.2f, avg solve time: %1.2es\n",
            (abip_float)p->tot_cg_its / (info->admm_iter + 1),
            p->total_solve_time / (info->admm_iter + 1) / 1e3);
    p->tot_cg_its = 0;
    p->total_solve_time = 0;
    return str;
}

/* M = inv(diag(RHO_Y * I + AA')) */
static void get_preconditioner
 (
  const ABIPMatrix *A,
  const ABIPSettings *stgs,
  ABIPLinSysWork *p
  )
{
    abip_int i;
    abip_int j;
    abip_int c1;
    abip_int c2;
    
    abip_float wrk;
    abip_float *M = p->M;
    
    memset(M, 0, A->m * sizeof(abip_float));
    
#if EXTRA_VERBOSE > 0
    abip_printf("getting pre-conditioner\n");
#endif
    
    for (i = 0; i < A->n; ++i)
    {
        c1 = A->p[i];
        c2 = A->p[i + 1];
        for (j = c1; j < c2; ++j)
        {
            wrk = A->x[j];
            M[A->i[j]] += wrk * wrk;
        }
    }
    
    for (i = 0; i < A->m; ++i)
    {
        M[i] = 1 / M[i];
        /* M[i] = 1; */
    }
    
#if EXTRA_VERBOSE > 0
    ABIP(print_array)(M, A->m, "M");
    abip_printf("norm M = %4f\n", ABIP(norm)(M, A->m));
    abip_printf("finished getting pre-conditioner\n");
#endif
}

static void transpose
 (
  const ABIPMatrix *A,
  ABIPLinSysWork *p
  )
{
    abip_int *Ci = p->At->i;
    abip_int *Cp = p->At->p;
    abip_float *Cx = p->At->x;
    
    abip_int m = A->m;
    abip_int n = A->n;
    
    abip_int *Ap = A->p;
    abip_int *Ai = A->i;
    abip_float *Ax = A->x;
    
    abip_int i;
    abip_int j;
    abip_int q;
    abip_int *z;
    abip_int c1;
    abip_int c2;
    
#if EXTRA_VERBOSE > 0
    ABIP(timer) transpose_timer;
    abip_printf("transposing A\n");
    ABIP(tic)(&transpose_timer);
#endif
    
    z = (abip_int *)abip_calloc(m, sizeof(abip_int));
    
    for (i = 0; i < Ap[n]; i++)
    {
        z[Ai[i]]++;           /* row counts */
    }
    
    ABIP(cumsum)(Cp, z, m);      /* row pointers */
    
    for (j = 0; j < n; j++)
    {
        c1 = Ap[j];
        c2 = Ap[j + 1];
        for (i = c1; i < c2; i++)
        {
            q = z[Ai[i]];
            Ci[q] = j;          /* place A(i,j) as entry C(j,i) */
            Cx[q] = Ax[i];
            z[Ai[i]]++;
        }
    }
    
    abip_free(z);
    
#if EXTRA_VERBOSE > 0
    abip_printf("finished transposing A, time: %1.2es\n",
                ABIP(tocq)(&transpose_timer) / 1e3);
#endif
}

void ABIP(free_lin_sys_work)
(
 ABIPLinSysWork *p
 )
{
    if (p)
    {
        if (p->p)
        {
            abip_free(p->p);
        }
        
        if (p->r)
        {
            abip_free(p->r);
        }
        
        if (p->Gp)
        {
            abip_free(p->Gp);
        }
        
        if (p->tmp)
        {
            abip_free(p->tmp);
        }
        
        if (p->At)
        {
            ABIPMatrix *At = p->At;
            
            if (At->i)
            {
                abip_free(At->i);
            }
            
            if (At->x)
            {
                abip_free(At->x);
            }
            
            if (At->p)
            {
                abip_free(At->p);
            }
            
            abip_free(At);
        }
        
        if (p->z)
        {
            abip_free(p->z);
        }
        
        if (p->M)
        {
            abip_free(p->M);
        }
        
        abip_free(p);
    }
}

/*y = (RHO_Y * I + AA')x */
static void mat_vec
 (
  const ABIPMatrix *A,
  const ABIPSettings *s,
  ABIPLinSysWork *p,
  const abip_float *x,
  abip_float *y
  )
{
    abip_float *tmp = p->tmp;
    memset(tmp, 0, A->n * sizeof(abip_float));
    ABIP(accum_by_Atrans)(A, p, x, tmp);
    memset(y, 0, A->m * sizeof(abip_float));
    ABIP(accum_by_A)(A, p, tmp, y);
    ABIP(add_scaled_array)(y, x, A->m, s->rho_y);
}

void ABIP(accum_by_Atrans)
(
 const ABIPMatrix *A,
 ABIPLinSysWork *p,
 const abip_float *x,
 abip_float *y
 )
{
    ABIP(_accum_by_Atrans)(A->n, A->x, A->i, A->p, x, y);
}

void ABIP(accum_by_A)
(
 const ABIPMatrix *A,
 ABIPLinSysWork *p,
 const abip_float *x,
 abip_float *y
 )
{
    ABIP(_accum_by_Atrans)(p->At->n, p->At->x, p->At->i, p->At->p, x, y);
}

void ABIP(normalize_A)
(
 ABIPMatrix *A,
 const ABIPSettings *stgs,
 ABIPScaling *scal)
{
    ABIP(_normalize_A)(A, stgs, scal);
}

void ABIP(un_normalize_A)
(
 ABIPMatrix *A,
 const ABIPSettings *stgs,
 const ABIPScaling *scal
 )
{
    ABIP(_un_normalize_A)(A, stgs, scal);
}

static void apply_pre_conditioner
 (
  abip_float *M,
  abip_float *z,
  abip_float *r,
  abip_int m,
  abip_float *ipzr
  )
{
    abip_int i;
    *ipzr = 0;
    
    for (i = 0; i < m; ++i)
    {
        z[i] = r[i] * M[i];
        *ipzr += z[i] * r[i];
    }
}

ABIPLinSysWork *ABIP(init_lin_sys_work)
(
 const ABIPMatrix *A,
 const ABIPSettings *stgs
 )
{
    ABIPLinSysWork *p = (ABIPLinSysWork *)abip_calloc(1, sizeof(ABIPLinSysWork));
    p->p = (abip_float *)abip_malloc((A->m) * sizeof(abip_float));
    p->r = (abip_float *)abip_malloc((A->m) * sizeof(abip_float));
    p->Gp = (abip_float *)abip_malloc((A->m) * sizeof(abip_float));
    p->tmp = (abip_float *)abip_malloc((A->n) * sizeof(abip_float));
    
    /* memory for A transpose */
    p->At = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
    p->At->m = A->n;
    p->At->n = A->m;
    p->At->i = (abip_int *)abip_malloc((A->p[A->n]) * sizeof(abip_int));
    p->At->p = (abip_int *)abip_malloc((A->m + 1) * sizeof(abip_int));
    p->At->x = (abip_float *)abip_malloc((A->p[A->n]) * sizeof(abip_float));
    transpose(A, p);
    
    /* preconditioner memory */
    p->z = (abip_float *)abip_malloc((A->m) * sizeof(abip_float));
    p->M = (abip_float *)abip_malloc((A->m) * sizeof(abip_float));
    get_preconditioner(A, stgs, p);
    
    p->total_solve_time = 0;
    p->tot_cg_its = 0;
    
    if (!p->p || !p->r || !p->Gp || !p->tmp || !p->At || !p->At->i || !p->At->p || !p->At->x)
    {
        ABIP(free_lin_sys_work)(p);
        return ABIP_NULL;
    }
    
    return p;
}

/* solves (I+AA')x = b, s warm start, solution stored in b */
static abip_int pcg
 (
  const ABIPMatrix *A,
  const ABIPSettings *stgs,
  ABIPLinSysWork *pr,
  const abip_float *s,
  abip_float *b,
  abip_int max_its,
  abip_float tol
  )
{
    abip_int i;
    abip_int m = A->m;
    
    abip_float ipzr;
    abip_float ipzr_old;
    abip_float alpha;
    
    abip_float *p = pr->p;   /* cg direction */
    abip_float *Gp = pr->Gp; /* updated CG direction */
    abip_float *r = pr->r;   /* cg residual */
    abip_float *z = pr->z;   /* for preconditioning */
    abip_float *M = pr->M;   /* inverse diagonal preconditioner */
    
    if (s == ABIP_NULL)
    {
        memcpy(r, b, m * sizeof(abip_float));
        memset(b, 0, m * sizeof(abip_float));
    }
    else
    {
        mat_vec(A, stgs, pr, s, r);
        ABIP(add_scaled_array)(r, b, m, -1);
        ABIP(scale_array)(r, -1, m);
        memcpy(b, s, m * sizeof(abip_float));
    }
    
    /* check to see if we need to run CG at all */
    if (ABIP(norm)(r, m) < MIN(tol, 1e-18))
    {
        return 0;
    }
    
    apply_pre_conditioner(M, z, r, m, &ipzr);
    memcpy(p, z, m * sizeof(abip_float));
    
    /* main loop */
    for (i = 0; i < max_its; ++i)
    {
        mat_vec(A, stgs, pr, p, Gp);
        alpha = ipzr / ABIP(dot)(p, Gp, m);
        ABIP(add_scaled_array)(b, p, m, alpha);
        ABIP(add_scaled_array)(r, Gp, m, -alpha);
        
        if (ABIP(norm)(r, m) < tol)
        {
#if EXTRA_VERBOSE > 0
            abip_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, ABIP(norm)(r, m), (long)i + 1);
#endif
            
            return i + 1;
        }
        
        ipzr_old = ipzr;
        apply_pre_conditioner(M, z, r, m, &ipzr);
        ABIP(scale_array)(p, ipzr / ipzr_old, m);
        ABIP(add_scaled_array)(p, z, m, 1);
    }
    
    return i;
}

abip_int ABIP(solve_lin_sys)
(       
 const ABIPMatrix *A,
 const ABIPSettings *stgs,
 ABIPLinSysWork *p,
 abip_float *b,
 const abip_float *s,
 abip_int iter
 ) 
{
    abip_int cg_its;
    ABIP(timer) linsys_timer;
    
    abip_float cg_tol = ABIP(norm)(b, A->m) * (iter < 0 ? CG_BEST_TOL
                                               : CG_MIN_TOL / POWF((abip_float)iter + 1, stgs->cg_rate));
    
    cg_tol = MAX(cg_tol, 1e-07);
    
    ABIP(tic)(&linsys_timer);
    
    /* solves Mx = b, for x but stores result in b */
    /* s contains warm-start (if available) */
    ABIP(accum_by_A)(A, p, &(b[A->m]), b);
    
    /* solves (I+AA')x = b, s warm start, solution stored in b */
    cg_its = pcg(A, stgs, p, s, b, A->m, MAX(cg_tol, CG_BEST_TOL));
    ABIP(scale_array)(&(b[A->m]), -1, A->n);
    ABIP(accum_by_Atrans)(A, p, b, &(b[A->m]));
    
    if (iter >= 0)
    {
        p->tot_cg_its += cg_its;
    }
    
    p->total_solve_time += ABIP(tocq)(&linsys_timer);
    
#if EXTRA_VERBOSE > 0
    abip_printf("linsys solve time: %1.2es\n", ABIP(tocq)(&linsys_timer) / 1e3);
#endif
    
    return 0;
}
