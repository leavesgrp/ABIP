#include "direct.h"

char *ABIP(get_lin_sys_method) ( const ABIPMatrix *A, const ABIPSettings *stgs ) {
    char *tmp = (char *) abip_malloc(sizeof(char) * 128);
#ifdef ABIP_PARDISO
    sprintf(tmp, "sparse-direct-intel-pardiso, nnz in A = %li", (long)A->p[A->n]);
#else
    sprintf(tmp, "sparse-direct, nnz in A = %li", (long)A->p[A->n]);
#endif
    return tmp;
}

char *ABIP(get_lin_sys_summary) ( ABIPLinSysWork *p, const ABIPInfo *info ) {
    char *str = (char *) abip_malloc(sizeof(char) * 128);
#ifdef ABIP_PARDISO
    sprintf(str, "Linear system avg solve time: %1.2es\n", p->total_solve_time / (info->admm_iter + 1) / 1e3);
#else
    abip_int n = p->L->n;
    sprintf(str, "\tLin-sys: nnz in L factor: %li, avg solve time: %1.2es\n",
            (long)(p->L->p[n] + n), p->total_solve_time / (info->admm_iter + 1) / 1e3);
#endif
    
    p->total_solve_time = 0; return str;
}

void ABIP(free_lin_sys_work) ( ABIPLinSysWork *p ) {
    if (p) {
        if (p->L)  { ABIP(cs_spfree)(p->L); }
        if (p->P)  { abip_free(p->P); }
        if (p->bp) { abip_free(p->bp); }
        if (p->D)  { abip_free(p->D); }
        abip_free(p);
    }
}

#ifdef ABIP_PARDISO
void ABIP(free_lin_sys_work_pds) ( ABIPLinSysWork *p, ABIPMatrix *A ) {
    if (p) {
        pardisoFree(p, A);
        if (p->P)  { abip_free(p->D); }
        if (p->D)  { abip_free(p->D); }
        abip_free(p);
    }
}
#endif

cs *form_kkt ( const ABIPMatrix *A, const ABIPSettings *s ) {
    
    abip_int j, k, kk; cs *K_cs;
    const abip_int Annz = A->p[A->n];
    const abip_int Knnzmax = A->m + A->n + Annz;
    
    cs *K = ABIP(cs_spalloc)(A->m + A->n, A->m + A->n, Knnzmax, 1, 1);
    
#if EXTRA_VERBOSE > 0
    abip_printf("forming KKT\n");
#endif
    
    if (!K)
    {
        return ABIP_NULL;
    }
    
    kk = 0;
    for (k = 0; k < A->m; k++)
    {
        K->i[kk] = k;
        K->p[kk] = k;
        K->x[kk] = s->rho_y;
        kk++;
    }
    
    for (j = 0; j < A->n; j++)
    {
        for (k = A->p[j]; k < A->p[j + 1]; k++)
        {
            K->p[kk] = j + A->m;
            K->i[kk] = A->i[k];
            K->x[kk] = A->x[k];
            kk++;
        }
    }
    for (k = 0; k < A->n; k++)
    {
        K->i[kk] = k + A->m;
        K->p[kk] = k + A->m;
        K->x[kk] = -1;
        kk++;
    }
    
    K->nnz = Knnzmax;
    K_cs = ABIP(cs_compress)(K);
    
#ifdef ABIP_PARDISO
    cs *KT = ABIP(cs_transpose)(K_cs, 1);
    ABIP(cs_spfree)(K); ABIP(cs_spfree)(K_cs);
    return (KT);
#else
    ABIP(cs_spfree)(K);
    return (K_cs);
#endif
}

abip_int _ldl_init ( cs *A, abip_int P[], abip_float **info ) {
#ifdef ABIP_PARDISO
    return 0;
#else
    *info = (abip_float *)abip_malloc(AMD_INFO * sizeof(abip_float));
    
#ifdef DLONG
    return (amd_l_order(A->n, A->p, A->i, P, (abip_float *) ABIP_NULL, *info));
#else
    return (amd_order(A->n, A->p, A->i, P, (abip_float *) ABIP_NULL, *info));
#endif
    
#endif
}

abip_int _ldl_factor ( cs *A, abip_int P[], abip_int Pinv[], cs **L, abip_float **D ) {
    
#ifdef ABIP_PARDISO
    return 0;
#else
    abip_int kk;
    abip_int n = A->n;
    
    abip_int *Parent = (abip_int *)abip_malloc(n * sizeof(abip_int));
    abip_int *Lnnz = (abip_int *)abip_malloc(n * sizeof(abip_int));
    abip_int *Flag = (abip_int *)abip_malloc(n * sizeof(abip_int));
    abip_int *Pattern = (abip_int *)abip_malloc(n * sizeof(abip_int));
    abip_float *Y = (abip_float *)abip_malloc(n * sizeof(abip_float));
    (*L)->p = (abip_int *)abip_malloc((1 + n) * sizeof(abip_int));
    
    /*abip_int Parent[n], Lnz[n], Flag[n], Pattern[n]; */
    /*abip_float Y[n]; */
    
    LDL_symbolic(n, A->p, A->i, (*L)->p, Parent, Lnnz, Flag, P, Pinv);
    
    (*L)->nnzmax = *((*L)->p + n);
    (*L)->x = (abip_float *)abip_malloc((*L)->nnzmax * sizeof(abip_float));
    (*L)->i = (abip_int *)abip_malloc((*L)->nnzmax * sizeof(abip_int));
    *D = (abip_float *)abip_malloc(n * sizeof(abip_float));
    
    if (!(*D) || !(*L)->i || !(*L)->x || !Y || !Pattern || !Flag || !Lnnz || !Parent)
    {
        abip_free(Pattern); abip_free(Y);
        return -1;
    }
    
#if EXTRA_VERBOSE > 0
    abip_printf("numeric factorization\n");
#endif
    
    kk = LDL_numeric(n, A->p, A->i, A->x, (*L)->p, Parent, Lnnz, (*L)->i, (*L)->x, *D, Y, Pattern, Flag, P, Pinv);
    
#if EXTRA_VERBOSE > 0
    abip_printf("finished numeric factorization\n");
#endif
    
    abip_free(Parent);
    abip_free(Lnnz);
    abip_free(Flag);
    abip_free(Pattern);
    abip_free(Y);
    return (kk - n);
#endif
    
}

void _ldl_solve ( abip_float *x, abip_float b[], cs *L, abip_float D[],
                  abip_int P[], abip_float *bp ) {
#ifdef ABIP_PARDISO
    return;
#else
    abip_int n = L->n;
    if (P == ABIP_NULL)
    {
        if (x != b)
        {
            memcpy(x, b, n * sizeof(abip_float));
        }
        
        LDL_lsolve(n, x, L->p, L->i, L->x);
        LDL_dsolve(n, x, D);
        LDL_ltsolve(n, x, L->p, L->i, L->x);
    }
    else
    {
        LDL_perm(n, bp, b, P);
        LDL_lsolve(n, bp, L->p, L->i, L->x);
        LDL_dsolve(n, bp, D);
        LDL_ltsolve(n, bp, L->p, L->i, L->x);
        LDL_permt(n, x, bp, P);
    }
#endif
}

void ABIP(accum_by_Atrans) ( const ABIPMatrix *A, ABIPLinSysWork *p,
                             const abip_float *x, abip_float *y ) {
    ABIP(_accum_by_Atrans)(A->n, A->x, A->i, A->p, x, y);
}

void ABIP(accum_by_A) ( const ABIPMatrix *A, ABIPLinSysWork *p,
                        const abip_float *x, abip_float *y ) {
    ABIP(_accum_by_A)(A->n, A->x, A->i, A->p, x, y);
}

void ABIP(normalize_A) ( ABIPMatrix *A, const ABIPSettings *stgs, ABIPScaling *scal ) {
    ABIP(_normalize_A)(A, stgs, scal);
}

void ABIP(un_normalize_A) ( ABIPMatrix *A, const ABIPSettings *stgs, const ABIPScaling *scal ) {
    ABIP(_un_normalize_A)(A, stgs, scal);
}

abip_int factorize ( const ABIPMatrix *A, const ABIPSettings *stgs, ABIPLinSysWork *p ) {

#ifdef ABIP_PARDISO
    abip_int ret_code = 0;
    cs *K = form_kkt(A, stgs);
    if (!K) { return -1; }
    ret_code = pardisoFactorize(p, K);
    ABIP(cs_spfree)(K);
    return ret_code;
#else
    abip_float *info;
    abip_int *Pinv;
    abip_int amd_status;
    abip_int ldl_status;
    
    cs *C;
    cs *K = form_kkt(A, stgs);
    if (!K)
    {
        return -1;
    }
    
    amd_status = _ldl_init(K, p->P, &info);
    if (amd_status < 0)
    {
        return (amd_status);
    }
    
#if EXTRA_VERBOSE > 0
    if (stgs->verbose)
    {
        abip_printf("Matrix factorization info:\n");
        
#ifdef DLONG
        amd_l_info(info);
#else
        amd_info(info);
#endif
    }
#endif
    
    Pinv = ABIP(cs_pinv)(p->P, A->m + A->n);
    C = ABIP(cs_symperm)(K, Pinv, 1);
    ldl_status = _ldl_factor(C, ABIP_NULL, ABIP_NULL, &p->L, &p->D);
    
    ABIP(cs_spfree)(C);
    ABIP(cs_spfree)(K);
    abip_free(Pinv);
    abip_free(info);
    
    return (ldl_status);
#endif
}


ABIPLinSysWork *ABIP(init_lin_sys_work)
(
 const ABIPMatrix *A,
 const ABIPSettings *stgs
 )
{
    ABIPLinSysWork *p = (ABIPLinSysWork *) abip_calloc(1, sizeof(ABIPLinSysWork));
    abip_int m_plus_n = A->m + A->n;
#ifdef ABIP_PARDISO
    p->P = (abip_int *) abip_malloc(sizeof(abip_int) * m_plus_n);
    p->D = (abip_float *) abip_malloc(sizeof(abip_float) * m_plus_n);
    if (factorize(A, stgs, p)) {
        ABIP(free_lin_sys_work)(p); return ABIP_NULL;
    }
#else
    p->P = (abip_int *) abip_malloc(sizeof(abip_int) * m_plus_n);
    p->L = (cs *) abip_malloc(sizeof(cs));
    p->bp = (abip_float *) abip_malloc(m_plus_n * sizeof(abip_float));
    p->L->m = m_plus_n;
    p->L->n = m_plus_n;
    p->L->nnz = -1;
    
    if (factorize(A, stgs, p) < 0)
    {
        ABIP(free_lin_sys_work)(p);
        return ABIP_NULL;
    }
#endif
    p->total_solve_time = 0.0;
    return p;
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
    ABIP(timer) linsys_timer;
    ABIP(tic)(&linsys_timer);
#ifdef ABIP_PARDISO
    pardisoSolve(p, (ABIPMatrix *) A, b);
#else
    _ldl_solve(b, b, p->L, p->D, p->P, p->bp);
#endif
    p->total_solve_time += ABIP(tocq)(&linsys_timer);
#if EXTRA_VERBOSE > 0
    abip_printf("linsys solve time: %1.2es\n", ABIP(tocq)(&linsys_timer) / 1e3);
#endif
    
    return 0;
}
