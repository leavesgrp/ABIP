#define _CRT_SECURE_NO_WARNINGS
#include "linsys.h"
#include "amd.h"

#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)

abip_int ABIP(copy_A_matrix)
(
    ABIPMatrix** dstp,
    const ABIPMatrix* src
    )
{
    abip_int Annz = src->p[src->n];
    ABIPMatrix* A = (ABIPMatrix*)abip_calloc(1, sizeof(ABIPMatrix));
    if (!A)
    {
        return 0;
    }
    A->n = src->n;
    A->m = src->m;
    A->x = (abip_float*)abip_malloc(sizeof(abip_float) * Annz);
    A->i = (abip_int*)abip_malloc(sizeof(abip_int) * Annz);
    A->p = (abip_int*)abip_malloc(sizeof(abip_int) * (src->n + 1));

    if (!A->x || !A->i || !A->p)
    {
        return 0;
    }

    memcpy(A->x, src->x, sizeof(abip_float) * Annz);
    memcpy(A->i, src->i, sizeof(abip_int) * Annz);
    memcpy(A->p, src->p, sizeof(abip_int) * (src->n + 1));

    *dstp = A;
    return 1;
}

char* ABIP(get_lin_sys_method)
(
    spe_problem *spe
 )
{
    char* tmp = (char*)abip_malloc(sizeof(char) * 128);

    if(spe->data->A == ABIP_NULL){
        sprintf(tmp, "This problem has no linear constraints");
        return tmp;
    }

    abip_int n = spe->data->A->p[spe->data->A->n];

    if (spe->stgs->linsys_solver == 0){
        sprintf(tmp, "sparse-direct using MKL-DSS, nnz in A = %li", (long)n);
    }
    else if (spe->stgs->linsys_solver == 1) {
        sprintf(tmp, "sparse-direct using QDLDL, nnz in A = %li", (long)n);
    }
    else if (spe->stgs->linsys_solver == 2) {
        sprintf(tmp, "sparse-direct using sparse cholesky, nnz in A = %li", (long)n);
    }
    else if (spe->stgs->linsys_solver == 3) {
        sprintf(tmp, "sparse-indirect using pcg, nnz in A = %li", (long)n);
    }
    else if (spe->stgs->linsys_solver == 4){
        sprintf(tmp, "sparse-direct using MKL-PARDISO, nnz in A = %li", (long)n);
    }
    else if (spe->stgs->linsys_solver == 5){
        sprintf(tmp, "dense-direct using dense cholesky, nnz in A = %li", (long)n);
    }
    else{
        sprintf(tmp,"\nlinsys solver type error\n");

    }
    return tmp;
}

char* ABIP(get_lin_sys_summary)
(
    spe_problem *self,
    ABIPInfo* info
 )
{
    char* str = (char*)abip_malloc(sizeof(char) * 128);

    
    abip_int n = self->L->nnz_LDL;
    info->avg_linsys_time = self->L->total_solve_time / (info->admm_iter + 1) / 1e3;


    if(self->stgs->linsys_solver == 3) {
        sprintf(str, "\tLin-sys: avg # CG iterations: %2.2f, avg solve time per admm iter: %1.2es\n",
            (abip_float)self->L->total_cg_iters / (info->admm_iter + 1), info->avg_linsys_time);
    }
    else{
        sprintf(str, "\tLin-sys: nnz in L factor: %li, avg solve time per admm iter: %1.2es\n",
        (long)n, info->avg_linsys_time);
    }
    
    info->avg_cg_iters = (abip_float)self->L->total_cg_iters / (info->admm_iter + 1);
    self->L->total_solve_time = 0;
    self->L->total_cg_iters = 0;
    
    return str;
}


abip_int ABIP(validate_lin_sys)
(
    const ABIPMatrix* A
    )
{
    abip_int i;
    abip_int r_max;
    abip_int Annz;

    if(A == ABIP_NULL){
        return 0;
    }

    if (!A->x || !A->i || !A->p)
    {
        abip_printf("ERROR: incomplete data!\n");
        return -1;
    }

    for (i = 0; i < A->n; ++i)
    {
        if (A->p[i] == A->p[i + 1])
        {
            abip_printf("WARN: the %li-th column empty!\n", (long)i);
        }
        else if (A->p[i] > A->p[i + 1])
        {
            abip_printf("ERROR: the column pointers decreases!\n");
            return -1;
        }
    }

    Annz = A->p[A->n];
    if (((abip_float)Annz / A->m > A->n) || (Annz <= 0))
    {
        abip_printf("ERROR: the number of nonzeros in A = %li, outside of valid range!\n", (long)Annz);
        return -1;
    }

    r_max = 0;
    for (i = 0; i < Annz; ++i)
    {
        if (A->i[i] > r_max)
        {
            r_max = A->i[i];
        }
    }
    if (r_max > A->m - 1)
    {
        abip_printf("ERROR: the number of rows in A is inconsistent with input dimension!\n");
        return -1;
    }

    return 0;
}

void ABIP(free_A_matrix)
(
    ABIPMatrix* A
    )
{
    if (A->x)
    {
        abip_free(A->x);
    }
    if (A->i)
    {
        abip_free(A->i);
    }
    if (A->p)
    {
        abip_free(A->p);
    }

    abip_free(A);
}

#if EXTRA_VERBOSE > 0

static void print_A_matrix
(
    const ABIPMatrix* A
)
{
    abip_int i;
    abip_int j;

    if (A->p[A->n] < 2500)
    {
        abip_printf("\n");
        for (i = 0; i < A->n; ++i)
        {
            abip_printf("Col %li: ", (long)i);
            for (j = A->p[i]; j < A->p[i + 1]; j++)
            {
                abip_printf("A[%li,%li] = %4f, ", (long)A->i[j], (long)i, A->x[j]);
            }
            abip_printf("norm col = %4f\n", ABIP(norm)(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
        }
        abip_printf("norm A = %4f\n", ABIP(norm)(A->x, A->p[A->n]));
    }
}
#endif

//y += A'*x
void ABIP(accum_by_Atrans)
(
    const ABIPMatrix* A,
    const abip_float* x,
    abip_float* y
    )
{
    abip_int p;
    abip_int j;

    abip_int c1;
    abip_int c2;
    abip_float yj;

#if EXTRA_VERBOSE > 0
    ABIP(timer) mult_by_Atrans_timer;
    ABIP(tic)(&mult_by_Atrans_timer);
#endif

#ifdef _OPENMP
#pragma omp parallel for private(p, c1, c2, yj)
#endif

    for (j = 0; j < A->n; j++)
    {
        yj = y[j];
        c1 = A->p[j];
        c2 = A->p[j + 1];
        for (p = c1; p < c2; p++)
        {
            yj += A->x[p] * x[A->i[p]];
        }
        y[j] = yj;
    }

#if EXTRA_VERBOSE > 0
    abip_printf("mult By A trans time: %1.2e seconds. \n", ABIP(tocq)(&mult_by_Atrans_timer) / 1e3);
#endif
}

//y += A*x
void ABIP(accum_by_A)
(
    const ABIPMatrix* A,
    const abip_float* x,
    abip_float* y
    )
{
    abip_int p;
    abip_int j;

    abip_int c1;
    abip_int c2;
    abip_float xj;

#if EXTRA_VERBOSE > 0
    ABIP(timer) mult_by_A_timer;
    ABIP(tic)(&mult_by_A_timer);
#endif

#ifdef _OPENMP
#pragma omp parallel for private(p, c1, c2, xj)
    for (j = 0; j < n; j++)
    {
        xj = x[j];
        c1 = A->p[j];
        c2 = A->p[j + 1];
        for (p = c1; p < c2; p++)
        {
#pragma omp atomic
            y[A->i[p]] += A->x[p] * xj;
        }
    }
#endif

    for (j = 0; j < A->n; j++)
    {
        xj = x[j];
        c1 = A->p[j];
        c2 = A->p[j + 1];
        for (p = c1; p < c2; p++)
        {
            y[A->i[p]] += A->x[p] * xj;
        }
    }

#if EXTRA_VERBOSE > 0
    abip_printf("mult By A time: %1.2e seconds \n", ABIP(tocq)(&mult_by_A_timer) / 1e3);
#endif
}


static abip_int _ldl_init(cs *A, abip_int *P, abip_float **info) {
  *info = (abip_float *)abip_calloc(AMD_INFO, sizeof(abip_float));
  return amd_order(A->n, A->p, A->i, P, (abip_float *)ABIP_NULL, *info);
}

cs *permute_kkt(spe_problem *spe) {
  abip_float *info;
  abip_int *Pinv, amd_status, *idx_mapping, i;
  cs *kkt = spe->L->K;
  cs *kkt_perm;
  if (!kkt) {
    return ABIP_NULL;
  }
  amd_status = _ldl_init(kkt, spe->L->P, &info);
  if (amd_status < 0) {
    abip_printf("AMD permutatation error.\n");
    return ABIP_NULL;
  }
#if VERBOSITY > 0
  abip_printf("Matrix factorization info:\n");
  amd_info(info);
#endif
  Pinv = cs_pinv(spe->L->P, spe->L->K->n);
  kkt_perm = cs_symperm(kkt, Pinv, 1);
  abip_free(Pinv);
  abip_free(info);
  return kkt_perm;
}

static void _ldl_perm(abip_int n, abip_float *x, abip_float *b, abip_int *P) {
  abip_int j;
  for (j = 0; j < n; j++)
    x[j] = b[P[j]];
}

static void _ldl_permt(abip_int n, abip_float *x, abip_float *b, abip_int *P) {
  abip_int j;
  for (j = 0; j < n; j++)
    x[P[j]] = b[j];
}

static void _ldl_solve(abip_float *b, cs *L, abip_float *Dinv, abip_int *P, abip_float *bp) {
  /* solves PLDL'P' x = b for x */
  abip_int n = L->n;
  _ldl_perm(n, bp, b, P);
  QDLDL_solve(n, L->p, L->i, L->x, Dinv, bp);
  _ldl_permt(n, b, bp, P);
}


_MKL_DSS_HANDLE_t init_mkl_work(cs* K)
{
    _INTEGER_t error;
    MKL_INT create_opt = MKL_DSS_ZERO_BASED_INDEXING;
    MKL_INT order_opt = MKL_DSS_DEFAULTS;
    MKL_INT sym = MKL_DSS_SYMMETRIC;
    MKL_INT type = MKL_DSS_INDEFINITE;

    _MKL_DSS_HANDLE_t handle;

    error = dss_create(handle, create_opt);
    if (error != MKL_DSS_SUCCESS)
        goto printError;


    error = dss_define_structure(handle, sym, K->p, K->m, K->n, K->i, K->p[K->n]);
    if (error != MKL_DSS_SUCCESS)
        goto printError;

    
    error = dss_reorder(handle, order_opt, 0);
    if (error != MKL_DSS_SUCCESS)
        goto printError;

    error = dss_factor_real(handle, type, K->x);
    if (error != MKL_DSS_SUCCESS)
        goto printError;

    return handle;
printError:
    printf("MKL-DSS returned error code %i\n", error);
    return -1;
}

abip_int mkl_solve_linsys(_MKL_DSS_HANDLE_t handle, abip_float *b, abip_int n)
{
    _INTEGER_t error;
    abip_int nrhs = 1;
    abip_float* solValues = (abip_float*)abip_malloc(n * sizeof(abip_float));
    MKL_INT opt = MKL_DSS_DEFAULTS;
    error = dss_solve_real(handle, opt, b, nrhs, solValues);
    if (error != MKL_DSS_SUCCESS) {
        printf("solve err");
        exit(1);
    }

    memcpy(b, solValues, n * sizeof(abip_float));
    abip_free(solValues);
    return 0;
}


abip_int init_pardiso(spe_problem *self)
{
    
    MKL_INT PARDISO_PARAMS_LDL[64] = {
    
    1, /* Non-default value */ 3, /* P Nested dissection */ 0, /* Reserved          */
    0, /* No CG             */ 0, /* No user permutation */ 0, /* Overwriting    */
    0, /* Refinement report */ 0, /* Auto ItRef step     */ 0, /* Reserved          */
    8, /* Perturb           */ 0, /* Disable scaling     */ 0, /* No transpose      */
    0, /* Disable matching  */ 0, /* Report on pivots    */ 0, /* Output            */
    0, /* Output            */ 0, /* Output              */-1, /* No report         */
    0, /* No report         */ 0, /* Output              */ 2, /* Pivoting          */
    0, /* nPosEigVals       */ 0, /* nNegEigVals         */ 0, /* Classic factorize */
    0,                         0,                           0, /* Matrix checker    */
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         1, /* 0-based solve       */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0, /* No diagonal         */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0
    };

    
    cs *K = self->L->K;
    MKL_INT n = K->n;
    MKL_INT *ia = K->i;
    MKL_INT *ja = K->p;
    abip_float *a = K->x;

    self->L->mtype = -2;       /* Real symmetric matrix */
    /* RHS and solution vectors. */
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    //void *pt[64];
    /* Pardiso control parameters. */
    //MKL_INT iparm[64];
    MKL_INT phase;
    /* Auxiliary variables. */
    MKL_INT i;

/* -------------------------------------*/
/* .. Setup Pardiso control parameters. */
/* -------------------------------------*/
    for ( i = 0; i < 64; i++ )
    {
      self->L->iparm[i] = 0;
    }

    for ( i = 0; i < 64; i++ )
    {
        self->L->iparm[i] = PARDISO_PARAMS_LDL[i];
    }



    self->L->maxfct = 1;           /* Maximum number of numerical factorizations. */
    self->L->mnum = 1;         /* Which factorization to use. */
    self->L->msglvl = 1;           /* Print statistical information in file */
    self->L->error = 0;            /* Initialize error flag */
/* ----------------------------------------------------------------*/
/* .. Initialize the internal solver memory pointer. This is only  */
/*   necessary for the FIRST call of the PARDISO solver.           */
/* ----------------------------------------------------------------*/
    for ( i = 0; i < 64; i++ )
    {
        self->L->pt[i] = 0;
    }
/* --------------------------------------------------------------------*/
/* .. Reordering and Symbolic Factorization. This step also allocates  */
/*    all memory that is necessary for the factorization.              */
/* --------------------------------------------------------------------*/
    phase = 11;
    PARDISO (self->L->pt, &(self->L->maxfct), &(self->L->mnum), &(self->L->mtype), &phase,
             &n, a, ja, ia, &(self->L->idum), &nrhs, self->L->iparm, &(self->L->msglvl), &(self->L->ddum), &(self->L->ddum), &(self->L->error));
    if ( self->L->error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %i" , self->L->error);
        exit (1);
    }
    printf ("\nReordering completed ... ");
    printf ("\nNumber of nonzeros in factors = %i", self->L->iparm[17]);
    printf ("\nNumber of factorization MFLOPS = %i", self->L->iparm[18]);
/* ----------------------------*/
/* .. Numerical factorization. */
/* ----------------------------*/
    phase = 22;
    PARDISO (self->L->pt, &(self->L->maxfct), &(self->L->mnum), &(self->L->mtype), &phase,
             &n, a, ja, ia, &(self->L->idum), &nrhs, self->L->iparm, &(self->L->msglvl), &(self->L->ddum), &(self->L->ddum), &(self->L->error));
    if ( self->L->error != 0 )
    {
        printf ("\nERROR during numerical factorization: %i", self->L->error);
        exit (2);
    }
    printf ("\nFactorization completed ... ");

    return 0;
}


abip_int pardiso_solve(spe_problem *self, abip_float *b, abip_int n)
{
/* -----------------------------------------------*/
/* .. Back substitution and iterative refinement. */
/* -----------------------------------------------*/
    cs *K = self->L->K;
    MKL_INT phase = 33;
    self->L->iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
    /* Set right hand side to one. */
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    abip_float *x = (abip_float*)abip_malloc(n*sizeof(abip_float));

    PARDISO (self->L->pt, &self->L->maxfct, &self->L->mnum, &self->L->mtype, &phase,
             &n, K->x, K->p, K->i, &self->L->idum, &nrhs, self->L->iparm, &self->L->msglvl, b, x, &self->L->error);
    if (self->L->error != 0 )
    {
        printf ("\nERROR during solution: %i", self->L->error);
        exit (3);
    }
    memcpy(b,x,n*sizeof(abip_float));
    abip_free(x);

    return 0;
}


abip_int pardiso_free(spe_problem *self)
{
    cs *K = self->L->K;
    MKL_INT phase = -1;
    MKL_INT nrhs = 1;
    MKL_INT n = K->n;
    PARDISO (self->L->pt, &self->L->maxfct, &self->L->mnum, &self->L->mtype, &phase,
             &n, &self->L->ddum, K->p, K->i, &self->L->idum, &nrhs,
             self->L->iparm, &self->L->msglvl, &self->L->ddum, &self->L->ddum, &self->L->error);

    cs_spfree(K);
    
    return 0;
}



abip_int LDL_factor(cs *A, cs **L, abip_float *Dinv){

  //data for elim tree calculation
  QDLDL_int *etree;
  QDLDL_int *Lnz;
  QDLDL_int  sumLnz;

  //working data for factorisation
  QDLDL_int   *iwork;
  QDLDL_bool  *bwork;
  QDLDL_float *fwork;
  


  /*--------------------------------
   * pre-factorisation memory allocations
   *---------------------------------*/

  //These can happen *before* the etree is calculated
  //since the sizes are not sparsity pattern specific

  //For the elimination tree
  etree = (QDLDL_int*)malloc(sizeof(QDLDL_int)*A->n);
  Lnz   = (QDLDL_int*)malloc(sizeof(QDLDL_int)*A->n);

  //For the L factors.   Li and Lx are sparsity dependent
  //so must be done after the etree is constructed
  
  QDLDL_float *D = (QDLDL_float*)malloc(sizeof(QDLDL_float)*A->n);
//   (*Dinv) = (QDLDL_float*)malloc(sizeof(QDLDL_float) * A->n);
  

  //Working memory.  Note that both the etree and factor
  //calls requires a working vector of QDLDL_int, with
  //the factor function requiring 3*An elements and the
  //etree only An elements.   Just allocate the larger
  //amount here and use it in both places
  iwork = (QDLDL_int*)malloc(sizeof(QDLDL_int)*(3*A->n));
  bwork = (QDLDL_bool*)malloc(sizeof(QDLDL_bool)*A->n);
  fwork = (QDLDL_float*)malloc(sizeof(QDLDL_float)*A->n);

  /*--------------------------------
   * elimination tree calculation
   *---------------------------------*/
  sumLnz = QDLDL_etree(A->n,A->p,A->i,iwork,Lnz,etree);

  if (sumLnz < 0){
    return sumLnz;//error
  }
  /*--------------------------------
   * LDL factorisation
   *---------------------------------*/

  //First allocate memory for Li and Lx
  (*L) = cs_spalloc(A->n, A->n, sumLnz, 1, 0);
//   L->p    = (QDLDL_int*)malloc(sizeof(QDLDL_int)*(A->n+1));
//   L->i    = (QDLDL_int*)malloc(sizeof(QDLDL_int)*sumLnz);
//   L->x    = (QDLDL_float*)malloc(sizeof(QDLDL_float)*sumLnz);
//   Dinv  = (QDLDL_float*)malloc(sizeof(QDLDL_float)*A->n);

  //now factor
  abip_int status = QDLDL_factor(A->n,A->p,A->i,A->x,(*L)->p,(*L)->i,(*L)->x,D,Dinv,Lnz,etree,bwork,iwork,fwork);
  
  free(D);
  free(etree);
  free(Lnz);
  free(iwork);
  free(bwork);
  free(fwork);

  return status;

}

abip_int abip_cholsol (spe_problem *self, abip_float *b, abip_int n)
{
    abip_float *x ;
    x = cs_malloc (n, sizeof (abip_float)) ;    /* get workspace */
    abip_int ok = (self->L->S && self->L->N && x) ;
    if (ok)
    {
        cs_ipvec (self->L->S->pinv, b, x, n) ;   /* x = P*b */
        cs_lsolve (self->L->N->L, x) ;           /* x = L\x */
        cs_ltsolve (self->L->N->L, x) ;          /* x = L'\x */
        cs_pvec (self->L->S->pinv, x, b, n) ;    /* b = P'*x */
    }
    cs_free (x) ;
    return (ok) ;
}


abip_int pcg(spe_problem *self, abip_float *b, abip_float *x, abip_float rho_x, abip_int max_iter, abip_float tol){
    /*
        x is used for warm start
        result overwrite b
    */

    abip_int m = self->p;
    abip_int n = self->q;

    abip_float* ATx = (abip_float*)abip_malloc(n * sizeof(abip_float));
    memset(ATx, 0, n * sizeof(abip_float));
    self->spe_AT_times(self,x,ATx);
    ABIP(scale_array)(ATx, -1, n);

    abip_float* r = (abip_float*)abip_malloc(m * sizeof(abip_float));
    memcpy(r, b, m * sizeof(abip_float));
    self->spe_A_times(self,ATx,r);
    ABIP(add_scaled_array)(r, x, m, -rho_x);

    abip_float* z = (abip_float*)abip_malloc(m * sizeof(abip_float));
    for (int k = 0; k < m; k++) {
        z[k] = self->L->M[k] * r[k];
    }

    abip_float* p = (abip_float*)abip_malloc(m * sizeof(abip_float));
    memcpy(p, z, m * sizeof(abip_float));

    abip_float ip = ABIP(dot)(r, z, m);

    abip_int i;
    abip_float* Ap = (abip_float*)abip_malloc(m * sizeof(abip_float));
    abip_float alpha, ipold, beta;

    memcpy(b, x, m*sizeof(abip_float));

    for (i = 0; i < max_iter; i++) {

        memset(ATx, 0, n * sizeof(abip_float));
        self->spe_AT_times(self,p,ATx);
        memcpy(Ap, p, m * sizeof(abip_float));
        ABIP(scale_array)(Ap, rho_x, m);
        self->spe_A_times(self,ATx,Ap);

        alpha = ip / (ABIP(dot)(Ap, p, m));

        // ABIP(add_scaled_array)(x, p, m, alpha);
        ABIP(add_scaled_array)(b, p, m, alpha);

        ABIP(add_scaled_array)(r, Ap, m, -alpha);

        if (ABIP(norm)(r, m) < tol) {
#if EXTRA_VERBOSE > 0
            abip_printf("CG took %d iterations to converge, resisdual %4f <= tolerance %4f\n", i, ABIP(norm)(r, m), tol);
#endif
            
            abip_free(ATx);
            abip_free(r);
            abip_free(z);
            abip_free(p);
            abip_free(Ap);
            
            return i + 1;
        }

        for (int j = 0; j < m; j++) {
            z[j] = self->L->M[j] * r[j];
        }

        ipold = ip;
        ip = ABIP(dot)(z, r, m);
        beta = ip / ipold;

        ABIP(scale_array)(p, beta, m);
        ABIP(add_scaled_array)(p, z, m, 1);
    }

    printf("CG did not converge within %d iterations, resisdual %4f > tolerance %4f\n", max_iter, ABIP(norm)(r, m), tol);

    abip_free(ATx);
    abip_free(r);
    abip_free(z);
    abip_free(p);
    abip_free(Ap);

    return i + 1;
}


/* y = (R_x + Q + A' R_y^{-1} A) x */
static void mat_vec(spe_problem *self, const abip_float *x, abip_float *y) {

    abip_int m = self->m;
    abip_int n = self->n;
    abip_int i;

    memcpy(y, x, n * sizeof(abip_float));
    
    for(i=0; i<n; i++){
        y[i] *= self->rho_dr[i + m];
    }

    if(self->Q != ABIP_NULL){
        ABIP(accum_by_A)(self->Q, x, y);
    }

    abip_float *tem = (abip_float*)abip_malloc(sizeof(abip_float) * m);
    memset(tem, 0, m * sizeof(abip_float));
    ABIP(accum_by_A)(self->A, x, tem);
    for(i=0; i<m; i++){
        tem[i] /= self->rho_dr[i];
    }

    ABIP(accum_by_Atrans)(self->A, tem, y);

    abip_free(tem);
}



abip_int qcp_pcg(spe_problem *self, abip_float *b, abip_float *x, abip_int max_iter, abip_float tol){

    /*
        x is used for warm start
        result overwrite b
    */

    abip_int m = self->m;
    abip_int n = self->n;
    abip_int i,j;

    abip_float ztr, ztr_prev, alpha;
    abip_float *p = (abip_float *)abip_calloc(n, sizeof(abip_float));   /* cg direction */
    abip_float *Gp = (abip_float *)abip_calloc(n, sizeof(abip_float)); /* updated CG direction */
    abip_float *r = (abip_float *)abip_calloc(n, sizeof(abip_float));   /* cg residual */
    abip_float *z = (abip_float *)abip_calloc(n, sizeof(abip_float));;   /* for preconditioning */

    if (x == ABIP_NULL) {
        /* no warm_start, take x = 0 */
        /* r = b */
        memcpy(r, b, n * sizeof(abip_float));
        /* b = 0 */
        memset(b, 0, n * sizeof(abip_float));
    } else {
        /* r = Mat * s */
        mat_vec(self, x, r);
        /* r = Mat * s - b */
        ABIP(add_scaled_array)(r, b, n, -1.);
        /* r = b - Mat * s */
        ABIP(scale_array)(r, -1., n);
        /* b = s */
        memcpy(b, x, n * sizeof(abip_float));
    }
    /* check to see if we need to run CG at all */
    if (ABIP(norm_inf)(r, n) < MAX(tol, 1e-12)) {

        abip_free(p);
        abip_free(Gp);
        abip_free(r);
        abip_free(z);
        return 0;
    }

    /* z = M r (M is inverse preconditioner) */
    memcpy(z, r, n * sizeof(abip_float));
    ABIP(c_dot)(z, self->L->M, n);

    /* ztr = z'r */
    ztr = ABIP(dot)(z, r, n);
    /* p = z */
    memcpy(p, z, n * sizeof(abip_float));

    for (i = 0; i < max_iter; ++i) {
        /* Gp = Mat * p */
        mat_vec(self, p, Gp);
        /* alpha = z'r / p'G p */
        alpha = ztr / ABIP(dot)(p, Gp, n);
        /* b += alpha * p */
        ABIP(add_scaled_array)(b, p, n, alpha);
        /* r -= alpha * G p */
        ABIP(add_scaled_array)(r, Gp, n, -alpha);

        if (ABIP(norm_inf)(r, n) < tol) {
            break;
        }
        /* z = M r (M is inverse preconditioner) */
        memcpy(z, r, n * sizeof(abip_float));
        ABIP(c_dot)(z, self->L->M, n);

        ztr_prev = ztr;
        /* ztr = z'r */
        ztr = ABIP(dot)(z, r, n);
        /* p = beta * p, where beta = ztr / ztr_prev */
        ABIP(scale_array)(p, ztr / ztr_prev, n);
        /* p = z + beta * p */
        ABIP(add_scaled_array)(p, z, n, 1.);
    }

    printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, ABIP(norm_inf)(r, n), (long)i + 1);
    if(i == max_iter - 1){
        printf("CG did not converge within %d iterations, resisdual %4f > tolerance %4f\n", max_iter, ABIP(norm)(r, m), tol);
    }

    abip_free(p);
    abip_free(Gp);
    abip_free(r);
    abip_free(z);

    return i + 1;

}



/* y = ( rhoy * I + A * H^-1 * A' ) x , where H = (rhox * I + Q)*/
static void svm_mat_vec(spe_problem *self, const abip_float *x, abip_float *y) {

    abip_int m = self->p;
    abip_int n = self->q;
    abip_int i;

    memcpy(y, x, m * sizeof(abip_float));
    
    for(i=0; i<m; i++){
        y[i] *= self->rho_dr[i];
    }

    abip_float *tem = (abip_float*)abip_malloc(sizeof(abip_float) * n);
    memset(tem, 0, n * sizeof(abip_float));
    self->spe_AT_times(self, x, tem);

    abip_float *H = (abip_float*)abip_malloc(sizeof(abip_float) * n);
    for(i=0; i<self->q; i++){
        if(i<self->n) H[i] = self->stgs->rho_x + self->Q->x[i];
        else H[i] = self->stgs->rho_x;
    }

    for(i=0; i<n; i++){
        // tem[i] /= self->H[i];
        tem[i] /= H[i];
    }

    self->spe_A_times(self, tem, y);

    abip_free(tem);
    abip_free(H);
}



abip_int svmqp_pcg(spe_problem *self, abip_float *b, abip_float *x, abip_int max_iter, abip_float tol){

    /*
        x is used for warm start
        result overwrite b
    */

    // abip_int m = self->p;
    abip_int n = self->p;
    abip_int i,j;

    abip_float ztr, ztr_prev, alpha;
    abip_float *p = (abip_float *)abip_calloc(n, sizeof(abip_float));   /* cg direction */
    abip_float *Gp = (abip_float *)abip_calloc(n, sizeof(abip_float)); /* updated CG direction */
    abip_float *r = (abip_float *)abip_calloc(n, sizeof(abip_float));   /* cg residual */
    abip_float *z = (abip_float *)abip_calloc(n, sizeof(abip_float));;   /* for preconditioning */

    if (x == ABIP_NULL) {
        /* no warm_start, take x = 0 */
        /* r = b */
        memcpy(r, b, n * sizeof(abip_float));
        /* b = 0 */
        memset(b, 0, n * sizeof(abip_float));
    } else {
        /* r = Mat * s */
        svm_mat_vec(self, x, r);
        /* r = Mat * s - b */
        ABIP(add_scaled_array)(r, b, n, -1.);
        /* r = b - Mat * s */
        ABIP(scale_array)(r, -1., n);
        /* b = s */
        memcpy(b, x, n * sizeof(abip_float));
    }
    /* check to see if we need to run CG at all */
    if (ABIP(norm_inf)(r, n) < MAX(tol, 1e-12)) {

        abip_free(p);
        abip_free(Gp);
        abip_free(r);
        abip_free(z);
        return 0;
    }

    /* z = M r (M is inverse preconditioner) */
    memcpy(z, r, n * sizeof(abip_float));
    ABIP(c_dot)(z, self->L->M, n);

    /* ztr = z'r */
    ztr = ABIP(dot)(z, r, n);
    /* p = z */
    memcpy(p, z, n * sizeof(abip_float));

    for (i = 0; i < max_iter; ++i) {
        /* Gp = Mat * p */
        svm_mat_vec(self, p, Gp);
        /* alpha = z'r / p'G p */
        alpha = ztr / ABIP(dot)(p, Gp, n);
        /* b += alpha * p */
        ABIP(add_scaled_array)(b, p, n, alpha);
        /* r -= alpha * G p */
        ABIP(add_scaled_array)(r, Gp, n, -alpha);

        if (ABIP(norm_inf)(r, n) < tol) {
            break;
        }
        /* z = M r (M is inverse preconditioner) */
        memcpy(z, r, n * sizeof(abip_float));
        ABIP(c_dot)(z, self->L->M, n);

        ztr_prev = ztr;
        /* ztr = z'r */
        ztr = ABIP(dot)(z, r, n);
        /* p = beta * p, where beta = ztr / ztr_prev */
        ABIP(scale_array)(p, ztr / ztr_prev, n);
        /* p = z + beta * p */
        ABIP(add_scaled_array)(p, z, n, 1.);
    }

    printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, ABIP(norm_inf)(r, n), (long)i + 1);
    if(i == max_iter - 1){
        printf("CG did not converge within %d iterations, resisdual %4f > tolerance %4f\n", max_iter, ABIP(norm)(r, n), tol);
    }

    abip_free(p);
    abip_free(Gp);
    abip_free(r);
    abip_free(z);

    return i + 1;

}



// MKL-LAPACK dense cholesky
abip_int init_dense_chol(spe_problem *spe){

    // K is CSC format and only stores upper triangle
    spe->L->U = ABIP(csc_to_dense)(spe->L->K, ColMajor);
    abip_int info;
    abip_int n = spe->L->K->n;

    info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U',  n, spe->L->U, n);
    return info; // 0 if successful
}


abip_int dense_chol_sol(spe_problem *spe, abip_float *b, abip_int n){

    abip_int info;

    info = LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'U', n, 1, spe->L->U, n, b, n);
    return info; // 0 if successful
}

abip_int dense_chol_free(spe_problem *spe){
    if(spe->L->U){
        abip_free(spe->L->U);
    }
    if(spe->L){
        abip_free(spe->L);
    }

    return 0;
}
/*-----------------------------------------------------------------*/



abip_int ABIP(init_linsys_work)(spe_problem *spe){

    if(spe->stgs->linsys_solver == 3){//pcg
        spe->L->S = ABIP_NULL ;              
        spe->L->N = ABIP_NULL; 
        spe->L->handle = ABIP_NULL;
        spe->L->Dinv = ABIP_NULL;
        spe->L->L = ABIP_NULL;
        spe->L->P = ABIP_NULL;
        spe->L->U = ABIP_NULL;
        

        spe->L->total_solve_time = 0.0;
        spe->L->total_cg_iters = 0;
        
        return 0;
    }

    cs *K = spe->L->K;
    spe->L->nnz_LDL =  K->nzmax; // K is NULL ptr if using pcg

    printf("\nStarting decomposition, with\nn = %d, m = %d , nnz = %d\n", K->m, K->n, K->nzmax);

    if(spe->stgs->linsys_solver == 0){//mkl_dss
        spe->L->handle = init_mkl_work(K);
        if(spe->L->handle == -1){
            printf("\nerror in LDL factorization using MKL-DSS\n");
            return -1;
        }
        cs_spfree(K);
        spe->L->Dinv = ABIP_NULL;                                                                                                                                                                                          
        spe->L->L = ABIP_NULL;
        spe->L->P = ABIP_NULL;
        spe->L->M = ABIP_NULL;
        spe->L->S = ABIP_NULL;
        spe->L->N = ABIP_NULL;
        spe->L->U = ABIP_NULL;

        spe->L->total_solve_time = 0.0;
        return 0;
    }

    else if(spe->stgs->linsys_solver == 1){//qdldl

        spe->L->Dinv = (abip_float*)abip_malloc(K->n*sizeof(abip_float));
        spe->L->bp = (abip_float*)abip_malloc(K->n*sizeof(abip_float));
        spe->L->P = (abip_int*)abip_malloc(K->n*sizeof(abip_int));
        cs *kkt_perm = permute_kkt(spe);
        cs_spfree(K);
        if(LDL_factor(kkt_perm, &(spe->L->L), spe->L->Dinv) < 0){
            spe->L->L = ABIP_NULL;
            printf("\nerror in LDL factorization using QDLDL\n");
            return -1;
        }
        cs_spfree(kkt_perm);

        spe->L->handle = ABIP_NULL;
        spe->L->M = ABIP_NULL;
        spe->L->S = ABIP_NULL;
        spe->L->N = ABIP_NULL;
        spe->L->U = ABIP_NULL;

        spe->L->total_solve_time = 0.0;
        return 0;

    }
    else if(spe->stgs->linsys_solver == 2){//cholesky
        spe->L->S = cs_schol (1, K) ;              
        spe->L->N = cs_chol (K, spe->L->S) ; 
        cs_spfree(K);
        spe->L->handle = ABIP_NULL;
        spe->L->Dinv = ABIP_NULL;
        spe->L->L = ABIP_NULL;
        spe->L->P = ABIP_NULL;
        spe->L->M = ABIP_NULL;
        spe->L->U = ABIP_NULL;

        spe->L->total_solve_time = 0.0;
        return 0;
    }
    else if(spe->stgs->linsys_solver == 4){//mkl_parsido
        init_pardiso(spe);
        spe->L->Dinv = ABIP_NULL;
        spe->L->handle = ABIP_NULL;
        spe->L->L = ABIP_NULL;
        spe->L->P = ABIP_NULL;
        spe->L->M = ABIP_NULL;
        spe->L->S = ABIP_NULL;
        spe->L->N = ABIP_NULL;
        spe->L->U = ABIP_NULL;

        spe->L->total_solve_time = 0.0;
        return 0;
    }
    else if(spe->stgs->linsys_solver == 5){//lapack dense chol

        init_dense_chol(spe);
        spe->L->handle = ABIP_NULL;
        spe->L->Dinv = ABIP_NULL;
        spe->L->L = ABIP_NULL;
        spe->L->P = ABIP_NULL;
        spe->L->M = ABIP_NULL;
        spe->L->S = ABIP_NULL;
        spe->L->N = ABIP_NULL;

        spe->L->total_solve_time = 0.0;
        return 0;
    }
 
    else{
        printf("\nlinsys solver type error\n");
        return -1;
    }
}


abip_int ABIP(solve_linsys)(spe_problem *spe, abip_float *b, abip_int n, abip_float *pcg_warm_start, abip_float pcg_tol){

    if(spe->stgs->linsys_solver == 0){//mkl_dss
        mkl_solve_linsys(spe->L->handle, b, n);
        return 0;
    }
    else if(spe->stgs->linsys_solver == 1){//qdldl
        // QDLDL_solve(spe->L->L->n,spe->L->L->p,spe->L->L->i,spe->L->L->x,spe->L->Dinv,b);
        // return 0;

        //for new ldl
        _ldl_solve(b, spe->L->L, spe->L->Dinv, spe->L->P, spe->L->bp);
        return 0;
    }
    else if(spe->stgs->linsys_solver == 2){//cholesky
        abip_cholsol(spe, b, n);
        return 0;
    }
    else if(spe->stgs->linsys_solver == 3){//PCG
        
        if(spe->stgs->prob_type == 3){
            // return qcp_pcg(spe, b, pcg_warm_start, n, pcg_tol);
            return qcp_pcg(spe, b, pcg_warm_start, n, pcg_tol);
        }
        else if(spe->stgs->prob_type == 4){
            return svmqp_pcg(spe, b, pcg_warm_start, n, pcg_tol);
        }
        else{
            return pcg(spe, b, pcg_warm_start, spe->stgs->rho_y, n, pcg_tol);
        }
    }
    else if(spe->stgs->linsys_solver == 4){//mkl_pardiso
        pardiso_solve(spe, b, n);
        return 0;
    }
    else if(spe->stgs->linsys_solver == 5){//lapack dense chol
        dense_chol_sol(spe, b, n);
        return 0;
    }
    else{
        printf("\nlinsys solver type error\n");
        return -1;
    }
}


abip_int ABIP(free_linsys)(spe_problem *spe){
    if(spe->L){
        if(spe->stgs->linsys_solver == 0){//mkl_dss
            MKL_INT opt = MKL_DSS_DEFAULTS;
            dss_delete(spe->L->handle, opt);
            return 0;
        }
        else if(spe->stgs->linsys_solver == 1){//qdldl
            if(spe->L->Dinv) abip_free(spe->L->Dinv);
            if(spe->L->L) cs_spfree(spe->L->L);
            if(spe->L->P) abip_free(spe->L->P);
            if(spe->L->bp) abip_free(spe->L->bp);
        return 0;
        }
        else if(spe->stgs->linsys_solver == 2){//cholesky
            if(spe->L->S) cs_sfree (spe->L->S) ;
            if(spe->L->N) cs_nfree (spe->L->N) ;
        return 0;
        }
        else if(spe->stgs->linsys_solver == 3){//pcg
            if(spe->L->M) abip_free(spe->L->M);
        return 0;
        }
        else if(spe->stgs->linsys_solver == 4){//mkl_pardiso
            pardiso_free(spe);
            return 0;
        }
        else if(spe->stgs->linsys_solver == 5){//lapack dense chol
            dense_chol_free(spe);
            return 0;
        }

        else{
        printf("\nlinsys solver type error\n");
        return -1;
        }
    }
}