#include "abip.h"
#include "cs.h"
#include "amd.h"
#include "ldl.h"
#include "direct.h"
#include "abip_pardiso.h"

static void printMat( cs *A ) {
    abip_int *Ap = A->p, *Ai = A->i; abip_float *Ax = A->x;
    abip_printf("Matrix dimension: %d by %d \n", A->m, A->n);
    abip_printf("%8s %8s %8s \n", "i", "j", "val");
    for (abip_int i = 0, j; i < A->m; ++i) {
        for (j = Ap[i]; j < Ap[i + 1]; ++j) {
            abip_printf("%8d %8d %8.2e \n", Ai[j], i, Ax[j]);
        }
    }
}

extern abip_int pardisoFactorize( ABIPLinSysWork *p, cs *A ) {
    /* Factorize the spsMat matrix */
    abip_int phase = PARDISO_SYM_FAC, error = PARDISO_OK;
    pardiso(p->pardiso_work, &maxfct, &mnum, &mtype, &phase, &A->n,
            A->x, A->p, A->i, p->P, &idummy, PARDISO_PARAMS_LDL,
            &msglvl, NULL, NULL, &error);
    
    if (!p->i) {
        p->i = (abip_int *) abip_calloc((A->m + 1), sizeof(abip_int));
    }
    
    if (!p->j) {
        p->j = (abip_int *) abip_calloc(A->p[A->m], sizeof(abip_int));
    }
    
    memcpy(p->i, A->p, sizeof(abip_int) * (A->m + 1));
    memcpy(p->j, A->i, sizeof(abip_int) * A->p[A->m]);
    
    if (error) {
        abip_printf("[Pardiso Error]: Matrix factorization failed."
               " Error code: %d \n", error);
    }
    
    return error;
}

extern void pardisoFree( ABIPLinSysWork *p, ABIPMatrix *A ) {
    /* Free the internal structure of pardiso */
    abip_int phase = PARDISO_FREE, error = PARDISO_OK, n = A->m + A->n, one = 1;
    pardiso(p->pardiso_work, &maxfct, &mnum, &mtype, &phase, &n,
            NULL, p->i, p->j, &idummy, &one,
            PARDISO_PARAMS_LDL, &msglvl, &ddummy, &ddummy, &error);
    abip_free(p->i); abip_free(p->j);
    if (error) {
        abip_printf("[Pardiso Error]: Pardiso free failed."
               " Error code %d \n", error);
    }
}

extern void pardisoSolve( ABIPLinSysWork *p, ABIPMatrix *A, abip_float *b ) {
    /* Solve the linear system S * X = B using Pardiso */
    // Invoke pardiso to perform solution
    abip_int phase = PARDISO_SOLVE, error = PARDISO_OK, nrhs = 1, n = A->m + A->n;
    pardiso(p->pardiso_work, &maxfct, &mnum, &mtype, &phase, &n,
            NULL, p->i, p->j, &idummy, &nrhs, PARDISO_PARAMS_LDL,
            &msglvl, b, p->D, &error);
    if (error) {
        abip_printf("[Pardiso Error]: Pardiso solve failed."
               " Error code %d \n", error);
    }
}

