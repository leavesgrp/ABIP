#include "potlp_cs.h"
#include "abip_mps.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

/* CSparse routines for reading inputs. Referenced from Tim Davis Suite Sparse */
void *pot_dcs_malloc (int n, size_t size) {
    
    return (malloc (POTLP_MAX (n, 1) * size)) ;
}

/* wrapper for calloc */
void *pot_dcs_calloc (int n, size_t size) {
    
    return (calloc (POTLP_MAX (n,1), size)) ;
}

/* wrapper for free */
void *pot_dcs_free (void *p) {
    
    if (p) free (p) ;       /* free p if it is not already NULL */
    return (NULL) ;         /* return NULL to simplify the use of dpot_dcs_free */
}

/* wrapper for realloc */
void *pot_dcs_realloc (void *p, int n, size_t size, int *ok) {
    
    void *pnew ;
    pnew = realloc (p, POTLP_MAX (n,1) * size) ; /* realloc the block */
    *ok = (pnew != NULL) ;                  /* realloc fails if pnew is NULL */
    return ((*ok) ? pnew : p) ;             /* return original p if failure */
}

pot_dcs *pot_dcs_spalloc (int m, int n, int nzmax, int values, int triplet) {
    
    pot_dcs *A = pot_dcs_calloc (1, sizeof (pot_dcs)) ;    /* allocate the pot_dcs struct */
    if (!A) return (NULL) ;                 /* out of memory */
    A->m = m ;                              /* define dimensions and nzmax */
    A->n = n ;
    A->nzmax = nzmax = POTLP_MAX (nzmax, 1) ;
    A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col */
    A->p = pot_dcs_malloc (triplet ? nzmax : n+1, sizeof (int)) ;
    A->i = pot_dcs_malloc (nzmax, sizeof (int)) ;
    A->x = values ? pot_dcs_malloc (nzmax, sizeof (double)) : NULL ;
    return ((!A->p || !A->i || (values && !A->x)) ? pot_dcs_spfree (A) : A) ;
}

/* change the max # of entries sparse matrix */
int pot_dcs_sprealloc (pot_dcs *A, int nzmax) {
    
    int ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = IS_CSC(A) ? (A->p [A->n]) : A->nz ;
    nzmax = POTLP_MAX (nzmax, 1) ;
    A->i = pot_dcs_realloc (A->i, nzmax, sizeof (int), &oki) ;
    if (IS_TRIPLET (A)) A->p = pot_dcs_realloc (A->p, nzmax, sizeof (int), &okj) ;
    if (A->x) A->x = pot_dcs_realloc (A->x, nzmax, sizeof (double), &okx) ;
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/* free a sparse matrix */
pot_dcs *pot_dcs_spfree (pot_dcs *A) {
    
    if (!A) return (NULL) ;     /* do nothing if A already NULL */
    pot_dcs_free (A->p) ;
    pot_dcs_free (A->i) ;
    pot_dcs_free (A->x) ;
    return ((pot_dcs *) pot_dcs_free (A)) ;   /* free the pot_dcs struct and return NULL */
}

/* free workspace and return a sparse matrix result */
pot_dcs *pot_dcs_done (pot_dcs *C, void *w, void *x, int ok) {
    
    pot_dcs_free (w) ;                       /* free workspace */
    pot_dcs_free (x) ;
    return (ok ? C : pot_dcs_spfree (C)) ;   /* return result if OK, else free it */
}

double pot_dcs_cumsum (int *p, int *c, int n) {
    
    int i, nz = 0 ;
    double nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++) {
        p [i] = nz ;
        nz += c [i] ;
        nz2 += c [i] ;              /* also in double to avoid int overflow */
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz2) ;                  /* return sum (c [0..n-1]) */
}

/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
int pot_dcs_entry (pot_dcs *T, int i, int j, double x) {
    
    if (!IS_TRIPLET (T) || i < 0 || j < 0) return (0) ;     /* check inputs */
    if (T->nz >= T->nzmax && !pot_dcs_sprealloc (T,2*(T->nzmax))) return (0);
    if (T->x) T->x [T->nz] = x ;
    T->i [T->nz] = i ;
    T->p [T->nz++] = j ;
    T->m = POTLP_MAX (T->m, i+1) ;
    T->n = POTLP_MAX (T->n, j+1) ;
    return (1) ;
}

pot_dcs *pot_dcs_compress (const pot_dcs *T) {
    
    int m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
    double *Cx, *Tx ;
    pot_dcs *C ;
    if (!IS_TRIPLET (T)) return (NULL) ;                /* check inputs */
    m = T->m ; n = T->n ; Ti = T->i ; Tj = T->p ; Tx = T->x ; nz = T->nz ;
    C = pot_dcs_spalloc (m, n, nz, Tx != NULL, 0) ;          /* allocate result */
    w = pot_dcs_calloc (n, sizeof (int)) ;                   /* get workspace */
    if (!C || !w) return (pot_dcs_done (C, w, NULL, 0)) ;    /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (k = 0 ; k < nz ; k++) w [Tj [k]]++ ;           /* column counts */
    pot_dcs_cumsum (Cp, w, n) ;                              /* column pointers */
    for (k = 0 ; k < nz ; k++) {
        Ci [p = w [Tj [k]]++] = Ti [k] ;    /* A(i,j) is the pth entry in C */
        if (Cx) Cx [p] = Tx [k] ;
    }
    return (pot_dcs_done (C, w, NULL, 1)) ;      /* success; free w and return C */
}

double pot_dcs_norm (const pot_dcs *A) {
    
    int p, j, n, *Ap ;
    double *Ax ;
    double nrm = 0, s ;
    if (!IS_CSC (A) || !A->x) return (-1) ;             /* check inputs */
    n = A->n ; Ap = A->p ; Ax = A->x ;
    for (j = 0 ; j < n ; j++) {
        for (s = 0, p = Ap [j] ; p < Ap [j+1] ; p++) s += fabs (Ax [p]) ;
        nrm = POTLP_MAX (nrm, s) ;
    }
    return (nrm) ;
}

int pot_dcs_print (const pot_dcs *A, int brief) {
    
    int p, j, m, n, nzmax, nz, *Ap, *Ai ;
    double *Ax ;
    if (!A) { printf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    if (nz < 0) {
        printf ("%g-by-%g, nzmax: %g nnz: %g, 1-norm: %g\n", (double) m,
            (double) n, (double) nzmax, (double) (Ap [n]), pot_dcs_norm (A)) ;
        for (j = 0 ; j < n ; j++) {
            printf ("    col %g : locations %g to %g\n", (double) j,
                (double) (Ap [j]), (double) (Ap [j+1]-1)) ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++) {
                printf ("      %g : ", (double) (Ai [p])) ;
                printf ("%50.50e \n", Ax ? Ax [p] : 1) ;
                if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
            }
        }
    }
    else {
        printf ("triplet: %g-by-%g, nzmax: %g nnz: %g\n", (double) m,
            (double) n, (double) nzmax, (double) nz) ;
        for (p = 0 ; p < nz ; p++) {
            printf ("    %g %g : ", (double) (Ai [p]), (double) (Ap [p])) ;
            printf ("%g\n", Ax ? Ax [p] : 1) ;
            if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
        }
    }
    return (1) ;
}
