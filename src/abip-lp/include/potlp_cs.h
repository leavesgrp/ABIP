/* ========================================================================== */
/* CXSparse/Include/cs.h file */
/* ========================================================================== */

/* This is the CXSparse/Include/cs.h file.  It has the same name (cs.h) as
   the CSparse/Include/cs.h file.  The 'make install' for SuiteSparse installs
   CXSparse, and this file, instead of CSparse.  The two packages have the same
   cs.h include filename, because CXSparse is a superset of CSparse.  Any user
   program that uses CSparse can rely on CXSparse instead, with no change to the
   user code.  The #include "cs.h" line will work for both versions, in user
   code, and the function names and user-visible typedefs from CSparse all
   appear in CXSparse.  For experimenting and changing the package itself, I
   recommend using CSparse since it's simpler and easier to modify.  For
   using the package in production codes, I recommend CXSparse since it has
   more features (support for complex matrices, and both int and long
   versions).
 */

/* ========================================================================== */

#ifndef _HDSDP_CS_
#define _HDSDP_CS_

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pot_cs_sparse {
    int nzmax;
    int m;         /* number of rows */
    int n;         /* number of columns */
    int *p;        /* column pointers (size n+1) or col indices (size nzmax) */
    int *i;        /* row indices, size nzmax */
    double *x;     /* numerical values, size nzmax */
    int nz;        /* # of entries in triplet matrix, -1 for compressed-col */
} pot_dcs ;

int pot_dcs_entry (pot_dcs *T, int i, int j, double x) ;
pot_dcs *pot_dcs_compress (const pot_dcs *T) ;
double pot_dcs_norm (const pot_dcs *A) ;
int pot_dcs_print (const pot_dcs *A, int brief) ;

/* utilities */
void *pot_dcs_calloc (int n, size_t size) ;
void *pot_dcs_free (void *p) ;
void *pot_dcs_realloc (void *p, int n, size_t size, int *ok) ;
pot_dcs *pot_dcs_spalloc (int m, int n, int nzmax, int values, int t) ;
pot_dcs *pot_dcs_spfree (pot_dcs *A) ;
int pot_dcs_sprealloc (pot_dcs *A, int nzmax) ;
void *pot_dcs_malloc (int n, size_t size) ;

/* utilities */
double pot_dcs_cumsum (int *p, int *c, int n) ;
pot_dcs *pot_dcs_done (pot_dcs *C, void *w, void *x, int ok) ;
int *pot_dcs_idone (int *p, pot_dcs *C, void *w, int ok) ;

#define IS_CSC(A) (A && (A->nz == -1))
#define IS_TRIPLET(A) (A && (A->nz >= 0))

#ifdef __cplusplus
}
#endif

#endif
