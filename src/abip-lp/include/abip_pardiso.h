#ifndef abip_pardiso_h
#define abip_pardiso_h

/* Implement the pardiso solver interface for IPM */

#include <stdio.h>
#include "glbopts.h"

#define PARDISO_OK       (  0)
#define PARDISOINDEX     ( 64)      // Pardiso working array length
#define PARDISO_SYM      ( 11)      // Pardiso symbolic analysis
#define PARDISO_FAC      ( 22)      // Pardiso numerical factorization
#define PARDISO_SYM_FAC  ( 12)      // Symbolic analysis and factorization
#define PARDISO_SOLVE    ( 33)      // Solve linear system
#define PARDISO_FORWARD  (331)      // Pardiso forward solve
#define PARDISO_BACKWARD (333)      // Pardiso backward solve
#define PARDISO_FREE     ( -1)      // Free internal data structure

// Pardiso default parameters
static abip_int maxfct = 1;  // Maximum number of factors
static abip_int mnum   = 1;  // The matrix used for the solution phase
static abip_int mtype  = -2; // Real and symmetric indefinite
static abip_int msglvl = 0;  // Print no information
static abip_int idummy = 0;  // Dummy variable for taking up space
static abip_float ddummy = 0.0;  // Dummy variable for taking up space

// Pardiso solver
#define SYMBOLIC 3
#define PIVOTING 2
#define FACTORIZE 0

static abip_int PARDISO_PARAMS_LDL[PARDISOINDEX] = {
    
    1, /* Non-default value */ SYMBOLIC, /* P Nested dissection */ 0, /* Reserved          */
    0, /* No CG             */ 0, /* No user permutation */ 1, /* Overwriting    */
    0, /* Refinement report */ 0, /* Auto ItRef step     */ 0, /* Reserved          */
    8, /* Perturb           */ 0, /* Disable scaling     */ 0, /* No transpose      */
    0, /* Disable matching  */ 0, /* Report on pivots    */ 0, /* Output            */
    0, /* Output            */ 0, /* Output              */-1, /* No report         */
    0, /* No report         */ 0, /* Output              */ PIVOTING, /* Pivoting          */
    0, /* nPosEigVals       */ 0, /* nNegEigVals         */ FACTORIZE, /* Classic factorize */
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

#ifdef __cplusplus
extern "C" {
#endif

extern void pardisoinit ( void *, const abip_int *, abip_int * );

extern void pardiso     ( void     *, abip_int    *, abip_int *, abip_int *, abip_int *, abip_int *,
                          double   *, abip_int    *, abip_int *, abip_int *, abip_int *, abip_int *,
                          abip_int *, double      *, double   *, abip_int * );

extern abip_int pardisoFactorize( ABIPLinSysWork *p, cs *A );
extern void pardisoFree         ( ABIPLinSysWork *p, ABIPMatrix *A );
extern void pardisoSolve        ( ABIPLinSysWork *p, ABIPMatrix *A, abip_float *b );
#ifdef __cplusplus
}
#endif


#endif /* abip_pardiso_h */
