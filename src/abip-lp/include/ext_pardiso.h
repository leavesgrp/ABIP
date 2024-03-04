#ifndef pardiso_h
#define pardiso_h

#define PARDISO_RET_OK         ( 0)
#define PARDISO_SYM_INDEFINITE (-2)
#define PARDISO_PHASE_SYM      (11)      // Pardiso symbolic analysis
#define PARDISO_PHASE_SYM_FAC  (12)      // Pardiso symbolic analysis
#define PARDISO_PHASE_FAC      (22)      // Pardiso numerical factorization
#define PARDISO_PHASE_SOLVE    (33)      // Solve linear system
#define PARDISO_PHASE_FREE     (-1)      // Free internal data structure

#define PARDISO_PARAM_NONDEFAULT    (0)
#define PARDISO_PARAM_SYMBOLIC      (1)
#define PARDISO_PARAM_SYMBOLIC_MMD  (0)
#define PARDISO_PARAM_SYMBOLIC_ND   (2)
#define PARDISO_PARAM_REFINEMENT    (7)
#define PARDISO_PARAM_INPLACE       (5)
#define PARDISO_PARAM_PERTURBATION  (9)
#define PARDISO_PARAM_SCALING      (10)
#define PARDISO_PARAM_MATCHING     (12)
#define PARDISO_PARAM_FACNNZ       (17)
#define PARDISO_PARAM_FACFLOP      (18)
#define PARDISO_PARAM_CNR          (33)
#define PARDISO_PARAM_INDEX        (34)
#define PARDISO_PARAM_INDEX_C       (1)

#define set_pardiso_param(iparm, param, val) iparm[param] = val
#define get_pardiso_output(iparm, param) iparm[param]

#ifdef POTLP64
#ifdef HEADERPATH
#include "include/pot_def.h"
#else
#include "pot_def.h"
#endif
extern void pardisoinit ( void *, int *, int * );
extern void pardiso     ( void     *, int    *, int *, int *, int *, int *,
                          double   *, int    *, int *, int *, int *, int *,
                          int *, double      *, double   *, int * );
#else
extern void pardisoinit ( void *, int *, int * );
extern void pardiso     ( void     *, int    *, int *, int *, int *, int *,
                          double   *, int    *, int *, int *, int *, int *,
                          int *, double      *, double   *, int * );
#endif

#endif /* pardiso_h */
