#ifndef COMMON_H_GUARD
#define COMMON_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "abip.h"
#include "amatrix.h"
#include "linsys.h"
#include "linalg.h"
#include "util.h"

void ABIP(_accum_by_Atrans)
(
	abip_int n, 
	abip_float *Ax, 
	abip_int *Ai, 
	abip_int *Ap, 
	const abip_float *x, 
	abip_float *y
);
 
void ABIP(_accum_by_A)
(
	abip_int n, 
	abip_float *Ax, 
	abip_int *Ai, 
	abip_int *Ap, 
	const abip_float *x, 
	abip_float *y
);

void ABIP(_normalize_A)
(
    ABIPMatrix *A, 
    const ABIPSettings *stgs, 
    ABIPScaling *scal
); 

void ABIP(_un_normalize_A)
(
    ABIPMatrix *A, 
    const ABIPSettings *stgs, 
    const ABIPScaling *scal
);

abip_float ABIP(cumsum)
(
	abip_int *p, 
	abip_int *c, 
	abip_int n
); 

#ifdef __cplusplus
}
#endif

#endif
