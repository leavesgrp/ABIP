#ifndef CONES_H_GUARD
#define CONES_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "abip.h"
#include "linalg.h"
#include "mkl.h"
#include "mkl_lapacke.h"



char* ABIP(get_cone_header)(const ABIPCone* k);
abip_int ABIP(validate_cones)(spe_problem *spe, const ABIPCone* k);

abip_int ABIP(get_cone_dims)(const ABIPCone* k);


void ABIP(soc_barrier_subproblem)(
    abip_float* x,
    abip_float* tmp,
    abip_float lambda,
    abip_int n
);

void ABIP(rsoc_barrier_subproblem)(
    abip_float* x,
    abip_float* tmp,
    abip_float lambda,
    abip_int n 
);

void ABIP(free_barrier_subproblem)(
    abip_float* x,
    abip_float* tmp,
    abip_float lambda,
    abip_int n 
);

void ABIP(zero_barrier_subproblem)(
    abip_float* x,
    abip_float* tmp,
    abip_float lambda,
    abip_int n 
);

void ABIP(positive_orthant_barrier_subproblem)(
    abip_float* x,
    abip_float* tmp,
    abip_float lambda,
    abip_int n 
);


#ifdef __cplusplus
}
#endif
#endif

