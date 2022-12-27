#ifndef LINSYS_H_GUARD
#define LINSYS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "abip.h"

ABIPLinSysWork *ABIP(init_lin_sys_work)(const ABIPMatrix *A, const ABIPSettings *stgs);

abip_int ABIP(solve_lin_sys)
(
	const ABIPMatrix *A, 
	const ABIPSettings *stgs,
    ABIPLinSysWork *p, 
    abip_float *b, 
    const abip_float *s,
    abip_int iter
);

void ABIP(free_lin_sys_work)
(
	ABIPLinSysWork *p
);

void ABIP(free_lin_sys_work_pds)
(
    ABIPLinSysWork *p,
    ABIPMatrix *A
);

/* forms y += A'*x */
void ABIP(accum_by_Atrans)
(
	const ABIPMatrix *A, 
	ABIPLinSysWork *p, 
	const abip_float *x,
    abip_float *y
);

/* forms y += A*x */
void ABIP(accum_by_A)
(
	const ABIPMatrix *A, 
	ABIPLinSysWork *p, 
	const abip_float *x,
    abip_float *y
);

abip_int ABIP(validate_lin_sys)
(
	const ABIPMatrix *A
);

char *ABIP(get_lin_sys_method)
(
	const ABIPMatrix *A, 
	const ABIPSettings *stgs
);

char *ABIP(get_lin_sys_summary)
(
	ABIPLinSysWork *p, 
	const ABIPInfo *info
);

void ABIP(normalize_A)
(
	ABIPMatrix *A, 
	const ABIPSettings *stgs, 
    ABIPScaling *scal
);

void ABIP(un_normalize_A)
(
	ABIPMatrix *A, 
	const ABIPSettings *stgs,
	const ABIPScaling *scal
);

void ABIP(free_A_matrix)
(
	ABIPMatrix *A
);

abip_int ABIP(copy_A_matrix)
(
	ABIPMatrix **dstp, 
	const ABIPMatrix *src
);

#ifdef __cplusplus
}
#endif

#endif
