#ifndef NORMALIZE_H_GUARD
#define NORMALIZE_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "abip.h"

void ABIP(normalize_b_c)
(
	ABIPWork *w
);

void ABIP(calc_scaled_resids)
(
	ABIPWork *w, 
	ABIPResiduals *r
);

void ABIP(normalize_warm_start)
(
	ABIPWork *w
);

void ABIP(un_normalize_sol)
(
	ABIPWork *w, 
	ABIPSolution *sol
);

#ifdef __cplusplus
}
#endif
#endif
