#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "abip.h"
#include <math.h>

void ABIP(set_as_scaled_array)
(
	abip_float *x, 
	const abip_float *a, 
	const abip_float b, 
	abip_int len
);

void ABIP(set_as_sqrt)
(
	abip_float *x, 
	const abip_float *v, 
	abip_int len
);

void ABIP(set_as_sq)
(
	abip_float *x, 
	const abip_float *v, 
	abip_int len
);

void ABIP(scale_array)
(
	abip_float *a, 
	const abip_float b, 
	abip_int len
);


abip_float ABIP(dot)
(
	const abip_float *x, 
	const abip_float *y, 
	abip_int len
);

abip_float ABIP(norm_sq)
(
	const abip_float *v, 
	abip_int len
);

abip_float ABIP(norm)
(
	const abip_float *v, 
	abip_int len
);

abip_float ABIP(norm_inf)
(
	const abip_float *a, 
	abip_int len
);


// Added by Kurt
abip_float ABIP(norm_one)
(
      const abip_float *v, 
      abip_int len
);

// Added by Kurt
abip_float ABIP(norm_one_sqrt)
(
      const abip_float *v, 
      abip_int len
);

// Added by Kurt
abip_float ABIP(norm_inf_sqrt)
(
      const abip_float *v, 
      abip_int len
);

abip_float ABIP(min_abs_sqrt)
(
      const abip_float *a,
      abip_int len,
      abip_float ref
);
// ------


void ABIP(add_array)
(
	abip_float *a, 
	const abip_float b, 
	abip_int len
);

void ABIP(add_scaled_array)
(
	abip_float *a, 
	const abip_float *b, 
	abip_int n, 
	const abip_float sc
);

abip_float ABIP(norm_diff)
(
	const abip_float *a, 
	const abip_float *b, 
	abip_int len
);

abip_float ABIP(norm_inf_diff)
(
	const abip_float *a, 
	const abip_float *b, 
	abip_int len
);

#ifdef __cplusplus
}
#endif
#endif
