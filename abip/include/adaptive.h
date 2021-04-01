#ifndef ADAPT_H_GUARD
#define ADAPT_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "abip.h"
#include "glbopts.h"
#include <math.h>

ABIPAdaptWork *ABIP(init_adapt)
(
	ABIPWork *w
);

void ABIP(free_adapt)
(
	ABIPAdaptWork *a
);

abip_int ABIP(adaptive)
(
	ABIPWork *w, 
	abip_int iter
);

char *ABIP(get_adapt_summary)
(
	const ABIPInfo *info, 
	ABIPAdaptWork *a
);

#ifdef __cplusplus
}
#endif
#endif
