#ifndef CS_H_GUARD
#define CS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"

typedef struct ABIP(cs_sparse) 
{
       abip_int nnzmax;  		
       abip_int m;  		
       abip_int n;  		
       abip_int *p;  		
       abip_int *i;  		
       abip_float *x;  		
       abip_int nnz;  		
} cs;

cs *ABIP(cs_compress)
(
	const cs *T
); 

cs *ABIP(cs_spalloc)
(
	abip_int m, 
	abip_int n, 
	abip_int nnzmax, 
	abip_int values, 
	abip_int triplet
);

cs *ABIP(cs_spfree)
(
	cs *A
);

abip_float ABIP(cs_cumsum)
(
	abip_int *p, 
	abip_int *c, 
	abip_int n
);

abip_int *ABIP(cs_pinv)
(
	abip_int const *p, 
	abip_int n
);

cs *ABIP(cs_symperm)
(
	const cs *A, 
	const abip_int *pinv, 
	abip_int values
);

#ifdef __cplusplus
}
#endif
#endif
