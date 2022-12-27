#ifndef AMATRIX_H_GUARD
#define AMATRIX_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"

/* A is supplied in column compressed format */
struct ABIP_A_DATA_MATRIX 
{
       abip_float *x;    /* A values, size: NNZ A */    
       abip_int *i;      /* A row index, size: NNZ A */      
       abip_int *p;      /* A column pointer, size: n+1 */   
       abip_int m;       /* m rows, n cols */
       abip_int n; 
};

#ifdef __cplusplus
}
#endif
#endif
