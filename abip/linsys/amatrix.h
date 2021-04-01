#ifndef AMATRIX_H_GUARD
#define AMATRIX_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"

struct ABIP_A_DATA_MATRIX 
{
       abip_float *x;        
       abip_int *i;            
       abip_int *p;         
       abip_int m; 
       abip_int n; 
};

#ifdef __cplusplus
}
#endif
#endif
