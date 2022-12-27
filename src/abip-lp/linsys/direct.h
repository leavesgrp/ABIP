#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "abip.h"
#include "cs.h"
#include "amd.h"
#include "ldl.h"
#include "common.h"
#include "abip_pardiso.h"

struct ABIP_LIN_SYS_WORK
{
  	cs *L;         		
  	abip_float *D;  		
  	abip_int *P;
    abip_int *i;
    abip_int *j;
  	abip_float *bp;
    
    void *pardiso_work[PARDISOINDEX];
  	abip_float total_solve_time;
};

#ifdef __cplusplus
}
#endif
#endif
