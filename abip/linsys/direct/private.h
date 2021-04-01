#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "abip.h"
#include "cs.h"
#include "external/amd.h"
#include "external/ldl.h"
#include "common.h"

struct ABIP_LIN_SYS_WORK 
{
  	cs *L;         		
  	abip_float *D;  		
  	abip_int *P;    		
  	abip_float *bp; 		
  	
  	abip_float total_solve_time;
};

#ifdef __cplusplus
}
#endif
#endif
