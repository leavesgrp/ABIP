#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include "common.h"
#include "glbopts.h"
#include "linalg.h"
#include "abip.h"

struct ABIP_LIN_SYS_WORK 
{
    abip_float *p;    /* cg iterate  */
    abip_float *r;    /* cg residual */
    abip_float *Gp;
    abip_float *tmp;
    ABIPMatrix *At;
    
    /* preconditioning */
    abip_float *z;
    abip_float *M;
    
    /* reporting */
    abip_int tot_cg_its;
    abip_float total_solve_time;
};

#ifdef __cplusplus
}
#endif
#endif
