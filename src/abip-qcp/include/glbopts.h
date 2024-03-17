#ifndef GLB_H_GUARD
#define GLB_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

//#define DLONG

#ifndef ABIP
#define ABIP(x) abip_##x
#endif

/* ABIP VERSION NUMBER ----------------------------------------------    */
#define ABIP_VERSION                                                            \
  ("2.0.0") /* string literals automatically null-terminated */

#define ABIP_INFEASIBLE_INACCURATE (-7)  
#define ABIP_UNBOUNDED_INACCURATE (-6)  
#define ABIP_SIGINT (-5) 
#define ABIP_FAILED (-4)  
#define ABIP_INDETERMINATE (-3)  
#define ABIP_INFEASIBLE (-2) 
#define ABIP_UNBOUNDED (-1) 
#define ABIP_UNFINISHED (0) 
#define ABIP_SOLVED (1) 
#define ABIP_SOLVED_INACCURATE (2)
#define ABIP_UNSOLVED (3) 
#define SIGMA (0.8)
#define GAMMA (1.6)

#define MAX_IPM_ITERS (500)
#define MAX_ADMM_ITERS (10000000)
#define EPS (1E-3)
#define ALPHA (1.8)
#define CG_RATE (2.0)
#define CG_BEST_TOL (1e-9)
#define CG_MIN_TOL (1e-5)
#define NORMALIZE (1)
#define SCALE (1.0)
#define SPARSITY_RATIO (0.01)
#define RHO_Y (1E-3)
#define ADAPTIVE (1)
#define EPS_COR (0.2)
#define EPS_PEN (0.1)
#define ADAPTIVE_LOOKBACK (20)
#define VERBOSE (1)
#define WARM_START (0)

#define CONE_TOL (1e-8)

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define abip_printf mexPrintf
#define _abip_free mxFree
#define _abip_malloc mxMalloc
#define _abip_calloc mxCalloc
#define _abip_realloc mxRealloc
#elif defined PYTHON
#include <Python.h>
#include <stdlib.h>
#define abip_printf(...)                                                            \
{                                                                                               \
      PyGILState_STATE gilstate = PyGILState_Ensure();       \
      PySys_WriteStdout(__VA_ARGS__);                                 \
      PyGILState_Release(gilstate);                                          \
}
#define _abip_free free
#define _abip_malloc malloc
#define _abip_calloc calloc
#define _abip_realloc realloc
#else
#include <stdio.h>
#include <stdlib.h>
#define abip_printf printf
#define _abip_free free
#define _abip_malloc malloc
#define _abip_calloc calloc
#define _abip_realloc realloc
#endif

#define abip_free(x)   \
      _abip_free(x);          \
      x = ABIP_NULL 
#define abip_malloc(x) _abip_malloc(x)
#define abip_calloc(x, y) _abip_calloc(x, y)
#define abip_realloc(x, y) _abip_realloc(x, y)

// //#ifdef DLONG
// //#ifdef _WIN64
// //typedef __int64 abip_int; 
// //#else
// //typedef long abip_int;
// //#endif
// //#else
// //typedef int abip_int;
// //#endif
// typedef int abip_int;


#ifdef DLONG
/*#ifdef _WIN64
#include <stdint.h>
typedef int64_t abip_int;
#else
typedef long abip_int;
#endif
*/
typedef long long abip_int;
#else
typedef int abip_int;
#endif


#ifndef SFLOAT
typedef double abip_float;
#ifndef NAN
#define NAN ((abip_float)0x7ff8000000000000)
#endif
#ifndef INFINITY
#define INFINITY NAN
#endif
#else
typedef float abip_float;
#ifndef NAN
#define NAN ((float)0x7fc00000)
#endif
#ifndef INFINITY
#define INFINITY NAN
#endif
#endif

#define ABIP_NULL 0

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(x) (((x) < 0) ? -(x) : (x))
#endif

#ifndef POWF
#ifdef SFLOAT
#define POWF powf
#else
#define POWF pow
#endif
#endif

#ifndef SQRTF
#ifdef SFLOAT
#define SQRTF sqrtf
#else
#define SQRTF sqrt
#endif
#endif

#if EXTRA_VERBOSE > 1
#if (defined _WIN32 || defined _WIN64 || defined _WINDLL)
#define __func__ __FUNCTION__
#endif
#define DEBUG_FUNC abip_printf("IN function: %s, time: %4f ms, file: %s, line: %i\n", __func__,  ABIP(tocq)(&global_timer), __FILE__, __LINE__);
#define RETURN 
      abip_printf("EXIT function: %s, time: %4f ms, file: %s, line: %i\n", __func__, ABIP(tocq)(&global_timer), __FILE__, __LINE__);   \
      return
#else
#define DEBUG_FUNC
#define RETURN return
#endif

#define EPS_TOL (1E-18)
#define SAFEDIV_POS(X, Y) ((Y) < EPS_TOL ? ((X) / EPS_TOL) : (X) / (Y))

#define CONVERGED_INTERVAL (1)
#define INDETERMINATE_TOL (1e-9)

#ifdef __cplusplus
}
#endif
#endif
