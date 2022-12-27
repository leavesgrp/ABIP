#ifndef ABIP_BLAS_H_GUARD
#define ABIP_BLAS_H_GUARD

#ifdef USE_LAPACK

#ifdef __cplusplus
extern "C" {
#endif

#ifndef BLASSUFFIX
#define BLASSUFFIX _
#endif

#if defined(NOBLASSUFFIX) && NOBLASSUFFIX > 0
#ifndef SFLOAT
#define BLAS(x) d##x
#else
#define BLAS(x) s##x
#endif
#else
#define stitch_(pre, x, post) pre##x##post
#define stitch__(pre, x, post) stitch_(pre, x, post)
#ifndef SFLOAT
#define BLAS(x) stitch__(d, x, BLASSUFFIX)
#else
#define BLAS(x) stitch__(s, x, BLASSUFFIX)
#endif
#endif

#ifdef MATLAB_MEX_FILE
typedef ptrdiff_t blas_int;
#elif defined BLAS64
#include <stdint.h>
typedef int64_t blas_int;
#else
typedef int blas_int;
#endif

#ifdef __cplusplus
}
#endif

#endif

#endif
