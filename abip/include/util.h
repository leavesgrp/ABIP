#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "abip.h"
#include <stdlib.h>
#include <stdio.h>

#if (defined NOTIMER)
typedef void *ABIP(timer);

#elif(defined _WIN32 || defined _WIN64 || defined _WINDLL)

#include <windows.h>
typedef struct ABIP(timer) 
{
      LARGE_INTEGER tic;
      LARGE_INTEGER toc;
      LARGE_INTEGER freq;
} ABIP(timer);

#elif(defined __APPLE__)

#include <mach/mach_time.h>
typedef struct ABIP(timer) 
{
      uint64_t tic;
      uint64_t toc;
      mach_timebase_info_data_t tinfo;
} ABIP(timer);

#else

#include <time.h>
typedef struct ABIP(timer) 
{
      struct timespec tic;
      struct timespec toc;
} ABIP(timer);

#endif

#if EXTRA_VERBOSE > 1
extern ABIP(timer) global_timer;
#endif

void ABIP(tic)(ABIP(timer) *t);
abip_float ABIP(toc)(ABIP(timer) *t);
abip_float ABIP(str_toc)(char *str, ABIP(timer) *t);
abip_float ABIP(tocq)(ABIP(timer) *t);

void ABIP(print_data)(const ABIPData *d);
void ABIP(print_work)(const ABIPWork *w);
void ABIP(print_array)(const abip_float *arr, abip_int n, const char *name);
void ABIP(set_default_settings)(ABIPData *d);
void ABIP(free_sol)(ABIPSolution *sol);
void ABIP(free_data)(ABIPData *d);

#ifdef __cplusplus
}
#endif
#endif
