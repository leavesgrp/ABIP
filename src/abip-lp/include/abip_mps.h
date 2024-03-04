#ifndef lp_mps_h
#define lp_mps_h

#include <stdlib.h>
#include <string.h>

#ifndef POTLP_FREE
#define POTLP_FREE(var) do {if (var) {free((var)); (var) = NULL;}} while (0)
#define POTLP_INIT(var, type, size) (var) = (type *) calloc(size, sizeof(type))
#define POTLP_REALLOC(var, type, size) (var) = (type *) realloc(var, sizeof(type) * (size))
#define POTLP_MEMCPY(dst, src, type, size) memcpy(dst, src, sizeof(type) * (size))
#define POTLP_ZERO(var, type, size) memset(var, 0, sizeof(type) * (size))
#define POTLP_MAX(x, y) ((x) >= (y) ? (x) : (y))
#define POTLP_MIN(x, y) ((x) <= (y) ? (x) : (y))
#define POTLP_INFINITY          (1e+30)
#define POTLP_DIV(x, y) ((x) / (y))
#define POT_CALL(func) retcode = (func);                      \
                         if (retcode != 0) {     \
                             goto exit_cleanup;                 \
                         }

#endif

/* Implement an LP mps file reader */
int potLpMpsRead( char *fname, char *name, int *pnRow, int *pnEqRow, int *pnInEqRow, int *pnCol, int *pnElem,
                  int **peqMatBeg,  int **peqMatIdx, double **peqMatElem, int **pIneqMatBeg,
                  int **pIneqMatIdx, double **pIneqMatElem, double **prowRHS, double **pcolObj,
                  int *pnColUb, int **pcolUbIdx, double **pcolUbElem );

#endif /* lp_mps_h */
