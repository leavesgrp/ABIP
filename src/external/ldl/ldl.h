/* ========================================================================== */
/* === ldl.h:  include file for the LDL package ============================= */
/* ========================================================================== */

/* Copyright (c) Timothy A Davis, http://www.suitesparse.com.
 * All Rights Reserved.  See LDL/Doc/License.txt for the License.
 */

#include "SuiteSparse_config.h"

#ifdef DLONG
#define LDL_int SuiteSparse_long
#define LDL_ID SuiteSparse_long_id

#define LDL_symbolic ldl_l_symbolic
#define LDL_numeric ldl_l_numeric
#define LDL_lsolve ldl_l_lsolve
#define LDL_dsolve ldl_l_dsolve
#define LDL_ltsolve ldl_l_ltsolve
#define LDL_perm ldl_l_perm
#define LDL_permt ldl_l_permt
#define LDL_valid_perm ldl_l_valid_perm
#define LDL_valid_matrix ldl_l_valid_matrix

#else
#define LDL_int int
#define LDL_ID "%d"

#define LDL_symbolic ldl_symbolic
#define LDL_numeric ldl_numeric
#define LDL_lsolve ldl_lsolve
#define LDL_dsolve ldl_dsolve
#define LDL_ltsolve ldl_ltsolve
#define LDL_perm ldl_perm
#define LDL_permt ldl_permt
#define LDL_valid_perm ldl_valid_perm
#define LDL_valid_matrix ldl_valid_matrix

#endif

/* ========================================================================== */
/* === int version ========================================================== */
/* ========================================================================== */

void ldl_symbolic (int n, int Ap [ ], int Ai [ ], int Lp [ ],
    int Parent [ ], int Lnz [ ], int Flag [ ], int P [ ], int Pinv [ ]) ;

int ldl_numeric (int n, int Ap [ ], int Ai [ ], abip_float Ax [ ],
    int Lp [ ], int Parent [ ], int Lnz [ ], int Li [ ], abip_float Lx [ ],
    abip_float D [ ], abip_float Y [ ], int Pattern [ ], int Flag [ ],
    int P [ ], int Pinv [ ]) ;

void ldl_lsolve (int n, abip_float X [ ], int Lp [ ], int Li [ ],
    abip_float Lx [ ]) ;

void ldl_dsolve (int n, abip_float X [ ], abip_float D [ ]) ;

void ldl_ltsolve (int n, abip_float X [ ], int Lp [ ], int Li [ ],
    abip_float Lx [ ]) ;

void ldl_perm  (int n, abip_float X [ ], abip_float B [ ], int P [ ]) ;
void ldl_permt (int n, abip_float X [ ], abip_float B [ ], int P [ ]) ;

int ldl_valid_perm (int n, int P [ ], int Flag [ ]) ;
int ldl_valid_matrix ( int n, int Ap [ ], int Ai [ ]) ;

/* ========================================================================== */
/* === long version ========================================================= */
/* ========================================================================== */

void ldl_l_symbolic (SuiteSparse_long n, SuiteSparse_long Ap [ ],
    SuiteSparse_long Ai [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Parent [ ], SuiteSparse_long Lnz [ ],
    SuiteSparse_long Flag [ ], SuiteSparse_long P [ ],
    SuiteSparse_long Pinv [ ]) ;

SuiteSparse_long ldl_l_numeric (SuiteSparse_long n, SuiteSparse_long Ap [ ],
    SuiteSparse_long Ai [ ], abip_float Ax [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Parent [ ], SuiteSparse_long Lnz [ ],
    SuiteSparse_long Li [ ], abip_float Lx [ ], abip_float D [ ], abip_float Y [ ],
    SuiteSparse_long Pattern [ ], SuiteSparse_long Flag [ ],
    SuiteSparse_long P [ ], SuiteSparse_long Pinv [ ]) ;

void ldl_l_lsolve (SuiteSparse_long n, abip_float X [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Li [ ], abip_float Lx [ ]) ;

void ldl_l_dsolve (SuiteSparse_long n, abip_float X [ ], abip_float D [ ]) ;

void ldl_l_ltsolve (SuiteSparse_long n, abip_float X [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Li [ ], abip_float Lx [ ]) ;

void ldl_l_perm  (SuiteSparse_long n, abip_float X [ ], abip_float B [ ],
    SuiteSparse_long P [ ]) ;
void ldl_l_permt (SuiteSparse_long n, abip_float X [ ], abip_float B [ ],
    SuiteSparse_long P [ ]) ;

SuiteSparse_long ldl_l_valid_perm (SuiteSparse_long n, SuiteSparse_long P [ ],
    SuiteSparse_long Flag [ ]) ;
SuiteSparse_long ldl_l_valid_matrix ( SuiteSparse_long n,
    SuiteSparse_long Ap [ ], SuiteSparse_long Ai [ ]) ;

/* ========================================================================== */
/* === LDL version ========================================================== */
/* ========================================================================== */

#define LDL_DATE "May 4, 2016"
#define LDL_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define LDL_MAIN_VERSION 2
#define LDL_SUB_VERSION 2
#define LDL_SUBSUB_VERSION 6
#define LDL_VERSION LDL_VERSION_CODE(LDL_MAIN_VERSION,LDL_SUB_VERSION)

