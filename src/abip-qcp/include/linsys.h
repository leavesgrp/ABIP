#ifndef LINSYS_H_GUARD
#define LINSYS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include "abip.h"
#include "amatrix.h"
#include "linalg.h"
#include "glbopts.h"
#include "cs.h"
#include "mkl.h"
#include "mkl_dss.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_lapacke.h"
#include "util.h"
#include "qdldl.h"


struct ABIP_LIN_SYS_WORK
{
	abip_float total_solve_time;

	//matrix to factorize
	cs *K;

	//mkl-dss
	_MKL_DSS_HANDLE_t handle;


	//mkl-pardiso
	void *pt[64];
	MKL_INT iparm[64];
    MKL_INT maxfct, mnum, error, msglvl;
	abip_float ddum;          /* Double dummy */
    MKL_INT idum; 
	MKL_INT mtype; 



	//pcg
	abip_float* M;     //preconditioner for pcg
	abip_int total_cg_iters;
	

	//sparse cholesky
	css *S ;
    csn *N ;


	//qdldl
	cs *L;         	  //matrix L of LDL' factorization	
  	abip_float *Dinv;    //the diagonal vector of the diagonal matrix D L of LDL' factorization		
	abip_int nnz_LDL;
	abip_int *P;  //permutation for KKT matrix
  	abip_float *bp;	


	//lapack dense cholesky
	// A = UTU
	abip_float *U;
};

/*   y += A'*x
		A in column compressed format
		parallelizes over columns (rows of A')
*/
void ABIP(accum_by_Atrans)
(
	const ABIPMatrix* A,
	const abip_float* x,
	abip_float* y
);

/*	y += A*x
	A in column compressed format
	this parallelizes over columns and uses
	pragma atomic to prevent concurrent writes to y
*/
void ABIP(accum_by_A)
(
	const ABIPMatrix* A,
	const abip_float* x,
	abip_float* y
);



abip_int ABIP(validate_lin_sys)
(
	const ABIPMatrix *A
);

char *ABIP(get_lin_sys_method)
(
	spe_problem *spe
);

char *ABIP(get_lin_sys_summary)
(
	spe_problem *self,
	ABIPInfo *info
);

void ABIP(free_A_matrix)
(
	ABIPMatrix *A
);

abip_int ABIP(copy_A_matrix)
(
	ABIPMatrix **dstp, 
	const ABIPMatrix *src
);

abip_int LDL_factor(cs *A, cs **L, abip_float *Dinv);

abip_int ABIP(init_linsys_work)(spe_problem *spe);


abip_int ABIP(solve_linsys)(spe_problem *spe, abip_float *b, abip_int n, abip_float *pcg_warm_start, abip_float pcg_tol);


abip_int ABIP(free_linsys)(spe_problem *spe);
#ifdef __cplusplus
}
#endif

#endif
