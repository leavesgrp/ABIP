#include "abip.h"
#include "linsys.h"
#include "matrix.h"
#include "mex.h"
#include "cones.h"

/*
This is mex file for creating matlab routine:

function [results, info] = abip_qp(data, cones, params)

    min       1/2x'Qx + c'x
subject to         Ax = b
                    x in K

where A \in R^m*n, b \in R^m, c \in R^n

data:
A:     m*n matrix used as above
b:     m dimensional column vector used as above
c:     n dimensional column vector used as above

K:     covex cone, now SOC, RSOC, free cone, zero cone and positive orthant are supported
*/

static void free_mex(ABIPData *d, ABIPCone *k) {

            
   if (k) {
        
        if(k->q){
            abip_free(k->q);
        }
        if(k->rq){
            abip_free(k->rq);
        }
        abip_free(k);
    }
    if (d) {

        if(d->A){

        abip_free(d->A->x);
 
        abip_free(d->A->i);
        abip_free(d->A->p);  

        }

        if (d->b) {
            abip_free(d->b);
        }
        if (d->c) {
            abip_free(d->c);
        }

        abip_free(d);
    }
}



#if !(DLONG > 0)
/* this memory must be freed */
static abip_int *cast_to_abip_int_arr(mwIndex *arr, abip_int len) {
    abip_int i;
    abip_int *arr_out = (abip_int *)abip_malloc(sizeof(abip_int) * len);
    for (i = 0; i < len; i++) {
        arr_out[i] = (abip_int)arr[i];
    }
    return arr_out;
}
#endif


/* this memory must be freed */
static abip_float *cast_to_abip_float_arr(double *arr, abip_int len) {
    abip_int i;
    abip_float *arr_out = (abip_float *)abip_malloc(sizeof(abip_float) * len);
    for (i = 0; i < len; i++) {
        arr_out[i] = (abip_float)arr[i];
    }
    return arr_out;
}

static double *cast_to_double_arr(abip_float *arr, abip_int len) {
    abip_int i;
    double *arr_out = (double *)abip_malloc(sizeof(double) * len);
    for (i = 0; i < len; i++) {
        arr_out[i] = (double)arr[i];
    }
    return arr_out;
}


static void set_output_field(mxArray **pout, abip_float *out, abip_int len) {
    *pout = mxCreateDoubleMatrix(0, 0, mxREAL);
#if SFLOAT > 0
    mxSetPr(*pout, cast_to_double_arr(out, len));
    abip_free(out);
#else
    mxSetPr(*pout, out);
#endif
    mxSetM(*pout, len);
    mxSetN(*pout, 1);
}


void  mexFunction(int  nlhs, mxArray* plhs[], int  nrhs, const  mxArray* prhs[])
{

    const mxArray *data;
    const mxArray *A_mex;
    const mxArray *Q_mex;
    const mxArray *b_mex;
    const mxArray *c_mex;
    const mxArray *cone;
    const mxArray *settings;
    const mxArray *kq;
	const mxArray *krq;
    const mxArray *kf;
    const mxArray *kz;
	const mxArray *kl;
	const double *q_mex;
  	const double *rq_mex;
	const size_t *q_dims;
  	const size_t *rq_dims;
	abip_int ns;

    mxArray *tmp;

    const mwSize one[1] = {1};
    const int num_info_fields = 14;
    const char *info_fields[] = {"ipm_iter", "admm_iter", "status", "pobj", "dobj","res_pri", "res_dual", 
                                    "gap", "status_val", "setup_time", "solve_time", "runtime", "lin_sys_time_per_iter",
                                    "avg_cg_iters"};

    const int num_sol_fields = 3;
    const char *sol_fields[] = {"x", "y", "s"};


    /* get data*/
	ABIPData* d = (ABIPData*)abip_malloc(sizeof(ABIPData));
    data = prhs[0];


    b_mex = (mxArray *)mxGetField(data, 0, "b");
    if (b_mex == ABIP_NULL) {
        d->b = ABIP_NULL;
        printf("ABIPData doesn't contain `b`\n");
    }
    else if (mxIsSparse(b_mex)) {
        abip_free(d);
        mexErrMsgTxt("Input vector b must be in dense format (pass in full(b))");
    }
    else{

        d->b = cast_to_abip_float_arr(mxGetPr(b_mex), MAX(mxGetM(b_mex), mxGetN(b_mex)));

    }

    c_mex = (mxArray *)mxGetField(data, 0, "c");
    if (c_mex == ABIP_NULL) {
        abip_free(d);
        mexErrMsgTxt("ABIPData struct must contain a `c` entry.");
    }
    else if (mxIsSparse(c_mex)) {
        abip_free(d);
        mexErrMsgTxt("Input vector c must be in dense format (pass in full(c))");
    }
    else{

        d->c = cast_to_abip_float_arr(mxGetPr(c_mex), MAX(mxGetM(c_mex), mxGetN(c_mex)));

    }
    
    A_mex = (mxArray *)mxGetField(data, 0, "A");
    if (A_mex == ABIP_NULL) {
        d->A = ABIP_NULL;
        d->m = 0;
        d->n = MAX(mxGetM(c_mex), mxGetN(c_mex));

        printf("ABIPData doesn't contain `A`\n");
    }
    else if (!mxIsSparse(A_mex)) {
        abip_free(d);
        mexErrMsgTxt("Input matrix A must be in sparse format (pass in sparse(A))");
    }
    else{
        d->A = (ABIPMatrix*)abip_malloc(sizeof(ABIPMatrix));
        ABIPMatrix *A = d->A;

        A->m = mxGetM(A_mex);
        A->n = mxGetN(A_mex);
        
        d->m = A->m;
        d->n = A->n;


        A->p = cast_to_abip_int_arr(mxGetJc(A_mex), A->n + 1);
        A->i = cast_to_abip_int_arr(mxGetIr(A_mex), A->p[A->n]);

        A->x = cast_to_abip_float_arr(mxGetPr(A_mex), A->p[A->n]);

    }
    
    Q_mex = (mxArray *)mxGetField(data, 0, "Q");
    if (Q_mex == ABIP_NULL){
        d->Q = ABIP_NULL;
        printf("ABIPData doesn't contain `Q`\n");
    }
    else if (!mxIsSparse(Q_mex)) {
        abip_free(d);
        mexErrMsgTxt("Input matrix Q must be in sparse format (pass in sparse(Q))");
    }
    else{
        d->Q = (ABIPMatrix*)abip_malloc(sizeof(ABIPMatrix));
        ABIPMatrix *Q = d->Q;

        Q->m = mxGetM(Q_mex);
        Q->n = mxGetN(Q_mex);


        Q->p = cast_to_abip_int_arr(mxGetJc(Q_mex), Q->n + 1);
        Q->i = cast_to_abip_int_arr(mxGetIr(Q_mex), Q->p[Q->n]);
        Q->x = cast_to_abip_float_arr(mxGetPr(Q_mex), Q->p[Q->n]);

    }

    /*get cone*/
    cone = prhs[1];

    ABIPCone *k = (ABIPCone*)abip_malloc(sizeof(ABIPCone));

	kq = mxGetField(cone, 0, "q");
	if (kq && !mxIsEmpty(kq)) {
		q_mex = mxGetPr(kq);
		ns = (abip_int)mxGetNumberOfDimensions(kq);
		q_dims = mxGetDimensions(kq);
		k->qsize = (abip_int)q_dims[0];
		if (ns > 1 && q_dims[0] == 1) {
		k->qsize = (abip_int)q_dims[1];
		}
		k->q = (abip_int *)mxMalloc(sizeof(abip_int) * k->qsize);
		for (abip_int i = 0; i < k->qsize; i++) {
		k->q[i] = (abip_int)q_mex[i];
		}
	} else {
		k->qsize = 0;
		k->q = ABIP_NULL;
	}

	krq = mxGetField(cone, 0, "rq");
	if (krq && !mxIsEmpty(krq)) {
		rq_mex = mxGetPr(krq);
		ns = (abip_int)mxGetNumberOfDimensions(krq);
		rq_dims = mxGetDimensions(krq);
		k->rqsize = (abip_int)rq_dims[0];
		if (ns > 1 && rq_dims[0] == 1) {
		k->rqsize = (abip_int)rq_dims[1];
		}
		k->rq = (abip_int *)mxMalloc(sizeof(abip_int) * k->rqsize);
		for (abip_int i = 0; i < k->rqsize; i++) {
		k->rq[i] = (abip_int)rq_mex[i];
		}
	} else {
		k->rqsize = 0;
		k->rq = ABIP_NULL;
	}

    kf = mxGetField(cone, 0, "f");
	if (kf && !mxIsEmpty(kf)) {
		k->f = (abip_int)*mxGetPr(kf);
	} else {
		k->f = 0;
	}

    kz = mxGetField(cone, 0, "z");
	if (kz && !mxIsEmpty(kz)) {
		k->z = (abip_int)*mxGetPr(kz);
	} else {
		k->z = 0;
	}

	kl = mxGetField(cone, 0, "l");
	if (kl && !mxIsEmpty(kl)) {
		k->l = (abip_int)*mxGetPr(kl);
	} else {
		k->l = 0;
	}

	/*get settings*/
    settings = prhs[2];
	d->stgs = (ABIPSettings*)abip_malloc(sizeof(ABIPSettings));
    ABIP(set_default_settings)(d);

    tmp = mxGetField(settings, 0, "alpha");
    if (tmp != ABIP_NULL) {
        d->stgs->alpha = (abip_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "cg_rate");
    if (tmp != ABIP_NULL) {
        d->stgs->cg_rate = (abip_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "eps");
    if (tmp != ABIP_NULL) {
        d->stgs->eps = (abip_float)*mxGetPr(tmp);
        d->stgs->eps_p = d->stgs->eps;
        d->stgs->eps_d = d->stgs->eps;
        d->stgs->eps_g = d->stgs->eps;
        d->stgs->eps_inf = d->stgs->eps;
        d->stgs->eps_unb = d->stgs->eps;
    }
    tmp = mxGetField(settings, 0, "eps_p");
    if (tmp != ABIP_NULL) {
        d->stgs->eps_p = (abip_float)*mxGetPr(tmp);
    }
    tmp = mxGetField(settings, 0, "eps_d");
    if (tmp != ABIP_NULL) {
        d->stgs->eps_d = (abip_float)*mxGetPr(tmp);
    }
    tmp = mxGetField(settings, 0, "eps_g");
    if (tmp != ABIP_NULL) {
        d->stgs->eps_g = (abip_float)*mxGetPr(tmp);
    }
        tmp = mxGetField(settings, 0, "eps_inf");
    if (tmp != ABIP_NULL) {
        d->stgs->eps_inf = (abip_float)*mxGetPr(tmp);
    }    tmp = mxGetField(settings, 0, "eps_unb");
    if (tmp != ABIP_NULL) {
        d->stgs->eps_unb = (abip_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "max_admm_iters");
    if (tmp != ABIP_NULL) {
        d->stgs->max_admm_iters = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "max_ipm_iters");
    if (tmp != ABIP_NULL) {
        d->stgs->max_ipm_iters = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "normalize");
    if (tmp != ABIP_NULL) {
        d->stgs->normalize = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "rho_y");
    if (tmp != ABIP_NULL) {
        d->stgs->rho_y = (abip_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "rho_x");
    if (tmp != ABIP_NULL) {
        d->stgs->rho_x = (abip_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "rho_tau");
    if (tmp != ABIP_NULL) {
        d->stgs->rho_tau = (abip_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "scale");
    if (tmp != ABIP_NULL) {
        d->stgs->scale = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "scale_bc");
    if (tmp != ABIP_NULL) {
        d->stgs->scale_bc = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "scale_E");
    if (tmp != ABIP_NULL) {
        d->stgs->scale_E = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "use_indirect");
    if (tmp != ABIP_NULL) {
        d->stgs->use_indirect = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "verbose");
    if (tmp != ABIP_NULL) {
        d->stgs->verbose = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "linsys_solver");
    if (tmp != ABIP_NULL) {
        d->stgs->linsys_solver = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "inner_check_period");
    if (tmp != ABIP_NULL) {
        d->stgs->inner_check_period = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "outer_check_period");
    if (tmp != ABIP_NULL) {
        d->stgs->outer_check_period = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "err_dif");
    if (tmp != ABIP_NULL) {
        d->stgs->err_dif = (abip_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "time_limit");
    if (tmp != ABIP_NULL) {
        d->stgs->time_limit = (abip_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "psi");
    if (tmp != ABIP_NULL) {
        d->stgs->psi = (abip_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "origin_scaling");
    if (tmp != ABIP_NULL) {
        d->stgs->origin_scaling = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "ruiz_scaling");
    if (tmp != ABIP_NULL) {
        d->stgs->ruiz_scaling = (abip_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(settings, 0, "pc_scaling");
    if (tmp != ABIP_NULL) {
        d->stgs->pc_scaling = (abip_int)*mxGetPr(tmp);
    }

    d->stgs->prob_type = QCP; //QCP

	ABIPSolution* sol = (ABIPSolution*)abip_malloc(sizeof(ABIPSolution));
	ABIPInfo* info = (ABIPInfo*)abip_malloc(sizeof(ABIPInfo));

	sol->x = ABIP_NULL;
	sol->y = ABIP_NULL;
	sol->s = ABIP_NULL;

	abip_int status = abip(d,sol,info,k);

    /* output sol */
    plhs[0] = mxCreateStructArray(1, one, num_sol_fields, sol_fields);

    set_output_field(&tmp, sol->x, d->n);
    mxSetField(plhs[0], 0, "x", tmp);

    set_output_field(&tmp, sol->y, d->m);
    mxSetField(plhs[0], 0, "y", tmp);

    set_output_field(&tmp, sol->s, d->n);
    mxSetField(plhs[0], 0, "s", tmp);

    /* output info */
    plhs[1] = mxCreateStructArray(1, one, num_info_fields, info_fields);

    /* if you add/remove fields here update the info_fields above */
    mxSetField(plhs[1], 0, "status", mxCreateString(info->status));

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "ipm_iter", tmp);
    *mxGetPr(tmp) = (abip_float)info->ipm_iter;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "admm_iter", tmp);
    *mxGetPr(tmp) = (abip_float)info->admm_iter;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "status_val", tmp);
    *mxGetPr(tmp) = (abip_float)info->status_val;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "pobj", tmp);
    *mxGetPr(tmp) = info->pobj;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "dobj", tmp);
    *mxGetPr(tmp) = info->dobj;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "res_pri", tmp);
    *mxGetPr(tmp) = info->res_pri;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "res_dual", tmp);
    *mxGetPr(tmp) = info->res_dual;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "gap", tmp);
    *mxGetPr(tmp) = info->rel_gap;

    /*return value in secs */
    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "setup_time", tmp);
    *mxGetPr(tmp) = info->setup_time / 1e3;

    /*return value in secs */
    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "solve_time", tmp);
    *mxGetPr(tmp) = info->solve_time / 1e3;

    /*return value in secs */
    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "runtime", tmp);
    *mxGetPr(tmp) = (info->setup_time + info->solve_time) / 1e3;

    /*return value in secs */
    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "lin_sys_time_per_iter", tmp);
    *mxGetPr(tmp) = info->avg_linsys_time / 1e3;
     
    //average cg iters per admm iter
    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[1], 0, "avg_cg_iters", tmp);
    *mxGetPr(tmp) = info->avg_cg_iters;


    free_mex(d, k);
    return;
}