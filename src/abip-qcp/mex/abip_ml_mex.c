#include "abip.h"
#include "linsys.h"
#include "matrix.h"
#include "mex.h"
#include "cones.h"


static void free_mex(ABIPData *d, ABIPCone *k) {

   /* if (k) {
        abip_free(k);
    }*/
    if (d) {

        if (d->stgs) {
            abip_free(d->stgs);
        }

        if (d->b) {
            abip_free(d->b);
        }
        if (d->c) {
            abip_free(d->c);
        }

        if (d->A) {

            if (d->A->p) {
                abip_free(d->A->p);
            }
            if (d->A->i) {
                abip_free(d->A->i);
            }
            if (d->A->x) {
                abip_free(d->A->x);
            }
            abip_free(d->A);
        }
        abip_free(d);
    }
}


/* this memory must be freed */
static abip_int *cast_to_abip_int_arr(mwIndex *arr, abip_int len) {
    abip_int i;
    abip_int *arr_out = (abip_int *)abip_malloc(sizeof(abip_int) * len);
    for (i = 0; i < len; i++) {
        arr_out[i] = (abip_int)arr[i];
    }
    return arr_out;
}


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




/*the input of matlab is data, cone, settings, output is sol, info
  data is struct containing X,y,lambda;
  matlab usage: [sol,info] = abip(data,settings)
*/
void  mexFunction(int  nlhs, mxArray* plhs[], int  nrhs, const  mxArray* prhs[])
{

    const mxArray *data;
    const mxArray *A_mex;
    const mxArray *b_mex;
    const mxArray *lambda_mex;


    const mxArray *settings;
    mxArray *tmp;

    const mwSize one[1] = {1};
    const int num_info_fields = 14;
    const char *info_fields[] = {"ipm_iter", "admm_iter", "status", "pobj", "dobj","res_pri", "res_dual", 
                                    "gap", "status_val", "setup_time", "solve_time", "runtime",  "lin_sys_time_per_iter",
                                    "avg_cg_iters"};

    const int svm_num_sol_fields = 3;
    const char *svm_sol_fields[] = {"w", "b", "xi"};

    const int lasso_num_sol_fields = 1;
    const char *lasso_sol_fields[] = {"x"};


    /* get data*/
	ABIPData* d = (ABIPData*)abip_malloc(sizeof(ABIPData));
    data = prhs[0];
    A_mex = (mxArray *)mxGetField(data, 0, "X");
    if (A_mex == ABIP_NULL) {
        abip_free(d);
        mexErrMsgTxt("ABIPData struct must contain a `X` entry.");
    }
    if (!mxIsSparse(A_mex)) {
        abip_free(d);
        mexErrMsgTxt("Input matrix X must be in sparse format (pass in sparse(X))");
    }
    b_mex = (mxArray *)mxGetField(data, 0, "y");
    if (b_mex == ABIP_NULL) {
        abip_free(d);
        mexErrMsgTxt("ABIPData struct must contain a `y` entry.");
    }
    if (mxIsSparse(b_mex)) {
        abip_free(d);
        mexErrMsgTxt("Input vector y must be in dense format (pass in full(y))");
    }
    lambda_mex = (mxArray *)mxGetField(data, 0, "lambda");
    if (lambda_mex == ABIP_NULL) {
        abip_free(d);
        mexErrMsgTxt("ABIPData struct must contain a `lambda` entry.");
    }


    d->lambda = (abip_float)*mxGetPr(lambda_mex);
    d->m = mxGetM(A_mex);
    d->n = mxGetN(A_mex);
    ABIPMatrix *A = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
    A->m = d->m;
	A->n = d->n;

    d->b = cast_to_abip_float_arr(mxGetPr(b_mex), d->m);

    A->p = cast_to_abip_int_arr(mxGetJc(A_mex), A->n + 1);
    A->i = cast_to_abip_int_arr(mxGetIr(A_mex), A->p[A->n]);
    A->x = cast_to_abip_float_arr(mxGetPr(A_mex), A->p[A->n]);

    d->A = A;
    d->c = ABIP_NULL;

	/*get settings*/
    settings = prhs[1];
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

    tmp = mxGetField(settings, 0, "prob_type");
    if (tmp != ABIP_NULL) {
        d->stgs->prob_type = (abip_int)*mxGetPr(tmp);
        if(d->stgs->prob_type != LASSO && d->stgs->prob_type != SVM && d->stgs->prob_type != SVMQP){
            mexErrMsgTxt("Invalid problem type");
        }
    }
    else{
        mexErrMsgTxt("Please input the machine learning problem type");
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


	ABIPCone* k = (ABIPCone*)abip_malloc(sizeof(ABIPCone));

    k->q = ABIP_NULL;
    k->qsize = 0;

    k->rqsize = 1;
    k->rq = (abip_int *)abip_malloc(sizeof(abip_int));
    
	k->f = 0;
    k->z = 0;

	if(d->stgs->prob_type == LASSO) {

        k->rq[0] = 2 + d->m;
        k->l = 2 * d->n; 
    } 
    else if(d->stgs->prob_type == SVM) {
        
        k->rq[0] = 2 + d->n;
        k->l = 2 + 2 * d->m + 2 * d->n;
    }
    else if(d->stgs->prob_type == SVMQP){
        k->rqsize = 0;
        k->rq = ABIP_NULL;
        k->f = d->n + 1;
        k->l = 2 * d->m;
    }
    else{
        mexErrMsgTxt("This type of machine learning problem is not supported yet");
    }

	ABIPSolution* sol = (ABIPSolution*)abip_malloc(sizeof(ABIPSolution));
	ABIPInfo* info = (ABIPInfo*)abip_malloc(sizeof(ABIPInfo));

	sol->x = ABIP_NULL;
	sol->y = ABIP_NULL;
	sol->s = ABIP_NULL;

	

	abip_int status = abip(d,sol,info,k);

    /* output sol */

    //SVM
    if(d->stgs->prob_type == 2 || d->stgs->prob_type == 4){

        plhs[0] = mxCreateStructArray(1, one, svm_num_sol_fields, svm_sol_fields);

        set_output_field(&tmp, sol->x, d->n);
        mxSetField(plhs[0], 0, "w", tmp);

        set_output_field(&tmp, sol->y, 1);
        mxSetField(plhs[0], 0, "b", tmp);

        set_output_field(&tmp, sol->s, d->m);
        mxSetField(plhs[0], 0, "xi", tmp);
    }
    //LASSO
    else{
        plhs[0] = mxCreateStructArray(1, one, lasso_num_sol_fields, lasso_sol_fields);

        set_output_field(&tmp, sol->x, d->n);
        mxSetField(plhs[0], 0, "x", tmp);
    }

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
    *mxGetPr(tmp) = (info->solve_time + info->setup_time) / 1e3;

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
