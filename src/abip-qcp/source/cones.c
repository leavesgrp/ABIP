#define _CRT_SECURE_NO_WARNINGS
#include "cones.h"



/* c = l + f + q[i] + s[i] */
abip_int ABIP(get_cone_dims)(const ABIPCone* k) {
    abip_int i, c = 0;
    
    if (k->l) {
        c += k->l;
    }
    if (k->z) {
        c += k->z;
    }
    if (k->f) {
        c += k->f;
    }
    if (k->qsize && k->q) {
        for (i = 0; i < k->qsize; ++i) {
            c += k->q[i];
        }
    }
    if (k->rqsize && k->rq) {
        for (i = 0; i < k->rqsize; ++i) {
            c += k->rq[i];
        }
    }

    return c;
}


abip_int ABIP(validate_cones)(spe_problem *spe, const ABIPCone* k) {
    abip_int i;
    if (ABIP(get_cone_dims)(k) != spe->q) {
        abip_printf("cone dimensions %li not equal to num rows in A = n = %li\n",
            (long)ABIP(get_cone_dims)(k), (long)spe->q);
        return -1;
    }
    if (k->l && k->l < 0) {
        abip_printf("lp cone error\n");
        return -1;
    }
    if (k->f && k->f < 0) {
        abip_printf("free cone error\n");
        return -1;
    }
    if (k->z && k->z < 0) {
        abip_printf("zero cone error\n");
        return -1;
    }
    if (k->qsize && k->q) {
        if (k->qsize < 0) {
            abip_printf("soc cone error\n");
            return -1;
        }
        for (i = 0; i < k->qsize; ++i) {
            if (k->q[i] < 0) {
                abip_printf("soc cone error\n");
                return -1;
            }
        }
    }
    if (k->rqsize && k->rq) {
        if (k->rqsize < 0) {
            abip_printf("rsoc cone error\n");
            return -1;
        }
        for (i = 0; i < k->rqsize; ++i) {
            if (k->rq[i] < 0) {
                abip_printf("rsoc cone error\n");
                return -1;
            }
        }
    }
   
    
    return 0;
}

char* ABIP(get_cone_header)(const ABIPCone* k) {

    abip_int i, soc_vars, rsoc_vars;
    char* tmp = (char*)abip_malloc(sizeof(char) * 512);
    sprintf(tmp, "Cones:");
   
    
    soc_vars = 0;
    if (k->qsize && k->q) {
        for (i = 0; i < k->qsize; i++) {
            soc_vars += k->q[i];
        }
        sprintf(tmp + strlen(tmp), "\tsoc vars: %li, soc blks: %li\n",
            (long)soc_vars, (long)k->qsize);
    }

    rsoc_vars = 0;
    if (k->rqsize && k->rq) {
        for (i = 0; i < k->rqsize; i++) {
            rsoc_vars += k->rq[i];
        }
        sprintf(tmp + strlen(tmp), "\trsoc vars: %li, rsoc blks: %li\n",
            (long)rsoc_vars, (long)k->rqsize);
    }

    if (k->f) {
        sprintf(tmp + strlen(tmp), "\tfree vars: %li\n", (long)k->f);
    }

    if (k->z) {
        sprintf(tmp + strlen(tmp), "\tzero vars: %li\n", (long)k->z);
    }

    if (k->l) {
        sprintf(tmp + strlen(tmp), "\tlinear vars: %li\n", (long)k->l);
    }
    return tmp;
}



/*  x in second order cone 
    K = {(t,x)| t>||x||}
*/
void ABIP(soc_barrier_subproblem)(
    abip_float* x,
    abip_float* tmp,
    abip_float lambda,
    abip_int n
 )
{
    abip_float a = tmp[0];
    abip_float tol = 1e-9;
    abip_float* b = (abip_float*)abip_malloc((n - 1) * sizeof(abip_float));
    memcpy(b, &tmp[1], (n - 1) * sizeof(abip_float));
    abip_float b_norm_sq = ABIP(norm_sq)(b, n - 1);
    if(ABS(a)<=tol){
        x[0] = SQRTF(2*lambda + b_norm_sq/4);
        memcpy(&x[1], b, (n-1) * sizeof(abip_float));
        ABIP(scale_array)(&x[1], 0.5, (n-1));
    }
    else{
        abip_float coef1 = -(POWF(a, 2) - ABIP(norm_sq)(b, n - 1)) / lambda;
        abip_float coef2 = -(4 * POWF(a, 2) + 4 * ABIP(norm_sq)(b, n - 1) + 16 * lambda) / lambda;

        abip_float r = 16*a*a/(8*lambda-a*a+b_norm_sq+SQRTF(POWF((8*lambda-a*a+b_norm_sq),2)+32*a*a*lambda));
        abip_float s1 = (r-SQRTF(r*(r+8)))/2;
        abip_float s2 = (r+SQRTF(r*(r+8)))/2;

        abip_float s = a > 0 ? s2 : s1;
        abip_float eta = (s+2) * a / s;
        ABIP(scale_array)(b, (s+2) / (s+4), n - 1);
        
        x[0] = eta;
        memcpy(&x[1], b, (n - 1) * sizeof(abip_float));
    }
    abip_free(b);
}


/*  K = {(t1,t2,x)| 2*t1*t2>||x||}
    n = len(t1,t2,x)
*/
void ABIP(rsoc_barrier_subproblem)(
    abip_float* x,
    abip_float* tmp,
    abip_float lambda,
    abip_int n 
)
{   
    abip_int nx = n-2;
    abip_float tol = 1e-9;

    abip_float zeta_eta = tmp[0];
    abip_float zeta_nu = tmp[1];

    abip_float* zeta_x = (abip_float*)abip_malloc(nx * sizeof(abip_float));
    memcpy(zeta_x, &tmp[2], nx * sizeof(abip_float));
    abip_float zeta_x_norm_sq = ABIP(norm_sq)(zeta_x, nx);

    if(zeta_eta+zeta_nu==0){
        x[1] = (-zeta_eta + SQRTF(zeta_eta*zeta_eta + 4*lambda + zeta_x_norm_sq )) / 2;
        x[0] = x[0] + zeta_eta;
        memcpy(&x[2], zeta_x, nx * sizeof(abip_float));
        ABIP(scale_array)(&x[2], 0.5, nx);
    }
    else{
        abip_float w;
        abip_float s;
        if(2*zeta_eta*zeta_nu - zeta_x_norm_sq < 0){
            w = (2 * POWF(zeta_eta+zeta_nu,2) / lambda) / (-(2*zeta_eta*zeta_nu - zeta_x_norm_sq)/(2*lambda)) / 
            (1 + 4/(-(2*zeta_eta*zeta_nu - zeta_x_norm_sq)/(2*lambda)) + SQRTF(1 + (4*(zeta_eta*zeta_eta+zeta_nu*zeta_nu+zeta_x_norm_sq)/lambda + 16) 
            / (-(2*zeta_eta*zeta_nu - zeta_x_norm_sq)/(2*lambda)) / (-(2*zeta_eta*zeta_nu - zeta_x_norm_sq)/(2*lambda))));
        }
        else{
            w = (2 * zeta_eta*zeta_nu - zeta_x_norm_sq) / (2*lambda) * (1 - 4/((2*zeta_eta*zeta_nu - zeta_x_norm_sq)/(2*lambda)) + 
            SQRTF(1 + (4*(zeta_eta*zeta_eta+zeta_nu*zeta_nu+zeta_x_norm_sq)/lambda + 16) / ((2*zeta_eta*zeta_nu - zeta_x_norm_sq)/(2*lambda)) / 
            ((2*zeta_eta*zeta_nu - zeta_x_norm_sq)/(2*lambda)))) / 2;
        }
        if(zeta_eta + zeta_nu > 0){
            s = (w + SQRTF(w*(w+4))) / 2;
            x[0] = (zeta_eta*POWF(s+1,2) + zeta_nu * (s+1)) / (s*(s+2));
            x[1] = (zeta_nu*POWF(s+1,2) + zeta_eta * (s+1)) / (s*(s+2));
            memcpy(&x[2], zeta_x, nx * sizeof(abip_float));
            ABIP(scale_array)(&x[2], (s+1)/(s+2), nx);
        }
        else{
            if(w > 10){
                s = 2 / (w + 2 + SQRTF(w * (w+4)));
                x[0] = (zeta_eta * POWF(s,2) + zeta_nu * s) / ((s-1)*(s+1));
                x[1] = (zeta_nu * POWF(s,2) + zeta_eta * s) / ((s-1)*(s+1));
                memcpy(&x[2], zeta_x, nx * sizeof(abip_float));
                ABIP(scale_array)(&x[2], s/(s+1), nx);
            }
            else{
                s = (w - SQRTF(w*(w+4))) / 2;
                x[0] = (zeta_eta*POWF(s+1,2) + zeta_nu * (s+1)) / (s*(s+2));
                x[1] = (zeta_nu*POWF(s+1,2) + zeta_eta * (s+1)) / (s*(s+2));
                memcpy(&x[2], zeta_x, nx * sizeof(abip_float));
                ABIP(scale_array)(&x[2], (s+1)/(s+2), nx);
            }
        }
    }

    for(int i=nx+2;i<n;i++){
        if(tmp[i]>=0){
            x[i] = (tmp[i] + SQRTF(tmp[i]*tmp[i] + 4 * lambda)) / 2;
        }
        else{
            x[i] = 2*lambda / (-tmp[i]*(1+SQRTF(1+4*lambda/POWF(tmp[i],2))));
        }
    }
    abip_free(zeta_x);
}


/*  K = {x | x free}
*/
void ABIP(free_barrier_subproblem)(
    abip_float* x,
    abip_float* tmp,
    abip_float lambda,
    abip_int n 
)
{   
    
    for(int i=0;i<n;i++){
        x[i] = tmp[i];
    }
}

/*  K = {x | x = 0}
*/
void ABIP(zero_barrier_subproblem)(
    abip_float* x,
    abip_float* tmp,
    abip_float lambda,
    abip_int n 
)
{   
    
    for(int i=0;i<n;i++){
        x[i] = 0;
    }
}


/*  K = {x | x >= 0}
*/
void ABIP(positive_orthant_barrier_subproblem)(
    abip_float* x,
    abip_float* tmp,
    abip_float lambda,
    abip_int n 
)
{   
    
    for(int i=0;i<n;i++){
        if(tmp[i]>=0){
            x[i] = (tmp[i] + SQRTF(tmp[i]*tmp[i] + 4 * lambda)) / 2;
        }
        else{
            x[i] = 2*lambda / (-tmp[i]*(1+SQRTF(1+4*lambda/POWF(tmp[i],2))));
        }
    }
}