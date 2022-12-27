#include "svm_qp_config.h"
#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)

abip_int init_svmqp(svmqp **self, ABIPData *d, ABIPSettings *stgs){

    svmqp *this_svm = (svmqp*)abip_malloc(sizeof(svmqp));
    *self = this_svm;

    this_svm->pro_type = SVMQP;
    this_svm->m = d->m;
    this_svm->n = d->n;
    this_svm->p = d->m;
    this_svm->q = 1 + d->n + 2*d->m;
    this_svm->lambda = d->lambda;
    abip_int m = this_svm->p;
    abip_int n = this_svm->q;

    this_svm->Q = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
    this_svm->Q->m = n;
    this_svm->Q->n = n;
    this_svm->Q->i = (abip_int*)abip_malloc(sizeof(abip_int) * d->n);
    this_svm->Q->p = (abip_int*)abip_malloc(sizeof(abip_int) * (n+1));
    this_svm->Q->x = (abip_float*)abip_malloc(sizeof(abip_float) * d->n);

    for(abip_int i=0; i<d->n; i++){
        this_svm->Q->i[i] = i;
        this_svm->Q->p[i] = i;
        this_svm->Q->x[i] = 1;
    }

    for(abip_int i=0; i<2*d->m + 2; i++){
        this_svm->Q->p[d->n + i] = d->n;
    }

    this_svm->sparsity = (((abip_float)d->A->p[d->n] / (d->m*d->n)) < 0.05);
    this_svm->rho_dr = (abip_float *)abip_malloc((m + n + 1) * sizeof(abip_float));
    for(int i=0;i<m+n+1;i++){
        if(i < m){
            this_svm->rho_dr[i] = stgs->rho_y;
        }
        else if(i < m +n){
            this_svm->rho_dr[i] = stgs->rho_x;
        }
        else{
            this_svm->rho_dr[i] = stgs->rho_tau;
        }
    }

  
    this_svm->L = (ABIPLinSysWork *)abip_malloc(sizeof(ABIPLinSysWork));
    this_svm->stgs = stgs;

    this_svm->data = (ABIPData*)abip_malloc(sizeof(ABIPData));

    abip_float *data_b = (abip_float*)abip_malloc(m*sizeof(abip_float));
    memset(data_b,0,m*sizeof(abip_float));
    for(int i=0;i<m;i++){
        data_b[i] = 1;
    }
    this_svm->data->b = data_b;

    abip_float *data_c = (abip_float*)abip_malloc(n*sizeof(abip_float));
    memset(data_c,0,this_svm->q*sizeof(abip_float));
    
    for(int i=0;i<d->m;i++){
        data_c[i+d->n+1] = 1.0 / (this_svm->m * this_svm->lambda);
    }
    this_svm->data->c = data_c;

    d->c = (abip_float*)abip_malloc(n*sizeof(abip_float));
    memcpy(d->c, data_c, n*sizeof(abip_float));

    ABIPMatrix *data_A = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
    data_A->m = d->m;
    data_A->n = d->n + 1;
    data_A->p = (abip_int*)abip_malloc((data_A->n + 1)*sizeof(abip_int));
    data_A->i = (abip_int*)abip_malloc((d->A->p[d->n] + d->m)*sizeof(abip_int));
    data_A->x = (abip_float*)abip_malloc((d->A->p[d->n] + d->m)*sizeof(abip_float));
    for(int i=0;i<d->A->p[d->n];i++){
        d->A->x[i] *= d->b[d->A->i[i]];
    }
    memcpy(data_A->p,d->A->p,(d->n + 1)*sizeof(abip_int));
    data_A->p[data_A->n] = d->A->p[d->n] + d->m;

    memcpy(data_A->i,d->A->i,d->A->p[d->n]*sizeof(abip_int));
    for(int i=d->A->p[d->n];i<d->A->p[d->n]+d->m;i++){
        data_A->i[i] = i - d->A->p[d->n];
    }

    memcpy(data_A->x,d->A->x,d->A->p[d->n]*sizeof(abip_float));
    for(int i=d->A->p[d->n];i<d->A->p[d->n]+d->m;i++){
        data_A->x[i] = d->b[i - d->A->p[d->n]];
    }
    this_svm->data->A = data_A;
    this_svm->A = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
    this_svm->b = (abip_float *)abip_malloc(this_svm->p*sizeof(abip_float));
    this_svm->c = (abip_float *)abip_malloc(this_svm->q*sizeof(abip_float));
    this_svm->D = (abip_float*)abip_malloc(this_svm->p*sizeof(abip_float));

    this_svm->E = (abip_float*)abip_malloc(this_svm->q*sizeof(abip_float));
    for(int i=0;i<this_svm->q;i++){
        this_svm->E[i] = 1.0;
    }

    this_svm->F = (abip_float*)abip_malloc(this_svm->m*sizeof(abip_float));
    this_svm->H = (abip_float*)abip_malloc(this_svm->q*sizeof(abip_float));


    this_svm->scaling_data = &scaling_svmqp_data;
    this_svm->un_scaling_sol = &un_scaling_svmqp_sol;
    this_svm->calc_residuals = &calc_svmqp_residuals;
    this_svm->init_spe_linsys_work = &init_svmqp_linsys_work;
    this_svm->solve_spe_linsys = &solve_svmqp_linsys;
    this_svm->free_spe_linsys_work = &free_svmqp_linsys_work;
    this_svm->spe_A_times = &svmqp_A_times;
    this_svm->spe_AT_times = &svmqp_AT_times;
    this_svm->inner_conv_check = &svmqp_inner_conv_check;

    return 0;
}

void svmqp_A_times(svmqp *self, const abip_float *x, abip_float* y){

    ABIP(accum_by_A)(self->A, x, y);

    for(int i=0; i<self->m; i++){
        y[i] += 1/self->D[i] * (x[self->n + 1 + i] - x[self->n + 1 + self->m + i]);
    }
}


void svmqp_AT_times(svmqp *self, const abip_float *x, abip_float* y){

    ABIP(accum_by_Atrans)(self->A, x, y);
    for(int i=0; i<self->m; i++){
        y[self->n + 1 + i] += 1/self->D[i] * x[i];
        y[self->n + 1 + self->m + i] -= 1/self->D[i] * x[i];
    }

}


abip_float svmqp_inner_conv_check(svmqp *self, ABIPWork *w){

    abip_int m = self->p;
    abip_int n = self->q;

    abip_float *Qu = (abip_float*)abip_malloc((m+n+1)*sizeof(abip_float));
    abip_float *Mu = (abip_float*)abip_malloc((m+n)*sizeof(abip_float));
    
    memset(Mu, 0, (m+n)*sizeof(abip_float));

    self->spe_A_times(self, &w->u[m], Mu);

    self->spe_AT_times(self, w->u, &Mu[m]);
    ABIP(scale_array)(&Mu[m], -1, n);

    if(self->Q != ABIP_NULL){
        ABIP(accum_by_A)(self->Q, &w->u[m], &Mu[m]);
    }

    memcpy(Qu, Mu, (m+n)*sizeof(abip_float));

    ABIP(add_scaled_array)(Qu, self->b, m, -w->u[m+n]);
    ABIP(add_scaled_array)(&Qu[m], self->c, n, w->u[m+n]);

    Qu[m+n] = -ABIP(dot)(w->u, Mu, m+n) / w->u[m+n] + ABIP(dot)(w->u, self->b, m) - ABIP(dot)(&w->u[m], self->c, n);


    abip_float *tem = (abip_float*)abip_malloc((m+n+1)*sizeof(abip_float));
    memcpy(tem, Qu, (m+n+1)*sizeof(abip_float));
    ABIP(add_scaled_array)(tem, w->v_origin, m+n+1, -1);

    abip_float error_inner = ABIP(norm)(tem, m+n+1) / (1 + ABIP(norm)(Qu, m+n+1) + ABIP(norm)(w->v_origin, m+n+1));

    abip_free(Qu);
    abip_free(Mu);
    abip_free(tem);

    return error_inner;
}


void scaling_svmqp_data(svmqp *self, ABIPCone *k){


    memcpy(self->b, self->data->b, self->p * sizeof(abip_float));
    
    memcpy(self->c, self->data->c, self->q * sizeof(abip_float));

    ABIP(copy_A_matrix)(&(self->A), self->data->A);

    abip_int m = self->p;
    abip_int n = self->n + 1;
    ABIPMatrix *A = self->A;
    ABIPMatrix *Q = self->Q;

    abip_float min_row_scale = MIN_SCALE * SQRTF((abip_float)n); 
    abip_float max_row_scale = MAX_SCALE * SQRTF((abip_float)n);
    abip_float min_col_scale = MIN_SCALE * SQRTF((abip_float)m); 
    abip_float max_col_scale = MAX_SCALE * SQRTF((abip_float)m);
    
    abip_float *E_hat = self->E;
    abip_float *D_hat = self->D;

    for(int i=0;i<n;i++){
        E_hat[i] = 1;
    }
    for(int i=0;i<m;i++){
        D_hat[i] = 1;
    }

    abip_float *E = (abip_float*)abip_malloc(n*sizeof(abip_float));
    memset(E,0,n*sizeof(abip_float));

    abip_float *E1 = (abip_float*)abip_malloc(n*sizeof(abip_float));
    memset(E1,0,n*sizeof(abip_float));

    abip_float *E2 = (abip_float*)abip_malloc(n*sizeof(abip_float));
    memset(E2,0,n*sizeof(abip_float));


    abip_float *D = (abip_float*)abip_malloc(m*sizeof(abip_float));
    memset(D,0,m*sizeof(abip_float));

    abip_int origin_scaling = self->stgs->origin_scaling;
    abip_int ruiz_scaling = self->stgs->ruiz_scaling;
    abip_int pc_scaling = self->stgs->pc_scaling;
    abip_int count;
    abip_float mean_E;


    if(A == ABIP_NULL && Q == ABIP_NULL){
        origin_scaling = 0;
        ruiz_scaling = 0;
        pc_scaling = 0;
    }

    if(ruiz_scaling){

        abip_int n_ruiz = 10;

        for(int ruiz_iter=0;ruiz_iter<n_ruiz;ruiz_iter++){

            count = 0;
            memset(E,0,n*sizeof(abip_float));
            memset(E1,0,n*sizeof(abip_float));
            memset(E2,0,n*sizeof(abip_float));
            memset(D,0,m*sizeof(abip_float));

            if(A != ABIP_NULL){
                for(int j=0;j<n;j++){
                    if(A->p[j] == A->p[j+1]){
                        E1[j] = 0;
                    }
                    else{
                        E1[j] = SQRTF(ABIP(norm_inf)(&A->x[A->p[j]], A->p[j+1] - A->p[j]));
                    }
                }
            }

            if(Q != ABIP_NULL){
                for(int j=0;j<n;j++){
                    if(Q->p[j] == Q->p[j+1]){
                        E2[j] = 0;
                    }
                    else{
                        E2[j] = SQRTF(ABIP(norm_inf)(&Q->x[Q->p[j]], Q->p[j+1] - Q->p[j]));
                    }
                }
            }

            for (int i = 0; i < n; i++){
                E[i] = E1[i] < E2[i] ? E2[i] : E1[i];
            }

            if(k->q){
                for(int i=0;i<k->qsize;i++){
                    mean_E = ABIP(vec_mean)(&E[count],k->q[i]);
                    for(int j=0;j<k->q[i];j++){
                        E[j+count] = mean_E;
                    }
                    count += k->q[i];
                }
            }

            if(k->rq){
                for(int i=0;i<k->rqsize;i++){
                    mean_E = ABIP(vec_mean)(&E[count],k->rq[i]);
                    for(int j=0;j<k->rq[i];j++){
                        E[j+count] = mean_E;
                    }
                    count += k->rq[i];
                }
            }

            if(A != ABIP_NULL){
                for(int i=0;i<A->p[n];i++){
                    if(D[A->i[i]] < ABS(A->x[i])){
                        D[A->i[i]] = ABS(A->x[i]);
                    }
                }
                for(int i=0;i<m;i++){
                    D[i] = SQRTF(D[i]);
                    if (D[i] < min_row_scale) D[i] = 1;
                    else if (D[i] > max_row_scale) D[i] = max_row_scale;
                }

                for(int i=0;i<n;i++){
                    if (E[i] < min_col_scale) E[i] = 1;
                    else if (E[i] > max_col_scale) E[i] = max_col_scale;
                    for(int j=A->p[i];j<A->p[i+1];j++){
                        A->x[j] /= E[i];
                    }
                }
            }

            if(Q != ABIP_NULL){

                for(int i=0;i<n;i++){
                    for(int j=Q->p[i];j<Q->p[i+1];j++){
                        Q->x[j] /= E[i];
                    }
                }
                for(int i=0;i<Q->p[n];i++){
                    Q->x[i] /= E[Q->i[i]];
                }

            }
        
            if(A != ABIP_NULL){         
                for(int i=0;i<A->p[n];i++){
                    A->x[i] /= D[A->i[i]];
                }
            }

            for(int i=0;i<n;i++){
                E_hat[i] *= E[i];
            }

            for(int i=0;i<m;i++){
                D_hat[i] *= D[i];
            }

        }
    }

    if(origin_scaling){
        
        memset(E,0,n*sizeof(abip_float));
        memset(E1,0,n*sizeof(abip_float));
        memset(E2,0,n*sizeof(abip_float));
        memset(D,0,m*sizeof(abip_float));

        count = 0;

        if(A != ABIP_NULL){         
            for (int i = 0; i < n; i++) {
                for (int j = A->p[i]; j < A->p[i + 1]; j++) {
                    E1[i] += A->x[j] * A->x[j];
                }
                E1[i] = SQRTF(E1[i]);
            }
        }

        if(Q != ABIP_NULL){

            for (int i = 0; i < n; i++) {
                for (int j = Q->p[i]; j < Q->p[i + 1]; j++) {
                    E2[i] += Q->x[j] * Q->x[j];
                }
                E2[i] = SQRTF(E2[i]);
            }

        }

        for (int i = 0; i < n; i++){
            E[i] = SQRTF(E1[i] < E2[i] ? E2[i] : E1[i]);
        }


        if(k->q){
            for(int i=0;i<k->qsize;i++){
                mean_E = ABIP(vec_mean)(&E[count],k->q[i]);
                for(int j=0;j<k->q[i];j++){
                    E[j+count] = mean_E;
                }
                count += k->q[i];
            }
        }

        if(k->rq){
            for(int i=0;i<k->rqsize;i++){
                mean_E = ABIP(vec_mean)(&E[count],k->rq[i]);
                for(int j=0;j<k->rq[i];j++){
                    E[j+count] = mean_E;
                }
                count += k->rq[i];
            }
        }

        if(A != ABIP_NULL){         
            for(int i=0;i<A->p[n];i++){
                D[A->i[i]] += A->x[i] * A->x[i];
            }
            for(int i=0;i<m;i++){
                D[i] = SQRTF(SQRTF(D[i]));
                if (D[i] < min_row_scale) D[i] = 1;
                else if (D[i] > max_row_scale) D[i] = max_row_scale;
            }

            for(int i=0;i<n;i++){
                if (E[i] < min_col_scale) E[i] = 1;
                else if (E[i] > max_col_scale) E[i] = max_col_scale;
                for(int j=A->p[i];j<A->p[i+1];j++){
                    A->x[j] /= E[i];
                }
            }
        }

        if(Q != ABIP_NULL){

            for(int i=0;i<n;i++){
                for(int j=Q->p[i];j<Q->p[i+1];j++){
                    Q->x[j] /= E[i];
                }
            }
            for(int i=0;i<Q->p[n];i++){
                Q->x[i] /= E[Q->i[i]];
            }
        }
    
        if(A != ABIP_NULL){         
            for(int i=0;i<A->p[n];i++){
                A->x[i] /= D[A->i[i]];
            }
        }

        for(int i=0;i<n;i++){
            E_hat[i] *= E[i];
        }

        for(int i=0;i<m;i++){
            D_hat[i] *= D[i];
        }
    }

    if(pc_scaling){

        memset(E,0,n*sizeof(abip_float));
        memset(E1,0,n*sizeof(abip_float));
        memset(E2,0,n*sizeof(abip_float));
        memset(D,0,m*sizeof(abip_float));
        count = 0;
        abip_float alpha_pc = 1;
        if(A != ABIP_NULL){         
            for (int i = 0; i < n; i++) {
                for (int j = A->p[i]; j < A->p[i + 1]; j++) {
                    E1[i] += POWF(ABS(A->x[j]), alpha_pc);
                }
                E1[i] = SQRTF(POWF(E1[i], 1/alpha_pc));
            }
        }
        if(Q != ABIP_NULL){

            for (int i = 0; i < n; i++) {
                for (int j = Q->p[i]; j < Q->p[i + 1]; j++) {
                    E2[i] += POWF(ABS(Q->x[j]), alpha_pc);
                }
                E2[i] = SQRTF(POWF(E2[i], 1/alpha_pc));
            }

        }

        for (int i = 0; i < n; i++){
            E[i] = E1[i] < E2[i] ? E2[i] : E1[i];
        }

        if(k->q){
            for(int i=0;i<k->qsize;i++){
                mean_E = ABIP(vec_mean)(&E[count],k->q[i]);
                for(int j=0;j<k->q[i];j++){
                    E[j+count] = mean_E;
                }
                count += k->q[i];
            }
        }

        if(k->rq){
            for(int i=0;i<k->rqsize;i++){
                mean_E = ABIP(vec_mean)(&E[count],k->rq[i]);
                for(int j=0;j<k->rq[i];j++){
                    E[j+count] = mean_E;
                }
                count += k->rq[i];
            }
        }

        if(A != ABIP_NULL){         
            for(int i=0;i<A->p[n];i++){
                D[A->i[i]] += POWF(ABS(A->x[i]), 2-alpha_pc);
            }
            for(int i=0;i<m;i++){
                D[i] = SQRTF(POWF(D[i], 1/(2-alpha_pc)));
                if (D[i] < min_row_scale) D[i] = 1;
                else if (D[i] > max_row_scale) D[i] = max_row_scale;
            }

            for(int i=0;i<n;i++){
                if (E[i] < min_col_scale) E[i] = 1;
                else if (E[i] > max_col_scale) E[i] = max_col_scale;
                for(int j=A->p[i];j<A->p[i+1];j++){
                    A->x[j] /= E[i];
                }
            }
        }

        if(Q != ABIP_NULL){

            for(int i=0;i<n;i++){
                for(int j=Q->p[i];j<Q->p[i+1];j++){
                    Q->x[j] /= E[i];
                }
            }
            for(int i=0;i<Q->p[n];i++){
                Q->x[i] /= E[Q->i[i]];
            }

        }
    
        if(A != ABIP_NULL){         
            for(int i=0;i<A->p[n];i++){
                A->x[i] /= D[A->i[i]];
            }
        }

        for(int i=0;i<n;i++){
            E_hat[i] *= E[i];
        }


        for(int i=0;i<m;i++){
            D_hat[i] *= D[i];
        }
    }


    abip_float sc = SQRTF(SQRTF(ABIP(norm_sq)(self->c, self->q) + ABIP(norm_sq)(self->b, m)));


    if(self->b != ABIP_NULL){
        for(int i=0;i<m;i++){
            self->b[i] /= D_hat[i];
        }
    }

    for(int i=0;i<n;i++){
        self->c[i] /= E_hat[i];
    }

    if(sc < MIN_SCALE) sc = 1;
    else if(sc > MAX_SCALE) sc = MAX_SCALE;
    self->sc_b = 1/sc;
    self->sc_c = 1/sc;

    if(self->b != ABIP_NULL){
        ABIP(scale_array)(self->b,self->sc_b * self->stgs->scale, m);
    }
    ABIP(scale_array)(self->c, self->sc_c * self->stgs->scale, self->q);
    

    for(int i=0; i<m; i++){
        self->F[i] = self->stgs->rho_y + 2/self->stgs->rho_x / POWF(self->D[i], 2);
    }

    for(int i=0; i<self->q; i++){
        if(i<self->n) self->H[i] = self->stgs->rho_x + self->Q->x[i];
        else self->H[i] = self->stgs->rho_x;
    }

    abip_free(E);
    abip_free(E1);
    abip_free(E2);
    abip_free(D);

}


void un_scaling_svmqp_sol(svmqp *self, ABIPSolution *sol){

    abip_int m = self->m;
    abip_int n = self->n;
    abip_float *x = sol->x;
    abip_float *y = sol->y;
    abip_float *s = sol->s;

    abip_float *w = (abip_float *)abip_malloc(n * sizeof(abip_float));
    abip_float *b = (abip_float *)abip_malloc(sizeof(abip_float));
    abip_float *xi = (abip_float *)abip_malloc(m * sizeof(abip_float));

    for (int i = 0; i < self->q; ++i)
    {
        sol->x[i] /= (self->E[i] * self->sc_b);
    }

    memcpy(w, x, n*sizeof(abip_float));
    b[0] = x[n];
    memcpy(xi, &x[n+1], m*sizeof(abip_float));

    abip_free(x);
    abip_free(y);
    abip_free(s);

    sol->x = w;
    sol->y = b;
    sol->s = xi;

}

void calc_svmqp_residuals(svmqp *self, ABIPWork *w, ABIPResiduals *r, abip_int ipm_iter, abip_int admm_iter)
{
    DEBUG_FUNC

    abip_int n = w->n;
    abip_int m = w->m;

    abip_float *y = (abip_float*)abip_malloc(m*sizeof(abip_float));
    abip_float *x = (abip_float*)abip_malloc(n*sizeof(abip_float));
    abip_float *s = (abip_float*)abip_malloc(n*sizeof(abip_float));

    abip_float this_pr;
    abip_float this_dr;
    abip_float this_gap;

    if (admm_iter && r->last_admm_iter == admm_iter)
    {
        RETURN;
    }

    r->last_ipm_iter = ipm_iter;
    r->last_admm_iter = admm_iter;

    r->tau = ABS(w->u[n + m]);
    r->kap = ABS(w->v_origin[n + m]) / (self->stgs->normalize ? (self->stgs->scale * self->sc_c * self->sc_b) : 1);

    memcpy(y, w->u, m*sizeof(abip_float));
    memcpy(x, &w->u[m], n*sizeof(abip_float));
    memcpy(s, &w->v_origin[m], n*sizeof(abip_float));

    ABIP(scale_array)(y, 1/r->tau, m);
    ABIP(scale_array)(x, 1/r->tau, n);
    ABIP(scale_array)(s, 1/r->tau, n);

    abip_float *Ax = (abip_float*)abip_malloc(m*sizeof(abip_float));
    abip_float *Ax_b = (abip_float*)abip_malloc(m*sizeof(abip_float));

    memset(Ax, 0, m*sizeof(abip_float));
    self->spe_A_times(self, x, Ax);

    memcpy(Ax_b, Ax, m*sizeof(abip_float));
    ABIP(add_scaled_array)(Ax_b, self->b, m, -1);

    r->Ax_b_norm = ABIP(norm_inf)(Ax_b, m);

    ABIP(c_dot)(Ax, self->D, m);
    ABIP(c_dot)(Ax_b, self->D, m);
    
    this_pr = ABIP(norm_inf)(Ax_b, m) / (self->sc_b + MAX(ABIP(norm_inf)(Ax, m), self->sc_b * w->nm_inf_b));
    
    abip_float *Qx = (abip_float*)abip_malloc(n*sizeof(abip_float));
    abip_float *ATy = (abip_float*)abip_malloc(n*sizeof(abip_float));
    abip_float *Qx_ATy_c_s = (abip_float*)abip_malloc(n*sizeof(abip_float));


    memset(Qx, 0, n*sizeof(abip_float));
    abip_float xQx_2 = 0;

    if(self->Q != ABIP_NULL){
        ABIP(accum_by_A)(self->Q, x, Qx);
        xQx_2 = ABIP(dot)(x, Qx, n) / (2 * self->sc_b *self->sc_c);
    }

    memset(ATy, 0, n*sizeof(abip_float));
    self->spe_AT_times(self, y, ATy);

    memcpy(Qx_ATy_c_s, Qx, n*sizeof(abip_float));
    ABIP(add_scaled_array)(Qx_ATy_c_s, ATy, n, -1);
    ABIP(add_scaled_array)(Qx_ATy_c_s, self->c, n, 1);
    ABIP(add_scaled_array)(Qx_ATy_c_s, s, n, -1);

    r->Qx_ATy_c_s_norm = ABIP(norm_inf)(Qx_ATy_c_s, n);

    ABIP(c_dot)(Qx, self->E, n);
    ABIP(c_dot)(ATy, self->E, n);
    ABIP(c_dot)(Qx_ATy_c_s, self->E, n);
    ABIP(c_dot)(s, self->E, n);

    this_dr = ABIP(norm_inf)(Qx_ATy_c_s, n) / (self->sc_c + MAX(self->sc_c * w->nm_inf_c, ABIP(norm_inf)(Qx, n)));

    abip_float cTx = ABIP(dot)(self->c, x, n) / (self->sc_b * self->sc_c);
    abip_float bTy = ABIP(dot)(self->b, y, m) / (self->sc_b * self->sc_c);

    this_gap = ABS(2 * xQx_2 + cTx - bTy) / (1 + MAX(2*xQx_2, MAX(ABS(cTx), ABS(bTy))));

    r->pobj = xQx_2 + cTx;
    r->dobj = -xQx_2 + bTy;

    r->res_dif = MAX(MAX(ABS(this_pr - r->res_pri),ABS(this_dr - r->res_dual)),ABS(this_gap - r->rel_gap));
    r->res_pri = this_pr;
    r->res_dual = this_dr;
    r->rel_gap = this_gap;
    r->error_ratio = MAX(r->res_pri/self->stgs->eps_p, MAX(r->res_dual/self->stgs->eps_d, r->rel_gap/self->stgs->eps_g));

    if(ABIP(dot)(self->c, &w->u[m], n) < 0){

        ABIP(scale_array)(Qx, r->tau, n);
        ABIP(scale_array)(Ax, r->tau, m);
        r->res_unbdd = MAX(ABIP(norm)(Qx, n), ABIP(norm)(Ax, m)) / (-ABIP(dot)(self->c, &w->u[m], n));
    }
    else{
        r->res_unbdd = INFINITY;
    }

    if(ABIP(dot)(self->b, w->u, m) > 0){

        ABIP(scale_array)(ATy, r->tau, n);
        ABIP(scale_array)(s, r->tau, n);
        ABIP(add_scaled_array)(ATy, s, n, 1);

        r->res_infeas = ABIP(norm)(ATy, n) / ABIP(dot)(self->b, w->u, m);
    }
    else{
        r->res_infeas = INFINITY;
    }

    abip_free(x);
    abip_free(y);
    abip_free(s);
    abip_free(Ax);
    abip_free(Ax_b);
    abip_free(Qx);
    abip_free(ATy);
    abip_free(Qx_ATy_c_s);

}


cs* form_svmqp_kkt
(
    svmqp *self
)
{
    
    abip_int n = self->n+1;
    abip_int m = self->m;
    cs *LTL;

    cs *B1 = cs_spalloc(m, n, self->A->p[n], 1, 0);
    memcpy(B1->i,self->A->i,self->A->p[n]*sizeof(abip_int));
    memcpy(B1->p,self->A->p,(n+1)*sizeof(abip_int));
    memcpy(B1->x,self->A->x,self->A->p[n]*sizeof(abip_float));

    cs *B2 = cs_spalloc(m, n, self->A->p[n], 1, 0);
    memcpy(B2->i, B1->i, B1->p[n]*sizeof(abip_int));
    memcpy(B2->p, B1->p, (n+1)*sizeof(abip_int));
    memcpy(B2->x, B1->x, B1->p[n]*sizeof(abip_float));


    if(m > n){

        cs *T1 = cs_spalloc (n, n, n, 1, 1) ;   


        for (int i = 0 ; i < n; i++){
            cs_entry (T1, i, i, self->H[i]) ;
        } 

        cs *eye = cs_compress (T1) ;             

        cs_spfree (T1);

        for(int i=0;i<B1->p[n];i++){
            B1->x[i] /= self->F[B1->i[i]];
        }

        LTL = cs_add(eye,cs_multiply(cs_transpose(B2,1), B1),1,1);

        cs_spfree(eye);
    }
    else{

        for(int i=0;i<B2->n;i++){
            for(int j=B2->p[i];j<B2->p[i+1];j++){
                B2->x[j] /= self->H[i];
            }
        }
        cs *diag = cs_spalloc (m, m, m, 1, 1) ;   

        for (int i = 0 ; i < m; i++){
            cs_entry (diag, i, i, self->F[i]);
        } 
        cs *diag_F = cs_compress(diag);
        cs_spfree(diag);

        LTL = cs_add(diag_F,cs_multiply(B2, cs_transpose(B1,1)),1,1);

    }
    cs_spfree(B2);
    cs_spfree(B1);

    for(int i=0;i<LTL->n;i++){
        for(int j=LTL->p[i];j<LTL->p[i+1];j++){
            if(LTL->i[j]>i) LTL->x[j] = 0;
        }
    }
    cs_dropzeros (LTL);
    return LTL;
}


void init_svmqp_precon(svmqp *self){

    abip_int i;

    self->L->M = (abip_float*)abip_malloc(self->p * sizeof(abip_float));
    memset(self->L->M, 0, self->p * sizeof(abip_float));

    for (i = 0; i < self->A->n; i++) {
        for (int j = self->A->p[i]; j < self->A->p[i + 1]; j++) {
            self->L->M[self->A->i[j]] += self->A->x[j] * self->A->x[j] / self->H[i];
        }
    }

    for(i=0; i<self->p; i++){
        self->L->M[i] += 1/(self->D[i] * self->D[i] * self->H[self->A->n + i]);
        self->L->M[i] += 1/(self->E[i] * self->E[i] * self->H[self->A->n + self->p + i]);
    }

    ABIP(add_scaled_array)(self->L->M, self->rho_dr, self->p, 1.0);

    for(i=0; i<self->p; i++){
        self->L->M[i] = 1.0 / self->L->M[i];

    }

}


abip_float get_svmqp_pcg_tol(abip_int k, abip_float error_ratio, abip_float norm_p){

    if(k == -1){
        return 1e-9 * norm_p; 
    } 
    else{
        return MAX(1e-9, 1e-5 * norm_p / POWF((k+1),2));       
    }
                                 
}


abip_int init_svmqp_linsys_work(svmqp *self){

    if(self->stgs->linsys_solver == 0){// mkl need lower triangle
        self->L->K = cs_transpose(form_svmqp_kkt(self),1);
    }
    else if(self->stgs->linsys_solver == 1){// qdldl need upper triangle
        self->L->K = form_svmqp_kkt(self);
    }
    else if(self->stgs->linsys_solver == 2){// cholesky need upper triangle
        self->L->K = form_svmqp_kkt(self);
    }
    else if(self->stgs->linsys_solver == 3){// pcg doesn't need kkt matrix
        init_svmqp_precon(self);
        self->L->K = ABIP_NULL;
    }
    else if(self->stgs->linsys_solver == 4){// mkl-pardiso need lower triangle
        self->L->K = cs_transpose(form_svmqp_kkt(self),1);
    }
    else if(self->stgs->linsys_solver == 5){//  dense cholesky need upper triangle
        self->L->K = form_svmqp_kkt(self);
    }
    else{
        printf("\nlinsys solver type error\n");
        return -1;
    }
    return ABIP(init_linsys_work)(self);
}


abip_int solve_svmqp_linsys(svmqp *self, abip_float *b, abip_float *pcg_warm_start, abip_int iter, abip_float error_ratio){

    ABIP(timer) linsys_timer;
    ABIP(tic)(&linsys_timer);
    abip_int n = self->n;
    abip_int m = self->m;


    abip_int p = self->p;
    abip_int q = self->q;

    if(self->stgs->linsys_solver == 3){  //pcg

        abip_int n = self->q;
        abip_int m = self->p;
        abip_int i;

        abip_float norm_p = ABIP(norm)(&b[m], n);

        abip_float *tem = (abip_float*)abip_malloc(sizeof(abip_float)*m);
        memcpy(tem, b, m*sizeof(abip_float));
        for(i=0; i<m; i++){
            tem[i] /= self->rho_dr[i];
        }
        self->spe_AT_times(self, tem, &b[m]);

        abip_free(tem);

        abip_float pcg_tol = get_svmqp_pcg_tol(iter, error_ratio, norm_p);
        abip_int cg_its = ABIP(solve_linsys)(self, &b[m], n, &pcg_warm_start[m], pcg_tol);

        if (iter >= 0)
        {
            self->L->total_cg_iters += cg_its;
        }

        ABIP(scale_array)(b, -1, m);
        self->spe_A_times(self, &b[m], b);
        for(i=0; i<m; i++){
            b[i] /= -self->rho_dr[i];
        }
    }    
    else{  // direct methods

        abip_float *b2 = (abip_float*)abip_malloc(p*sizeof(abip_float));
        memcpy(b2, b, p*sizeof(abip_float));
        abip_float *tmp = (abip_float*)abip_malloc(q*sizeof(abip_float));
        memcpy(tmp, &b[p], sizeof(abip_float)*q);

        for(int i=0; i<q; i++){
            tmp[i] /= -self->H[i];
        }

        self->spe_A_times(self, tmp, b2);

        if(m > n + 1 && self->stgs->linsys_solver != 3){

            for(int i=0;i<p;i++){
                b2[i] /= self->F[i];
            }

            abip_float *tmp1 = (abip_float*)abip_malloc((n+1)*sizeof(abip_float));
            memset(tmp1,0,(n+1)*sizeof(abip_float));

            
            ABIP(accum_by_Atrans)(self->A, b2, tmp1);


            ABIP(solve_linsys)(self, tmp1, n+1, ABIP_NULL, 0);
            
            abip_float *tmp2 = (abip_float*)abip_malloc(m*sizeof(abip_float));
            memset(tmp2,0,m*sizeof(abip_float));

            ABIP(accum_by_A)(self->A, tmp1, tmp2);

            for(int i=0;i<m;i++){
                tmp2[i] /= self->F[i];
            }

            ABIP(add_scaled_array)(b2, tmp2, m,-1);

            abip_free(tmp1);
            abip_free(tmp2);
        }
        else{
        
            abip_int cg_its = ABIP(solve_linsys)(self, b2, m, pcg_warm_start, error_ratio);
            if (iter >= 0)
            {
                self->L->total_cg_iters += cg_its;
            }
        }

        memcpy(b, b2, m*sizeof(abip_float));
        abip_free(b2);
        abip_free(tmp);

    } // if pcg


    self->spe_AT_times(self, b, &b[p]);

    for(int i=0; i<q; i++){
         b[p + i] /= self->H[i];
    }

    self->L->total_solve_time += ABIP(tocq)(&linsys_timer);

    return 0;
}


void free_svmqp_linsys_work(svmqp *self)
{
    ABIP(free_linsys)(self);
}