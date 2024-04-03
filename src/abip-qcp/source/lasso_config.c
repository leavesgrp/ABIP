#include "lasso_config.h"
#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)

/**
@brief Initialize the lasso problem structure
*/
abip_int init_lasso(lasso **self, ABIPData *d, ABIPSettings *stgs) {
  lasso *this_lasso = (lasso *)abip_malloc(sizeof(lasso));
  *self = this_lasso;

  this_lasso->m = d->m;
  this_lasso->n = d->n;
  this_lasso->p = d->m + 1;
  this_lasso->q = 2 + 2 * d->n + d->m;
  abip_int m = this_lasso->p;
  abip_int n = this_lasso->q;

  this_lasso->lambda = d->lambda;
  this_lasso->Q = ABIP_NULL;
  this_lasso->sparsity = (((abip_float)d->A->p[d->n] / (d->m * d->n)) < 0.1);

  // non-identity DR scaling
  this_lasso->rho_dr =
      (abip_float *)abip_malloc((m + n + 1) * sizeof(abip_float));
  for (int i = 0; i < m + n + 1; i++) {
    if (i < m) {
      this_lasso->rho_dr[i] = stgs->rho_y;
    } else if (i < m + n) {
      this_lasso->rho_dr[i] = stgs->rho_x;
    } else {
      this_lasso->rho_dr[i] = stgs->rho_tau;
    }
  }

  if (this_lasso->sparsity) {
    this_lasso->sc = 2;
    this_lasso->sc_c = 1 / d->lambda;
    this_lasso->sc_cone2 = d->lambda / d->m * 80;
    this_lasso->sc_cone1 = 0.8 / this_lasso->sc_c / this_lasso->sc_cone2;
    this_lasso->sc_b = this_lasso->sc_c * 300 * d->lambda / d->m;
  } else {
    if (d->m < d->n)
      this_lasso->sc = 4;
    else
      this_lasso->sc = 1;
    this_lasso->sc_c = 1 / d->lambda;
    this_lasso->sc_b = this_lasso->sc_c;
    this_lasso->sc_cone2 = 0.8;
    this_lasso->sc_cone1 = 1 / this_lasso->sc_c;
  }

  this_lasso->L = (ABIPLinSysWork *)abip_malloc(sizeof(ABIPLinSysWork));
  this_lasso->pro_type = LASSO;
  this_lasso->stgs = stgs;
  this_lasso->data = d;

  abip_float *data_b =
      (abip_float *)abip_malloc(this_lasso->p * sizeof(abip_float));
  data_b[0] = 1;
  memcpy(&data_b[1], d->b, d->m * sizeof(abip_float));
  this_lasso->data->b = data_b;

  abip_float *data_c =
      (abip_float *)abip_malloc(this_lasso->q * sizeof(abip_float));
  data_c[0] = 0;
  data_c[1] = 2;
  memset(&data_c[2], 0, d->m * sizeof(abip_float));
  for (int i = 0; i < 2 * d->n; i++) {
    data_c[i + 2 + d->m] = d->lambda;
  }
  this_lasso->data->c = data_c;
  /* for lasso problem, no need for inputing c*/
  d->c = data_c;

  this_lasso->A = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
  this_lasso->b = (abip_float *)abip_malloc(this_lasso->p * sizeof(abip_float));
  this_lasso->c = (abip_float *)abip_malloc(this_lasso->q * sizeof(abip_float));
  this_lasso->D = (abip_float *)abip_malloc(this_lasso->m * sizeof(abip_float));
  this_lasso->E = (abip_float *)abip_malloc(this_lasso->n * sizeof(abip_float));

  this_lasso->scaling_data = &scaling_lasso_data;
  this_lasso->un_scaling_sol = &un_scaling_lasso_sol;
  this_lasso->calc_residuals = &calc_lasso_residuals;
  this_lasso->init_spe_linsys_work = &init_lasso_linsys_work;
  this_lasso->solve_spe_linsys = &solve_lasso_linsys;
  this_lasso->free_spe_linsys_work = &free_lasso_linsys_work;
  this_lasso->spe_A_times = &lasso_A_times;
  this_lasso->spe_AT_times = &lasso_AT_times;
  this_lasso->inner_conv_check = &lasso_inner_conv_check;

  return 0;
}

/**
@brief Customized matrix-vector multiplication for the lasso problem with A
untransposed
*/
void lasso_A_times(lasso *self, const abip_float *x, abip_float *y) {
  y[0] += x[0];

  for (int i = 1; i < self->p; i++) {
    y[i] += self->D[i - 1] * SQRTF(self->sc_cone2) * x[i + 1];
  }

  ABIP(accum_by_A)(self->A, &x[self->m + 2], &y[1]);
  ABIP(scale_array)(&y[1], -1, self->m);
  ABIP(accum_by_A)(self->A, &x[self->m + self->n + 2], &y[1]);
  ABIP(scale_array)(&y[1], -1, self->m);
}

/**
@brief Customized matrix-vector multiplication for the lasso problem with A
transposed
*/
void lasso_AT_times(lasso *self, const abip_float *x, abip_float *y) {
  y[0] += x[0];
  for (int i = 2; i < self->m + 2; i++) {
    y[i] += x[i - 1] * self->D[i - 2] * SQRTF(self->sc_cone2);
  }

  ABIP(accum_by_Atrans)(self->A, &x[1], &y[self->m + 2]);
  ABIP(scale_array)(&y[self->m + 2 + self->n], -1, self->n);
  ABIP(accum_by_Atrans)(self->A, &x[1], &y[self->m + 2 + self->n]);
  ABIP(scale_array)(&y[self->m + 2 + self->n], -1, self->n);
}

/**
@brief Customized scaling procedure for the lasso problem
*/
void scaling_lasso_data(lasso *self, ABIPCone *k) {
  if (!ABIP(copy_A_matrix)(&(self->A), self->data->A)) {
    abip_printf("ERROR: copy A matrix failed\n");
    RETURN;
  }

  abip_int m = self->m;
  abip_int n = self->n;
  ABIPMatrix *A = self->A;

  abip_float min_row_scale = MIN_SCALE * SQRTF((abip_float)n);
  abip_float max_row_scale = MAX_SCALE * SQRTF((abip_float)n);
  abip_float min_col_scale = MIN_SCALE * SQRTF((abip_float)m);
  abip_float max_col_scale = MAX_SCALE * SQRTF((abip_float)m);

  abip_float *E = self->E;
  memset(E, 0, n * sizeof(abip_float));

  abip_float *D = self->D;
  memset(D, 0, m * sizeof(abip_float));

  abip_float avg = 0;
  abip_float avg1 = 0;

  if (self->stgs->scale_E) {
    if (self->sparsity) {
      for (int i = 0; i < n; i++) {
        for (int j = A->p[i]; j < A->p[i + 1]; j++) {
          E[i] += A->x[j] * A->x[j];
        }

        avg += SQRTF(E[i]);
      }

      avg /= n;

      for (int i = 0; i < n; i++) {
        E[i] = avg / SQRTF(E[i] + 1e-4) / self->sc;
        if (E[i] > 1000 * SQRTF(m)) {
          E[i] = 1000 * SQRTF(m);
        }
        if (E[i] < 0.001 * SQRTF(m)) {
          E[i] = 1;
        }
        if (E[i] > 50) {
          E[i] = 50;
        }

        avg1 += E[i];
      }

      avg1 /= n;

      for (int i = 0; i < n; i++) {
        E[i] = avg1 / E[i] / self->sc;
      }
    } else {
      for (int i = 0; i < n; i++) {
        for (int j = A->p[i]; j < A->p[i + 1]; j++) {
          E[i] += A->x[j] * A->x[j];
        }

        E[i] = SQRTF(E[i]);

        if (E[i] > 1000 * SQRTF(m)) {
          E[i] = 1000 * SQRTF(m);
        }
        if (E[i] < 0.001 * SQRTF(m)) {
          E[i] = 1;
        }
        if (E[i] > 7) {
          E[i] = 7;
        }

        E[i] = 1 / (E[i] * self->sc);
      }
    }

    for (int i = 0; i < n; i++) {
      for (int j = A->p[i]; j < A->p[i + 1]; j++) {
        A->x[j] *= E[i];
      }
    }
  }

  for (int i = 0; i < A->p[n]; i++) {
    D[A->i[i]] += A->x[i] * A->x[i];
  }

  avg = 0;
  for (int i = 0; i < m; i++) {
    avg += SQRTF(2 * D[i] + self->sc_cone2);
  }

  avg /= m;

  for (int i = 0; i < m; i++) {
    D[i] = avg / SQRTF(2 * D[i] + self->sc_cone2);
  }

  for (int i = 0; i < A->p[n]; i++) {
    A->x[i] *= D[A->i[i]];
  }

  memcpy(self->b, self->data->b, self->p * sizeof(abip_float));

  self->b[0] = self->sc_cone1;

  for (int i = 1; i < m + 1; i++) {
    self->b[i] *= D[i - 1];
  }
  ABIP(scale_array)(self->b, self->sc_b, self->p);

  self->c[0] = 0;
  self->c[1] = self->sc_cone1 * self->sc_cone2;
  memset(&self->c[2], 0, self->m * sizeof(abip_float));

  for (int i = 0; i < self->n; i++) {
    self->c[i + self->m + 2] = E[i] * self->lambda;
    self->c[i + self->m + 2 + self->n] = E[i] * self->lambda;
  }
  ABIP(scale_array)(self->c, self->sc_c, self->q);

  self->D_hat = (abip_float *)abip_malloc(m * sizeof(abip_float));

  for (int i = 0; i < m; i++) {
    self->D_hat[i] =
        self->D[i] * self->D[i] * self->sc_cone2 + self->stgs->rho_y;
  }
}

/**
@brief Get the unscaled solution x of the lasso qcp problem
*/
void get_unscaled_x(lasso *self, abip_float *x, abip_float *us_x) {
  memcpy(us_x, x, self->q * sizeof(abip_float));
  for (int i = 0; i < self->m; i++) {
    us_x[i + 2] /= (self->sc_b / SQRTF(self->sc_cone2));
  }
  for (int i = 0; i < self->n; i++) {
    us_x[self->m + 2 + i] /= (self->E[i] * self->sc_b);
    us_x[self->m + self->n + 2 + i] /= (self->E[i] * self->sc_b);
  }
}

/**
@brief Get the unscaled solution y of the lasso qcp problem
*/
void get_unscaled_y(lasso *self, abip_float *y, abip_float *us_y) {
  memcpy(us_y, y, self->p * sizeof(abip_float));
  for (int i = 0; i < self->m; i++) {
    us_y[i + 1] *= (self->D[i] / self->sc_c);
  }
}

/**
@brief Get the unscaled solution s of the lasso qcp problem
*/
void get_unscaled_s(lasso *self, abip_float *s, abip_float *us_s) {
  memcpy(us_s, s, self->q * sizeof(abip_float));
  for (int i = 0; i < self->m; i++) {
    us_s[i + 2] /= (self->sc_c / SQRTF(self->sc_cone2));
  }
  for (int i = 0; i < self->n; i++) {
    us_s[self->m + 2 + i] /= (self->E[i] * self->sc_c);
    us_s[self->m + self->n + 2 + i] /= (self->E[i] * self->sc_c);
  }
}

/**
@brief Get the unscaled solution of the original lasso problem
*/
void un_scaling_lasso_sol(lasso *self, ABIPSolution *sol) {
  abip_int m = self->m;
  abip_int n = self->n;

  abip_float *beta = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(beta, &sol->x[m + 2], n * sizeof(abip_float));
  ABIP(add_scaled_array)(beta, &sol->x[m + n + 2], n, -1);
  ABIP(c_dot)(beta, self->E, n);
  ABIP(scale_array)(beta, 1 / self->sc_b, n);

  abip_free(sol->x);
  abip_free(sol->y);
  abip_free(sol->s);

  sol->x = beta;
}

/**
@brief Check whether the inner loop of the lasso problem has converged
*/
abip_float lasso_inner_conv_check(lasso *self, ABIPWork *w) {
  abip_int m = self->p;
  abip_int n = self->q;

  abip_float *Qu = (abip_float *)abip_malloc((m + n + 1) * sizeof(abip_float));
  abip_float *Mu = (abip_float *)abip_malloc((m + n) * sizeof(abip_float));

  memset(Mu, 0, (m + n) * sizeof(abip_float));

  self->spe_A_times(self, &w->u[m], Mu);

  self->spe_AT_times(self, w->u, &Mu[m]);
  ABIP(scale_array)(&Mu[m], -1, n);

  if (self->Q != ABIP_NULL) {
    ABIP(accum_by_A)(self->Q, &w->u[m], &Mu[m]);
  }

  memcpy(Qu, Mu, (m + n) * sizeof(abip_float));

  ABIP(add_scaled_array)(Qu, self->b, m, -w->u[m + n]);
  ABIP(add_scaled_array)(&Qu[m], self->c, n, w->u[m + n]);

  Qu[m + n] = -ABIP(dot)(w->u, Mu, m + n) / w->u[m + n] +
              ABIP(dot)(w->u, self->b, m) - ABIP(dot)(&w->u[m], self->c, n);

  abip_float *tem = (abip_float *)abip_malloc((m + n + 1) * sizeof(abip_float));
  memcpy(tem, Qu, (m + n + 1) * sizeof(abip_float));
  ABIP(add_scaled_array)(tem, w->v_origin, m + n + 1, -1);

  abip_float error_inner =
      ABIP(norm)(tem, m + n + 1) /
      (1 + ABIP(norm)(w->u, m + n + 1) + ABIP(norm)(w->v_origin, m + n + 1));

  abip_free(Qu);
  abip_free(Mu);
  abip_free(tem);

  return error_inner;
}

/**
@brief Calculate the residuals of the lasso qcp problem
*/
void calc_lasso_residuals(lasso *self, ABIPWork *w, ABIPResiduals *r,
                          abip_int ipm_iter, abip_int admm_iter) {
  abip_int p = w->m;
  abip_int q = w->n;
  abip_int m = self->m;
  abip_int n = self->n;
  abip_float this_pr;
  abip_float this_dr;
  abip_float this_gap;

  r->tau = w->u[p + q];

  abip_float *x = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(x, &(w->u[m + 3]), m * sizeof(abip_float));
  ABIP(scale_array)(x, SQRTF(self->sc_cone2) / (r->tau * self->sc_b), m);

  abip_float *beta_plus = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(beta_plus, &(w->u[2 * m + 3]), n * sizeof(abip_float));
  for (int i = 0; i < n; i++) {
    beta_plus[i] *= self->E[i] / (r->tau * self->sc_b);
  }

  abip_float *beta_minus = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(beta_minus, &(w->u[2 * m + 3 + n]), n * sizeof(abip_float));
  for (int i = 0; i < n; i++) {
    beta_minus[i] *= self->E[i] / (r->tau * self->sc_b);
  }

  abip_float *z = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(z, &(w->u[1]), m * sizeof(abip_float));
  for (int i = 0; i < m; i++) {
    z[i] *= self->D[i] / (r->tau * self->sc_c);
  }

  abip_float *s1 = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(s1, &(w->v[2 * m + 3]), n * sizeof(abip_float));
  for (int i = 0; i < n; i++) {
    s1[i] /= self->E[i] * r->tau * self->sc_c;
  }

  abip_float *s2 = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(s2, &(w->v[2 * m + n + 3]), n * sizeof(abip_float));
  for (int i = 0; i < n; i++) {
    s2[i] /= self->E[i] * r->tau * self->sc_c;
  }

  abip_float *pr = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(pr, x, m * sizeof(abip_float));
  ABIP(accum_by_A)(self->data->A, beta_plus, pr);
  ABIP(scale_array)(pr, -1, m);
  ABIP(accum_by_A)(self->data->A, beta_minus, pr);
  ABIP(scale_array)(pr, -1, m);
  ABIP(add_scaled_array)(pr, &self->data->b[1], m, -1);
  this_pr = ABIP(norm)(pr, m) / MAX(ABIP(norm)(&self->data->b[1], m), 1);

  abip_float *dr1 = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(dr1, s1, n * sizeof(abip_float));
  ABIP(accum_by_Atrans)(self->data->A, z, dr1);
  for (int i = 0; i < n; i++) {
    dr1[i] -= self->lambda;
  }

  abip_float *dr2 = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(dr2, s2, n * sizeof(abip_float));
  ABIP(scale_array)(dr2, -1, n);
  ABIP(accum_by_Atrans)(self->data->A, z, dr2);
  ABIP(scale_array)(dr2, -1, n);
  for (int i = 0; i < n; i++) {
    dr2[i] -= self->lambda;
  }
  abip_float *dr = (abip_float *)abip_malloc(2 * n * sizeof(abip_float));
  memcpy(dr, dr1, n * sizeof(abip_float));
  memcpy(&dr[n], dr2, n * sizeof(abip_float));
  this_dr = ABIP(norm)(dr, 2 * n) / (SQRTF(2 * n) * self->lambda);

  abip_float *lambda_ones = (abip_float *)abip_malloc(n * sizeof(abip_float));
  for (int i = 0; i < n; i++) {
    lambda_ones[i] = self->lambda;
  }
  this_gap =
      ABS(0.5 * ABIP(dot)(x, x, m) + ABIP(dot)(lambda_ones, beta_plus, n) +
          ABIP(dot)(lambda_ones, beta_minus, n) + 0.5 * ABIP(dot)(z, z, m) -
          ABIP(dot)(&self->data->b[1], z, m)) /
      (1 + ABS(0.5 * ABIP(dot)(x, x, m) + ABIP(dot)(lambda_ones, beta_plus, n) +
               ABIP(dot)(lambda_ones, beta_minus, n)));

  r->last_ipm_iter = ipm_iter;
  r->last_admm_iter = admm_iter;

  r->dobj = (-0.5 * ABIP(dot)(z, z, m) + ABIP(dot)(&self->data->b[1], z, m));
  r->pobj = (0.5 * ABIP(dot)(x, x, m) + ABIP(dot)(lambda_ones, beta_plus, n) +
             ABIP(dot)(lambda_ones, beta_minus, n));

  r->res_dif = MAX(MAX(ABS(this_pr - r->res_pri), ABS(this_dr - r->res_dual)),
                   ABS(this_gap - r->rel_gap));
  r->res_pri = this_pr;
  r->res_dual = this_dr;
  r->rel_gap = this_gap;
  r->error_ratio =
      MAX(r->res_pri / self->stgs->eps_p,
          MAX(r->res_dual / self->stgs->eps_d, r->rel_gap / self->stgs->eps_g));

  if (ABIP(dot)(self->c, &w->u[p], q) < 0) {
    abip_float *Ax = (abip_float *)abip_malloc(p * sizeof(abip_float));
    memset(Ax, 0, p * sizeof(abip_float));
    self->spe_A_times(self, &w->u[p], Ax);
    r->res_unbdd = ABIP(norm)(Ax, p) / (-ABIP(dot)(self->c, &w->u[p], q));
    abip_free(Ax);
  } else {
    r->res_unbdd = INFINITY;
  }

  if (ABIP(dot)(self->b, w->u, p) > 0) {
    abip_float *ATy_s = (abip_float *)abip_malloc(q * sizeof(abip_float));
    memset(ATy_s, 0, q * sizeof(abip_float));
    self->spe_AT_times(self, w->u, ATy_s);
    ABIP(add_scaled_array)(ATy_s, &w->v_origin[p], q, 1);

    r->res_infeas = ABIP(norm)(ATy_s, q) / ABIP(dot)(self->b, w->u, p);
    abip_free(ATy_s);
  } else {
    r->res_infeas = INFINITY;
  }

  abip_free(x);
  abip_free(beta_plus);
  abip_free(beta_minus);
  abip_free(z);
  abip_free(s1);
  abip_free(s2);
  abip_free(pr);
  abip_free(dr);
  abip_free(dr1);
  abip_free(dr2);
  abip_free(lambda_ones);
}

/**
@brief Formulate the qcp KKT matrix of the lasso problem
*/
cs *form_lasso_kkt(lasso *self) {
  abip_int n = self->n;
  abip_int m = self->m;
  cs *LTL;

  cs *Y1 = cs_spalloc(m, n, self->A->p[n], 1, 0);
  memcpy(Y1->i, self->A->i, self->A->p[n] * sizeof(abip_int));
  memcpy(Y1->p, self->A->p, (n + 1) * sizeof(abip_int));
  memcpy(Y1->x, self->A->x, self->A->p[n] * sizeof(abip_float));

  if (m > n) {
    cs *Y2 = cs_spalloc(m, n, self->A->p[n], 1, 0);
    memcpy(Y2->i, self->A->i, self->A->p[n] * sizeof(abip_int));
    memcpy(Y2->p, self->A->p, (n + 1) * sizeof(abip_int));
    memcpy(Y2->x, self->A->x, self->A->p[n] * sizeof(abip_float));

    cs *T1 = cs_spalloc(n, n, n, 1, 1); /* create triplet identity matrix */

    for (int i = 0; i < n; i++) {
      cs_entry(T1, i, i, 0.5);
    }

    cs *half_eye = cs_compress(T1);

    cs_spfree(T1);

    for (int i = 0; i < Y1->p[n]; i++) {
      Y1->x[i] /= self->D_hat[Y1->i[i]];
    }

    LTL = cs_add(half_eye, cs_multiply(cs_transpose(Y2, 1), Y1), 1, 1);

    cs_spfree(Y1);
    cs_spfree(Y2);
    cs_spfree(half_eye);
  } else {
    cs *YYT = cs_multiply(Y1, cs_transpose(Y1, 1));
    cs_spfree(Y1);

    cs *diag = cs_spalloc(m, m, m, 1, 1);

    for (int i = 0; i < m; i++) {
      cs_entry(diag, i, i, self->D_hat[i]);
    }
    cs *diag_D_hat = cs_compress(diag);
    cs_spfree(diag);

    LTL = cs_add(diag_D_hat, YYT, 1, 2);
    cs_spfree(YYT);
  }

  for (int i = 0; i < LTL->n; i++) {
    for (int j = LTL->p[i]; j < LTL->p[i + 1]; j++) {
      if (LTL->i[j] > i) LTL->x[j] = 0;
    }
  }
  cs_dropzeros(LTL);
  return LTL;
}

/**
@brief Initialize the preconditioner of conjugate gradient method for the lasso
problem
*/
void init_lasso_precon(lasso *self) {
  self->L->M = (abip_float *)abip_malloc(self->p * sizeof(abip_float));
  memset(self->L->M, 0, self->p * sizeof(abip_float));

  abip_float *M = self->L->M;

  M[0] = 1 / (self->stgs->rho_y + 1);

  for (int i = 0; i < self->A->p[self->A->n]; i++) {
    M[self->A->i[i] + 1] += 2 * POWF(self->A->x[i], 2);
  }

  for (int i = 1; i < self->p; i++) {
    M[i] = 1 / (self->stgs->rho_y + M[i] +
                self->sc_cone2 * POWF(self->D[i - 1], 2));
  }
}

/**
@brief Get the tolerance of the conjugate gradient method for the lasso problem
*/
abip_float get_lasso_pcg_tol(abip_int k, abip_float error_ratio,
                             abip_float norm_p) {
  if (k == -1) {
    return 1e-9 * norm_p;
  } else {
    if (error_ratio > 100000) {
      return MAX(1e-9, 1.2e-2 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 30000) {
      return MAX(1e-9, 8e-3 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 10000) {
      return MAX(1e-9, 6e-3 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 3000) {
      return MAX(1e-9, 5e-3 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 1000) {
      return MAX(1e-9, 3e-3 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 300) {
      return MAX(1e-9, 2e-3 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 100) {
      return MAX(1e-9, 1.5e-3 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 30) {
      return MAX(1e-9, 8e-4 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 10) {
      return MAX(1e-9, 6e-4 * norm_p / POWF((k + 1), 2));
    } else {
      return MAX(1e-9, 5e-4 * norm_p / POWF((k + 1), 2));
    }
  }
}

/**
@brief Initialize the linear system solver work space for the lasso problem
*/
abip_int init_lasso_linsys_work(lasso *self) {
  if (self->stgs->linsys_solver == 0) {
    self->L->K = cs_transpose(form_lasso_kkt(self), 1);
  } else if (self->stgs->linsys_solver == 1) {
    self->L->K = form_lasso_kkt(self);
  } else if (self->stgs->linsys_solver == 2) {
    self->L->K = form_lasso_kkt(self);
  } else if (self->stgs->linsys_solver == 3) {
    init_lasso_precon(self);
    self->L->K = ABIP_NULL;
  } else if (self->stgs->linsys_solver == 4) {
    self->L->K = cs_transpose(form_lasso_kkt(self), 1);
  } else if (self->stgs->linsys_solver == 5) {
    self->L->K = form_lasso_kkt(self);
  } else {
    printf("\nlinsys solver type error\n");
    return -1;
  }
  return ABIP(init_linsys_work)(self);
}

/**
@brief Customized linear system solver for the lasso problem
*/
abip_int solve_lasso_linsys(lasso *self, abip_float *b,
                            abip_float *pcg_warm_start, abip_int iter,
                            abip_float error_ratio) {
  ABIP(timer) linsys_timer;
  ABIP(tic)(&linsys_timer);

  ABIP(scale_array)(&b[self->p], -1, self->q);

  if (self->stgs->linsys_solver == 3) {  // pcg

    abip_int p = self->p;
    abip_float norm_p = ABIP(norm)(b, p);

    self->spe_A_times(self, &b[p], b);
    abip_float pcg_tol = get_lasso_pcg_tol(iter, error_ratio, norm_p);
    abip_int cg_its = ABIP(solve_linsys)(self, b, p, pcg_warm_start, pcg_tol);

    if (iter >= 0) {
      self->L->total_cg_iters += cg_its;
    }
  }

  else {  // direct methods
    abip_float *b2 = (abip_float *)abip_malloc(self->p * sizeof(abip_float));
    memcpy(b2, b, self->p * sizeof(abip_float));
    self->spe_A_times(self, &b[self->p], b2);
    b[0] = b2[0] / (1 + self->stgs->rho_y);
    abip_int n = self->n;
    abip_int m = self->m;

    if (m > n) {
      for (int i = 0; i < m; i++) {
        b2[i + 1] /= self->D_hat[i];
        b[i + 1] = b2[i + 1];
      }

      abip_float *tmp = (abip_float *)abip_malloc(n * sizeof(abip_float));
      memset(tmp, 0, n * sizeof(abip_float));
      ABIP(accum_by_Atrans)(self->A, &b2[1], tmp);

      ABIP(solve_linsys)(self, tmp, n, ABIP_NULL, 0);

      abip_float *tmp2 = (abip_float *)abip_malloc(m * sizeof(abip_float));
      memset(tmp2, 0, m * sizeof(abip_float));

      ABIP(accum_by_A)(self->A, tmp, tmp2);
      for (int i = 0; i < m; i++) {
        tmp2[i] /= self->D_hat[i];
      }

      ABIP(add_scaled_array)(&b[1], tmp2, m, -1);

      abip_free(tmp);
      abip_free(tmp2);
    } else {
      memcpy(&b[1], &b2[1], m * sizeof(abip_float));

      ABIP(solve_linsys)(self, &b[1], m, ABIP_NULL, 0);
    }

    abip_free(b2);
  }

  ABIP(scale_array)(&b[self->p], -1, self->q);
  self->spe_AT_times(self, b, &b[self->p]);

  self->L->total_solve_time += ABIP(tocq)(&linsys_timer);

  return 0;
}

/**
@brief Free the linear system solver work space for the lasso problem
*/
void free_lasso_linsys_work(lasso *self) { ABIP(free_linsys)(self); }