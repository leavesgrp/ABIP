#include "svm_config.h"
#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)

/**
@brief Initialize the svm socp formulation structure
*/
abip_int init_svm(svm **self, ABIPData *d, ABIPSettings *stgs) {
  svm *this_svm = (svm *)abip_malloc(sizeof(svm));
  *self = this_svm;

  this_svm->m = d->m;
  this_svm->n = d->n;
  this_svm->p = d->m + d->n + 1;
  this_svm->q = 4 + 3 * d->n + 2 * d->m;
  this_svm->lambda = d->lambda;
  abip_int m = this_svm->p;
  abip_int n = this_svm->q;
  this_svm->Q = ABIP_NULL;
  this_svm->sparsity = (((abip_float)d->A->p[d->n] / (d->m * d->n)) < 0.05);

  // non-identity DR scaling
  this_svm->rho_dr =
      (abip_float *)abip_malloc((m + n + 1) * sizeof(abip_float));
  for (int i = 0; i < m + n + 1; i++) {
    if (i < m) {
      this_svm->rho_dr[i] = stgs->rho_y;
    } else if (i < m + n) {
      this_svm->rho_dr[i] = stgs->rho_x;
    } else {
      this_svm->rho_dr[i] = stgs->rho_tau;
    }
  }

  this_svm->L = (ABIPLinSysWork *)abip_malloc(sizeof(ABIPLinSysWork));
  this_svm->pro_type = SVM;
  this_svm->stgs = stgs;

  this_svm->data = (ABIPData *)abip_malloc(sizeof(ABIPData));

  abip_float *data_b =
      (abip_float *)abip_malloc(this_svm->p * sizeof(abip_float));
  memset(data_b, 0, this_svm->p * sizeof(abip_float));
  for (int i = 0; i < d->m + 1; i++) {
    data_b[i] = 1;
  }
  this_svm->data->b = data_b;

  abip_float *data_c =
      (abip_float *)abip_malloc(this_svm->q * sizeof(abip_float));
  memset(data_c, 0, this_svm->q * sizeof(abip_float));

  data_c[1] = 1;

  for (int i = 0; i < d->m; i++) {
    data_c[i + 4 + 3 * d->n] = 1;
  }

  this_svm->data->c = data_c;

  /* for svm problem, no need for inputing c*/
  d->c = (abip_float *)abip_malloc(this_svm->q * sizeof(abip_float));
  memcpy(d->c, data_c, this_svm->q * sizeof(abip_float));

  if ((d->m < 10 * d->n) && (10 * d->m > d->n)) {
    this_svm->sc = 1;
    this_svm->sc_c = MAX(0.45, POWF(7.5, (-log(2 * d->lambda) / log(10))) * 2);
    this_svm->sc_b = 1;
    this_svm->sc_cone1 = MAX(3, log(2 * d->lambda) / log(10) * 4 + 4);
    this_svm->sc_cone2 = this_svm->sc_cone1;
  } else if (10 * d->m < d->n) {
    this_svm->sc = 1;
    this_svm->sc_b = 1;
    this_svm->sc_cone2 = MAX(3, log(2 * d->lambda) / log(10) * 2 + 2);
    if (d->lambda >= 1) {
      this_svm->sc_c =
          MAX(0.2, POWF(0.2, (log(2 * d->lambda) / log(10))) * 7.5);
      this_svm->sc_cone1 = this_svm->sc_cone2;
    } else {
      this_svm->sc_c = POWF(0.3, log(2 * d->lambda) / log(10)) * 3;
      this_svm->sc_cone1 = MAX(0.4, log(2 * d->lambda) / log(10) * 0.2 + 0.8);
    }
  } else if (d->m > 10 * d->n) {
    if (d->n < 10) {
      this_svm->sc = 1;
      this_svm->sc_c = 1 / this_svm->lambda;
      // this_svm->sc_c = 1.0/d->n;
      this_svm->sc_b = 1;
      this_svm->sc_cone1 = 6;
      if (d->lambda < 0.002) {
        this_svm->sc_cone2 =
            this_svm->sc_cone2 - 3 * log(this_svm->lambda * 500) / log(10);
      }
    } else if (d->lambda >= 1) {
      this_svm->sc = 1;
      this_svm->sc_c = 1 / this_svm->lambda;
      this_svm->sc_b = 1;
      this_svm->sc_cone1 = 6;
      this_svm->sc_cone2 = this_svm->lambda;
    } else {
      this_svm->sc = 1;
      this_svm->sc_c = MIN(POWF(5, (-log(5 * d->lambda) / log(10))) * 4, 300);
      this_svm->sc_b = MAX(0.1, log(5 * d->lambda) / log(10) * 0.2 + 0.9);
      this_svm->sc_cone1 = MAX(0.05, log(5 * d->lambda) / log(10) * 0.3 + 0.7);
      this_svm->sc_cone2 = -log(5 * d->lambda) / log(10) * 2 + 6;
      if (d->lambda < 0.002) {
        this_svm->sc_cone2 =
            this_svm->sc_cone2 - 3 * log(this_svm->lambda * 500) / log(10);
      }
    }
  }

  ABIPMatrix *data_A = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
  data_A->m = d->m;
  data_A->n = d->n + 1;
  data_A->p = (abip_int *)abip_malloc((data_A->n + 1) * sizeof(abip_int));
  data_A->i =
      (abip_int *)abip_malloc((d->A->p[d->n] + d->m) * sizeof(abip_int));
  data_A->x =
      (abip_float *)abip_malloc((d->A->p[d->n] + d->m) * sizeof(abip_float));
  for (int i = 0; i < d->A->p[d->n]; i++) {
    d->A->x[i] *= d->b[d->A->i[i]];
  }
  memcpy(data_A->p, d->A->p, (d->n + 1) * sizeof(abip_int));
  data_A->p[data_A->n] = d->A->p[d->n] + d->m;

  memcpy(data_A->i, d->A->i, d->A->p[d->n] * sizeof(abip_int));
  for (int i = d->A->p[d->n]; i < d->A->p[d->n] + d->m; i++) {
    data_A->i[i] = i - d->A->p[d->n];
  }

  memcpy(data_A->x, d->A->x, d->A->p[d->n] * sizeof(abip_float));
  for (int i = d->A->p[d->n]; i < d->A->p[d->n] + d->m; i++) {
    data_A->x[i] = d->b[i - d->A->p[d->n]];
  }
  this_svm->data->A = data_A;
  this_svm->A = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
  this_svm->b = (abip_float *)abip_malloc(this_svm->p * sizeof(abip_float));
  this_svm->c = (abip_float *)abip_malloc(this_svm->q * sizeof(abip_float));
  this_svm->sc_D = (abip_float *)abip_malloc(this_svm->m * sizeof(abip_float));
  this_svm->sc_E =
      (abip_float *)abip_malloc((this_svm->n + 1) * sizeof(abip_float));
  this_svm->sc_F =
      (abip_float *)abip_malloc((this_svm->n) * sizeof(abip_float));

  this_svm->wA = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
  this_svm->wy = (abip_float *)abip_malloc(this_svm->m * sizeof(abip_float));
  this_svm->wB = (abip_float *)abip_malloc(this_svm->m * sizeof(abip_float));
  this_svm->wC = (abip_float *)abip_malloc(this_svm->m * sizeof(abip_float));
  this_svm->wD = (abip_float *)abip_malloc((this_svm->n) * sizeof(abip_float));
  this_svm->wE = (abip_float *)abip_malloc((this_svm->n) * sizeof(abip_float));
  this_svm->wF = (abip_float *)abip_malloc(this_svm->m * sizeof(abip_float));
  this_svm->wG = (abip_float *)abip_malloc((this_svm->n) * sizeof(abip_float));
  this_svm->wH =
      (abip_float *)abip_malloc((this_svm->n + 1) * sizeof(abip_float));
  this_svm->wX = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));

  this_svm->scaling_data = &scaling_svm_data;
  this_svm->un_scaling_sol = &un_scaling_svm_sol;
  this_svm->calc_residuals = &calc_svm_residuals;
  this_svm->init_spe_linsys_work = &init_svm_linsys_work;
  this_svm->solve_spe_linsys = &solve_svm_linsys;
  this_svm->free_spe_linsys_work = &free_svm_linsys_work;
  this_svm->spe_A_times = &svm_A_times;
  this_svm->spe_AT_times = &svm_AT_times;
  this_svm->inner_conv_check = &svm_inner_conv_check;

  return 0;
}

/**
@brief Customized matrix-vector multiplication for the svm socp formulation with
A untransposed
*/
void svm_A_times(svm *self, const abip_float *x, abip_float *y) {
  y[0] += x[0];
  abip_int m = self->m;
  abip_int n = self->n;

  abip_float *tmp = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(tmp, &x[n + 2], n * sizeof(abip_float));
  ABIP(add_scaled_array)(tmp, &x[2 * n + 3], n, -1);
  ABIP(accum_by_A)(self->wA, tmp, &y[1]);

  ABIP(add_scaled_array)(&y[1], self->wy, m, (x[2 * n + 2] - x[3 * n + 3]));
  for (int i = 0; i < m; i++) {
    y[i + 1] +=
        (self->wB[i] * x[i + 3 * n + 4] - self->wC[i] * x[i + 3 * n + 4 + m]);
  }

  for (int i = 0; i < n; i++) {
    y[i + 1 + m] += (self->wD[i] * x[i + 2] - self->wE[i] * tmp[i]);
  }

  abip_free(tmp);
}

/**
@brief Customized matrix-vector multiplication for the svm socp formulation with
A transposed
*/
void svm_AT_times(svm *self, const abip_float *x, abip_float *y) {
  y[0] += x[0];
  abip_int m = self->m;
  abip_int n = self->n;

  for (int i = 0; i < n; i++) {
    y[i + 2] += self->wD[i] * x[i + m + 1];
  }

  abip_float *tmp = (abip_float *)abip_malloc(n * sizeof(abip_float));
  for (int i = 0; i < n; i++) {
    tmp[i] = -self->wE[i] * x[i + m + 1];
  }
  ABIP(accum_by_Atrans)(self->wA, &x[1], tmp);

  ABIP(add_scaled_array)(&y[n + 2], tmp, n, 1);
  ABIP(add_scaled_array)(&y[2 * n + 3], tmp, n, -1);

  y[2 * n + 2] += ABIP(dot)(self->wy, &x[1], m);
  y[3 * n + 3] -= ABIP(dot)(self->wy, &x[1], m);

  for (int i = 0; i < m; i++) {
    y[i + 3 * n + 4] += self->wB[i] * x[i + 1];
    y[i + 3 * n + 4 + m] -= self->wC[i] * x[i + 1];
  }

  abip_free(tmp);
}

/**
@brief Check whether the inner loop of the svm socp formulation has converged
*/
abip_float svm_inner_conv_check(svm *self, ABIPWork *w) {
  DEBUG_FUNC

  abip_int m = w->m;
  abip_int n = w->n;
  abip_int l = m + n + 1;

  abip_float *y = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(y, w->u, m * sizeof(abip_float));
  abip_float *x = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(x, &w->u[m], n * sizeof(abip_float));
  abip_float *s = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(s, &w->v_origin[m], n * sizeof(abip_float));

  abip_float tau = w->u[l - 1];
  abip_float kap = w->v_origin[l - 1];

  abip_float *row1 = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(row1, self->b, m * sizeof(abip_float));
  ABIP(scale_array)(row1, -tau, m);
  self->spe_A_times(self, x, row1);

  abip_float *row2 = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(row2, s, n * sizeof(abip_float));
  ABIP(add_scaled_array)(row2, self->c, n, -tau);
  self->spe_AT_times(self, y, row2);
  ABIP(scale_array)(row2, -1, n);

  abip_float *Pu_v = (abip_float *)abip_malloc(l * sizeof(abip_float));
  memcpy(Pu_v, row1, m * sizeof(abip_float));
  memcpy(&Pu_v[m], row2, n * sizeof(abip_float));
  Pu_v[l - 1] = ABIP(dot)(self->b, y, m) - ABIP(dot)(self->c, x, n) - kap;

  abip_float err_inner =
      ABIP(norm)(Pu_v, l) /
      (1 + SQRTF(ABIP(norm_sq)(w->u, l) + ABIP(norm_sq)(w->v_origin, l)));

  abip_free(x);
  abip_free(y);
  abip_free(s);
  abip_free(row1);
  abip_free(row2);
  abip_free(Pu_v);
  return err_inner;
}

/**
@brief Customized scaling procedure for the svm socp formulation
*/
void scaling_svm_data(svm *self, ABIPCone *k) {
  if (!ABIP(copy_A_matrix)(&(self->A), self->data->A)) {
    abip_printf("ERROR: copy A matrix failed\n");
    RETURN;
  }

  abip_int m = self->m;
  abip_int n = self->n + 1;
  ABIPMatrix *A = self->A;

  abip_float min_row_scale = MIN_SCALE * SQRTF((abip_float)n);
  abip_float max_row_scale = MAX_SCALE * SQRTF((abip_float)n);
  abip_float min_col_scale = MIN_SCALE * SQRTF((abip_float)m);
  abip_float max_col_scale = MAX_SCALE * SQRTF((abip_float)m);

  abip_float *E = self->sc_E;
  memset(E, 0, n * sizeof(abip_float));

  abip_float *D = self->sc_D;
  memset(D, 0, m * sizeof(abip_float));

  abip_float avg = 0;

  if (self->stgs->scale_E) {
    for (int i = 0; i < n; i++) {
      for (int j = A->p[i]; j < A->p[i + 1]; j++) {
        E[i] += A->x[j] * A->x[j];
      }
      E[i] = SQRTF(E[i]);
      avg += E[i];
    }

    avg /= n;

    for (int i = 0; i < n; i++) {
      E[i] = avg / E[i];
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
    avg += SQRTF(D[i]);
  }

  avg /= m;

  for (int i = 0; i < m; i++) {
    D[i] = avg / SQRTF(D[i]);
  }

  for (int i = 0; i < A->p[n]; i++) {
    A->x[i] *= D[A->i[i]];
  }

  for (int i = 0; i < n - 1; i++) {
    self->sc_F[i] = 1 / SQRTF(1 + 2 * E[i] * E[i]);
  }

  memcpy(self->b, self->data->b, self->p * sizeof(abip_float));

  self->b[0] = self->sc_cone2;

  for (int i = 1; i < m + 1; i++) {
    self->b[i] = D[i - 1];
  }
  ABIP(scale_array)(self->b, self->sc_b, self->p);

  memcpy(self->c, self->data->c, self->q * sizeof(abip_float));
  self->c[0] = 0;
  self->c[1] = self->sc_c * self->sc_cone1 * self->sc_cone2;

  for (int i = 0; i < m; i++) {
    self->c[i + 4 + 3 * self->n] = self->lambda * self->sc_c / self->sc;
  }

  n -= 1;

  ABIP(copy_A_matrix)(&(self->wA), self->A);
  self->wA->n -= 1;

  memcpy(self->wy, &self->A->x[self->A->p[n]], m * sizeof(abip_float));

  memcpy(self->wB, self->sc_D, m * sizeof(abip_float));
  ABIP(scale_array)(self->wB, 1 / self->sc, m);

  memcpy(self->wC, self->sc_D, m * sizeof(abip_float));

  memcpy(self->wD, self->sc_F, n * sizeof(abip_float));
  ABIP(scale_array)(self->wD, -SQRTF(self->sc_cone1), n);

  for (int i = 0; i < n; i++) {
    self->wE[i] = self->sc_E[i] * self->sc_F[i];
  }

  for (int i = 0; i < m; i++) {
    self->wF[i] = self->wB[i] * self->wB[i] + self->wC[i] * self->wC[i] +
                  self->stgs->rho_y;
  }

  for (int i = 0; i < n; i++) {
    self->wG[i] = self->wD[i] * self->wD[i] + 2 * self->wE[i] * self->wE[i] +
                  self->stgs->rho_y;
  }

  for (int i = 0; i < n; i++) {
    self->wH[i] = 2 - 4 / self->wG[i] * self->wE[i] * self->wE[i];
  }
  self->wH[n] = 2;

  ABIP(copy_A_matrix)(&(self->wX), self->wA);
  for (int i = 0; i < n; i++) {
    for (int j = self->wX->p[i]; j < self->wX->p[i + 1]; j++) {
      self->wX->x[j] *= (-2 * self->wE[i]);
    }
  }
}

/**
@brief Get the unscaled solution of the original svm problem
*/
void un_scaling_svm_sol(svm *self, ABIPSolution *sol) {
  abip_int m = self->m;
  abip_int n = self->n;
  abip_float *x = sol->x;
  abip_float *y = sol->y;
  abip_float *s = sol->s;

  abip_float *w = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memcpy(w, &x[n + 2], n * sizeof(abip_float));
  ABIP(add_scaled_array)(w, &x[2 * n + 3], n, -1);
  ABIP(c_dot)(w, self->sc_E, n);
  ABIP(scale_array)(w, 1 / self->sc_b, n);

  abip_float *b = (abip_float *)abip_malloc(sizeof(abip_float));
  b[0] = (x[2 * n + 2] - x[3 * n + 3]) * self->sc_E[n] / self->sc_b;

  abip_float *xi = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(xi, &x[3 * n + 4], m * sizeof(abip_float));
  ABIP(scale_array)(xi, 1 / (self->sc_b * self->sc_c), m);

  abip_free(x);
  abip_free(y);
  abip_free(s);

  sol->x = w;
  sol->y = b;
  sol->s = xi;
}

/**
@brief Calculate the residuals of the svm socp formulation
*/
void calc_svm_residuals(svm *self, ABIPWork *w, ABIPResiduals *r,
                        abip_int ipm_iter, abip_int admm_iter) {
  abip_int p = w->m;
  abip_int q = w->n;
  abip_int m = self->m;
  abip_int n = self->n;
  abip_float C = self->lambda;
  abip_float this_pr;
  abip_float this_dr;
  abip_float this_gap;

  r->tau = w->u[p + q];

  abip_float *w1 = (abip_float *)abip_malloc((n + 1) * sizeof(abip_float));
  memcpy(w1, &w->u[m + 2 * n + 3], n * sizeof(abip_float));
  ABIP(add_scaled_array)(w1, &w->u[m + 3 * n + 4], n, -1);
  for (int i = 0; i < n; i++) {
    w1[i] *= (self->sc_E[i] / (r->tau * self->sc_b));
  }

  abip_float b = (w->u[m + 3 * n + 3] - w->u[m + 4 * n + 4]) / r->tau *
                 self->sc_E[n] / self->sc_b;
  w1[n] = b;

  abip_float *xi = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(xi, &w->u[m + 4 * n + 5], m * sizeof(abip_float));
  ABIP(scale_array)(xi, 1 / (r->tau * self->sc * self->sc_b), m);

  abip_float *t = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(t, &w->u[2 * m + 4 * n + 5], m * sizeof(abip_float));
  ABIP(scale_array)(t, 1 / (r->tau * self->sc_b), m);

  abip_float *y = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(y, &w->u[1], m * sizeof(abip_float));
  for (int i = 0; i < m; i++) {
    y[i] = y[i] / r->tau * self->sc_D[i] / self->sc_c;
  }

  abip_float *s2 = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(s2, &w->v[2 * m + 4 * n + 5], m * sizeof(abip_float));
  ABIP(scale_array)(s2, 1 / (r->tau * self->sc_c), m);

  abip_float *s1 = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(s1, &w->v[m + 4 * n + 5], m * sizeof(abip_float));
  ABIP(scale_array)(s1, self->sc / (r->tau * self->sc_c), m);

  abip_float *pr = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memcpy(pr, xi, m * sizeof(abip_float));
  ABIP(accum_by_A)(self->data->A, w1, pr);
  for (int i = 0; i < m; i++) {
    pr[i] -= (t[i] + 1);
  }
  this_pr = ABIP(norm)(pr, m) / SQRTF(m);

  abip_float *dr = (abip_float *)abip_malloc(2 * m * sizeof(abip_float));
  for (int i = 0; i < m; i++) {
    dr[i] = y[i] - s2[i];
    dr[i + m] = y[i] + s1[i] - C;
  }
  this_dr = ABIP(norm)(dr, 2 * m) / (SQRTF(m) * C);

  ABIPMatrix *B = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
  ABIP(copy_A_matrix)(&B, self->data->A);
  B->n -= 1;
  abip_float *BTy = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memset(BTy, 0, n * sizeof(abip_float));
  ABIP(accum_by_Atrans)(B, y, BTy);

  r->dobj = 0;
  r->pobj = 0;
  for (int i = 0; i < m; i++) {
    r->dobj += y[i];
    r->pobj += C * xi[i];
  }

  r->pobj += 0.5 * ABIP(dot)(w1, w1, n);
  r->dobj -= 0.5 * ABIP(dot)(BTy, BTy, n);

  this_gap = r->dobj - r->pobj;
  this_gap = ABS(this_gap) / (1 + ABS(r->pobj));

  r->last_ipm_iter = ipm_iter;
  r->last_admm_iter = admm_iter;

  r->res_infeas = NAN;
  r->res_unbdd = NAN;

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

  ABIP(free_A_matrix)(B);
  abip_free(w1);
  abip_free(xi);
  abip_free(t);
  abip_free(y);
  abip_free(s1);
  abip_free(s2);
  abip_free(pr);
  abip_free(dr);
}

/**
@brief Formulate the qcp KKT matrix of the svm socp formulation
*/
cs *form_svm_kkt(svm *self) {
  abip_int n = self->n;
  abip_int m = self->m;
  cs *LTL;

  cs *N1 = cs_spalloc(m, n + 1, self->A->p[n + 1], 1, 0);
  memcpy(N1->i, self->A->i, self->A->p[n + 1] * sizeof(abip_int));
  memcpy(N1->p, self->A->p, (n + 2) * sizeof(abip_int));
  memcpy(N1->x, self->A->x, self->A->p[n + 1] * sizeof(abip_float));

  cs *N2 = cs_spalloc(m, n + 1, self->A->p[n + 1], 1, 0);
  memcpy(N2->i, self->A->i, self->A->p[n + 1] * sizeof(abip_int));
  memcpy(N2->p, self->A->p, (n + 2) * sizeof(abip_int));
  memcpy(N2->x, self->A->x, self->A->p[n + 1] * sizeof(abip_float));

  if (m > n + 1) {
    for (int i = 0; i < N1->p[n + 1]; i++) {
      N1->x[i] /= self->wF[N1->i[i]];
    }

    cs *diag = cs_spalloc(n + 1, n + 1, n + 1, 1, 1);

    for (int i = 0; i < n + 1; i++) {
      cs_entry(diag, i, i, 1 / self->wH[i]);
    }
    cs *diag_H = cs_compress(diag);
    cs_spfree(diag);

    LTL = cs_add(diag_H, cs_multiply(cs_transpose(N2, 1), N1), 1, 1);
    cs_spfree(diag_H);
  } else {
    for (int i = 0; i < N1->n; i++) {
      for (int j = N1->p[i]; j < N1->p[i + 1]; j++) {
        N1->x[j] *= self->wH[i];
      }
    }
    cs *diag = cs_spalloc(m, m, m, 1, 1);

    for (int i = 0; i < m; i++) {
      cs_entry(diag, i, i, self->wF[i]);
    }
    cs *diag_F = cs_compress(diag);
    cs_spfree(diag);

    LTL = cs_add(diag_F, cs_multiply(N1, cs_transpose(N2, 1)), 1, 1);

    cs_spfree(diag_F);
  }

  cs_spfree(N1);
  cs_spfree(N2);

  for (int i = 0; i < LTL->n; i++) {
    for (int j = LTL->p[i]; j < LTL->p[i + 1]; j++) {
      if (LTL->i[j] > i) LTL->x[j] = 0;
    }
  }
  cs_dropzeros(LTL);
  return LTL;
}

/**
@brief Initialize the preconditioner of conjugate gradient method for the svm
socp formulation
*/
void init_svm_precon(svm *self) {
  self->L->M = (abip_float *)abip_malloc(self->p * sizeof(abip_float));
  memset(self->L->M, 0, self->p * sizeof(abip_float));

  abip_float *M = self->L->M;

  M[0] = 1 / (self->stgs->rho_y + 1);

  for (int i = 0; i < self->A->p[self->A->n]; i++) {
    M[self->A->i[i] + 1] += 2 * POWF(self->A->x[i], 2);
  }

  for (int i = 1; i < self->m + 1; i++) {
    M[i] = 1 / (self->stgs->rho_y + M[i] + POWF(self->wB[i - 1], 2) +
                POWF(self->sc_D[i - 1], 2));
  }

  for (int i = self->m + 1; i < self->m + self->n + 1; i++) {
    M[i] = 1 / (self->stgs->rho_y + 2 * POWF(self->wE[i - self->m + 1], 2) +
                POWF(self->wD[i - self->m + 1], 2));
  }
}

/**
@brief Get the tolerance of the conjugate gradient method for the svm socp
formulation
*/
abip_float get_svm_pcg_tol(abip_int k, abip_float error_ratio,
                           abip_float norm_p) {
  if (k == -1) {
    return 1e-9 * norm_p;
  } else {
    if (error_ratio > 100000) {
      return MAX(1e-9, 3e-2 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 30000) {
      return MAX(1e-9, 3e-2 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 10000) {
      return MAX(1e-9, 3e-2 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 3000) {
      return MAX(1e-9, 2.5e-2 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 1000) {
      return MAX(1e-9, 2e-2 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 300) {
      return MAX(1e-9, 1.6e-2 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 100) {
      return MAX(1e-9, 1.3e-2 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 30) {
      return MAX(1e-9, 1e-2 * norm_p / POWF((k + 1), 2));
    } else if (error_ratio > 10) {
      return MAX(1e-9, 7e-3 * norm_p / POWF((k + 1), 2));
    } else {
      return MAX(1e-9, 4e-3 * norm_p / POWF((k + 1), 2));
    }
  }
}

/**
@brief Initialize the linear system solver work space for the svm socp
formulation
*/
abip_int init_svm_linsys_work(svm *self) {
  if (self->stgs->linsys_solver == 0) {  // mkl-dss need lower triangle
    self->L->K = cs_transpose(form_svm_kkt(self), 1);
  } else if (self->stgs->linsys_solver == 1) {  // qdldl need upper triangle
    self->L->K = form_svm_kkt(self);
  } else if (self->stgs->linsys_solver == 2) {  // cholesky need upper triangle
    self->L->K = form_svm_kkt(self);
  } else if (self->stgs->linsys_solver == 3) {  // pcg doesn't need kkt matrix
    init_svm_precon(self);
    self->L->K = ABIP_NULL;
  } else if (self->stgs->linsys_solver ==
             4) {  // mkl-pardiso need lower triangle
    self->L->K = cs_transpose(form_svm_kkt(self), 1);
  } else if (self->stgs->linsys_solver ==
             5) {  //  dense cholesky need upper triangle
    self->L->K = form_svm_kkt(self);
  } else {
    printf("\nlinsys solver type error\n");
    return -1;
  }
  return ABIP(init_linsys_work)(self);
}

/**
@brief Customized linear system solver for the svm socp formulation
*/
abip_int solve_svm_linsys(svm *self, abip_float *b, abip_float *pcg_warm_start,
                          abip_int iter, abip_float error_ratio) {
  ABIP(timer) linsys_timer;
  ABIP(tic)(&linsys_timer);

  ABIP(scale_array)(&b[self->p], -1, self->q);

  if (self->stgs->linsys_solver == 3) {  // pcg

    abip_int p = self->p;
    abip_float norm_p = ABIP(norm)(b, p);

    self->spe_A_times(self, &b[p], b);
    abip_float pcg_tol = get_svm_pcg_tol(iter, error_ratio, norm_p);
    abip_int cg_its = ABIP(solve_linsys)(self, b, p, pcg_warm_start, pcg_tol);

    if (iter >= 0) {
      self->L->total_cg_iters += cg_its;
    }
  } else {  // direct methods
    abip_int n = self->n;
    abip_int m = self->m;
    abip_float *b2 = (abip_float *)abip_malloc(self->p * sizeof(abip_float));
    memcpy(b2, b, self->p * sizeof(abip_float));
    self->spe_A_times(self, &b[self->p], b2);
    b[0] = b2[0] / (1 + self->stgs->rho_y);

    for (int i = 0; i < n; i++) {
      b[i + m + 1] = b2[i + m + 1] / self->wG[i];
    }

    abip_float *b3 = (abip_float *)abip_malloc(m * sizeof(abip_float));
    memcpy(b3, &b2[1], m * sizeof(abip_float));
    ABIP(scale_array)(b3, -1, m);
    ABIP(accum_by_A)(self->wX, &b[m + 1], b3);
    ABIP(scale_array)(b3, -1, m);

    abip_float *tmp = (abip_float *)abip_malloc((n + 1) * sizeof(abip_float));
    if (m > n + 1) {
      for (int i = 0; i < m; i++) {
        b3[i] /= self->wF[i];
      }
      memset(tmp, 0, (n + 1) * sizeof(abip_float));
      ABIP(accum_by_Atrans)(self->A, b3, tmp);
      ABIP(solve_linsys)(self, tmp, n + 1, ABIP_NULL, 0);

      abip_float *tmp2 = (abip_float *)abip_malloc(m * sizeof(abip_float));
      memset(tmp2, 0, m * sizeof(abip_float));

      ABIP(accum_by_A)(self->A, tmp, tmp2);
      for (int i = 0; i < m; i++) {
        tmp2[i] /= self->wF[i];
      }
      memcpy(&b[1], b3, m * sizeof(abip_float));
      ABIP(add_scaled_array)(&b[1], tmp2, m, -1);
      abip_free(tmp2);
    } else {
      ABIP(solve_linsys)(self, b3, m, ABIP_NULL, 0);
      memcpy(&b[1], b3, m * sizeof(abip_float));
    }

    memset(tmp, 0, (n + 1) * sizeof(abip_float));
    ABIP(accum_by_Atrans)(self->wX, &b[1], tmp);
    for (int i = 0; i < n; i++) {
      tmp[i] /= self->wG[i];
    }
    ABIP(add_scaled_array)(&b[m + 1], tmp, n, -1);

    abip_free(b2);
    abip_free(b3);
    abip_free(tmp);
  }

  ABIP(scale_array)(&b[self->p], -1, self->q);
  self->spe_AT_times(self, b, &b[self->p]);

  self->L->total_solve_time += ABIP(tocq)(&linsys_timer);

  return 0;
}

/**
@brief Free the linear system solver work space for the svm socp formulation
*/
void free_svm_linsys_work(svm *self) { ABIP(free_linsys)(self); }