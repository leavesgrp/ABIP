#include "qcp_config.h"
#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)

/**
@brief Initialize the qcp problem structure
*/
abip_int init_qcp(qcp **self, ABIPData *d, ABIPSettings *stgs) {
  qcp *this_qcp = (qcp *)abip_malloc(sizeof(qcp));
  *self = this_qcp;

  this_qcp->m = d->m;
  this_qcp->n = d->n;
  this_qcp->p = d->m;
  this_qcp->q = d->n;
  abip_int m = this_qcp->p;
  abip_int n = this_qcp->q;

  if (d->A == ABIP_NULL) {
    this_qcp->sparsity = 0;
  } else {
    this_qcp->sparsity = ((d->A->p[d->n] / (d->m * d->n)) < 0.05);
  }

  // non-identity DR scaling
  this_qcp->rho_dr =
      (abip_float *)abip_malloc((m + n + 1) * sizeof(abip_float));
  for (int i = 0; i < m + n + 1; i++) {
    if (i < m) {
      this_qcp->rho_dr[i] = stgs->rho_y;
    } else if (i < m + n) {
      this_qcp->rho_dr[i] = stgs->rho_x;
    } else {
      this_qcp->rho_dr[i] = stgs->rho_tau;
    }
  }

  this_qcp->L = (ABIPLinSysWork *)abip_malloc(sizeof(ABIPLinSysWork));
  this_qcp->pro_type = QCP;
  this_qcp->stgs = stgs;
  this_qcp->data = d;
  if (d->A != ABIP_NULL) {
    this_qcp->A = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
  }
  if (d->Q != ABIP_NULL) {
    this_qcp->Q = (ABIPMatrix *)abip_malloc(sizeof(ABIPMatrix));
  }
  if (d->b != ABIP_NULL) {
    this_qcp->b = (abip_float *)abip_malloc(this_qcp->p * sizeof(abip_float));
  }
  this_qcp->c = (abip_float *)abip_malloc(this_qcp->q * sizeof(abip_float));
  this_qcp->D = (abip_float *)abip_malloc(this_qcp->m * sizeof(abip_float));
  this_qcp->E = (abip_float *)abip_malloc(this_qcp->n * sizeof(abip_float));

  this_qcp->scaling_data = &scaling_qcp_data;
  this_qcp->un_scaling_sol = &un_scaling_qcp_sol;
  this_qcp->calc_residuals = &calc_qcp_residuals;
  this_qcp->init_spe_linsys_work = &init_qcp_linsys_work;
  this_qcp->solve_spe_linsys = &solve_qcp_linsys;
  this_qcp->free_spe_linsys_work = &free_qcp_linsys_work;
  this_qcp->spe_A_times = &qcp_A_times;
  this_qcp->spe_AT_times = &qcp_AT_times;
  this_qcp->inner_conv_check = &qcp_inner_conv_check;

  return 0;
}

/**
@brief Matrix-vector multiplication for the general qcp problem with A
untransposed
*/
void qcp_A_times(qcp *self, const abip_float *x, abip_float *y) {
  if (self->A != ABIP_NULL) {
    ABIP(accum_by_A)(self->A, x, y);
  }
}

/**
@brief Matrix-vector multiplication for the general qcp problem with A
transposed
*/
void qcp_AT_times(qcp *self, const abip_float *x, abip_float *y) {
  if (self->A != ABIP_NULL) {
    ABIP(accum_by_Atrans)(self->A, x, y);
  }
}

/**
@brief Scale the data for the qcp problem
*/
void scaling_qcp_data(qcp *self, ABIPCone *k) {
  if (self->data->b == ABIP_NULL) {
    self->b = ABIP_NULL;
  } else {
    memcpy(self->b, self->data->b, self->m * sizeof(abip_float));
  }

  memcpy(self->c, self->data->c, self->n * sizeof(abip_float));

  if (self->data->A == ABIP_NULL) {
    self->A = ABIP_NULL;
  } else if (!ABIP(copy_A_matrix)(&(self->A), self->data->A)) {
    abip_printf("ERROR: copy A matrix failed\n");
    RETURN;
  }

  if (self->data->Q == ABIP_NULL) {
    self->Q = ABIP_NULL;
  } else if (!ABIP(copy_A_matrix)(&(self->Q), self->data->Q)) {
    abip_printf("ERROR: copy Q matrix failed\n");
    RETURN;
  }

  abip_int m = self->m;
  abip_int n = self->n;
  ABIPMatrix *A = self->A;
  ABIPMatrix *Q = self->Q;

  abip_float min_row_scale = MIN_SCALE * SQRTF((abip_float)n);
  abip_float max_row_scale = MAX_SCALE * SQRTF((abip_float)n);
  abip_float min_col_scale = MIN_SCALE * SQRTF((abip_float)m);
  abip_float max_col_scale = MAX_SCALE * SQRTF((abip_float)m);

  abip_float *E_hat = self->E;
  abip_float *D_hat = self->D;

  for (int i = 0; i < n; i++) {
    E_hat[i] = 1;
  }
  for (int i = 0; i < m; i++) {
    D_hat[i] = 1;
  }

  abip_float *E = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memset(E, 0, n * sizeof(abip_float));

  abip_float *E1 = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memset(E1, 0, n * sizeof(abip_float));

  abip_float *E2 = (abip_float *)abip_malloc(n * sizeof(abip_float));
  memset(E2, 0, n * sizeof(abip_float));

  abip_float *D = (abip_float *)abip_malloc(m * sizeof(abip_float));
  memset(D, 0, m * sizeof(abip_float));

  abip_int origin_scaling = self->stgs->origin_scaling;
  abip_int ruiz_scaling = self->stgs->ruiz_scaling;
  abip_int pc_scaling = self->stgs->pc_scaling;
  abip_int count;
  abip_float mean_E;

  if (A == ABIP_NULL && Q == ABIP_NULL) {
    origin_scaling = 0;
    ruiz_scaling = 0;
    pc_scaling = 0;
  }

  if (ruiz_scaling) {
    abip_int n_ruiz = 10;

    for (int ruiz_iter = 0; ruiz_iter < n_ruiz; ruiz_iter++) {
      count = 0;
      memset(E, 0, n * sizeof(abip_float));
      memset(E1, 0, n * sizeof(abip_float));
      memset(E2, 0, n * sizeof(abip_float));
      memset(D, 0, m * sizeof(abip_float));

      if (A != ABIP_NULL) {
        for (int j = 0; j < n; j++) {
          if (A->p[j] == A->p[j + 1]) {
            E1[j] = 0;
          } else {
            E1[j] =
                SQRTF(ABIP(norm_inf)(&A->x[A->p[j]], A->p[j + 1] - A->p[j]));
          }
        }
      }

      if (Q != ABIP_NULL) {
        for (int j = 0; j < n; j++) {
          if (Q->p[j] == Q->p[j + 1]) {
            E2[j] = 0;
          } else {
            E2[j] =
                SQRTF(ABIP(norm_inf)(&Q->x[Q->p[j]], Q->p[j + 1] - Q->p[j]));
          }
        }
      }

      for (int i = 0; i < n; i++) {
        E[i] = E1[i] < E2[i] ? E2[i] : E1[i];
      }

      if (k->q) {
        for (int i = 0; i < k->qsize; i++) {
          mean_E = ABIP(vec_mean)(&E[count], k->q[i]);
          for (int j = 0; j < k->q[i]; j++) {
            E[j + count] = mean_E;
          }
          count += k->q[i];
        }
      }

      if (k->rq) {
        for (int i = 0; i < k->rqsize; i++) {
          mean_E = ABIP(vec_mean)(&E[count], k->rq[i]);
          for (int j = 0; j < k->rq[i]; j++) {
            E[j + count] = mean_E;
          }
          count += k->rq[i];
        }
      }

      if (A != ABIP_NULL) {
        for (int i = 0; i < A->p[n]; i++) {
          if (D[A->i[i]] < ABS(A->x[i])) {
            D[A->i[i]] = ABS(A->x[i]);
          }
        }
        for (int i = 0; i < m; i++) {
          D[i] = SQRTF(D[i]);
          if (D[i] < min_row_scale)
            D[i] = 1;
          else if (D[i] > max_row_scale)
            D[i] = max_row_scale;
        }

        for (int i = 0; i < n; i++) {
          if (E[i] < min_col_scale)
            E[i] = 1;
          else if (E[i] > max_col_scale)
            E[i] = max_col_scale;
          for (int j = A->p[i]; j < A->p[i + 1]; j++) {
            A->x[j] /= E[i];
          }
        }
      }

      if (Q != ABIP_NULL) {
        for (int i = 0; i < n; i++) {
          for (int j = Q->p[i]; j < Q->p[i + 1]; j++) {
            Q->x[j] /= E[i];
          }
        }
        for (int i = 0; i < Q->p[n]; i++) {
          Q->x[i] /= E[Q->i[i]];
        }
      }

      if (A != ABIP_NULL) {
        for (int i = 0; i < A->p[n]; i++) {
          A->x[i] /= D[A->i[i]];
        }
      }

      for (int i = 0; i < n; i++) {
        E_hat[i] *= E[i];
      }

      for (int i = 0; i < m; i++) {
        D_hat[i] *= D[i];
      }
    }
  }

  if (origin_scaling) {
    memset(E, 0, n * sizeof(abip_float));
    memset(E1, 0, n * sizeof(abip_float));
    memset(E2, 0, n * sizeof(abip_float));
    memset(D, 0, m * sizeof(abip_float));

    count = 0;

    if (A != ABIP_NULL) {
      for (int i = 0; i < n; i++) {
        for (int j = A->p[i]; j < A->p[i + 1]; j++) {
          E1[i] += A->x[j] * A->x[j];
        }
        E1[i] = SQRTF(E1[i]);
      }
    }

    if (Q != ABIP_NULL) {
      for (int i = 0; i < n; i++) {
        for (int j = Q->p[i]; j < Q->p[i + 1]; j++) {
          E2[i] += Q->x[j] * Q->x[j];
        }
        E2[i] = SQRTF(E2[i]);
      }
    }

    for (int i = 0; i < n; i++) {
      E[i] = SQRTF(E1[i] < E2[i] ? E2[i] : E1[i]);
    }

    if (k->q) {
      for (int i = 0; i < k->qsize; i++) {
        mean_E = ABIP(vec_mean)(&E[count], k->q[i]);
        for (int j = 0; j < k->q[i]; j++) {
          E[j + count] = mean_E;
        }
        count += k->q[i];
      }
    }

    if (k->rq) {
      for (int i = 0; i < k->rqsize; i++) {
        mean_E = ABIP(vec_mean)(&E[count], k->rq[i]);
        for (int j = 0; j < k->rq[i]; j++) {
          E[j + count] = mean_E;
        }
        count += k->rq[i];
      }
    }

    if (A != ABIP_NULL) {
      for (int i = 0; i < A->p[n]; i++) {
        D[A->i[i]] += A->x[i] * A->x[i];
      }
      for (int i = 0; i < m; i++) {
        D[i] = SQRTF(SQRTF(D[i]));
        if (D[i] < min_row_scale)
          D[i] = 1;
        else if (D[i] > max_row_scale)
          D[i] = max_row_scale;
      }

      for (int i = 0; i < n; i++) {
        if (E[i] < min_col_scale)
          E[i] = 1;
        else if (E[i] > max_col_scale)
          E[i] = max_col_scale;
        for (int j = A->p[i]; j < A->p[i + 1]; j++) {
          A->x[j] /= E[i];
        }
      }
    }

    if (Q != ABIP_NULL) {
      for (int i = 0; i < n; i++) {
        for (int j = Q->p[i]; j < Q->p[i + 1]; j++) {
          Q->x[j] /= E[i];
        }
      }
      for (int i = 0; i < Q->p[n]; i++) {
        Q->x[i] /= E[Q->i[i]];
      }
    }

    if (A != ABIP_NULL) {
      for (int i = 0; i < A->p[n]; i++) {
        A->x[i] /= D[A->i[i]];
      }
    }

    for (int i = 0; i < n; i++) {
      E_hat[i] *= E[i];
    }

    for (int i = 0; i < m; i++) {
      D_hat[i] *= D[i];
    }
  }

  if (pc_scaling) {
    memset(E, 0, n * sizeof(abip_float));
    memset(E1, 0, n * sizeof(abip_float));
    memset(E2, 0, n * sizeof(abip_float));
    memset(D, 0, m * sizeof(abip_float));
    count = 0;
    abip_float alpha_pc = 1;
    if (A != ABIP_NULL) {
      for (int i = 0; i < n; i++) {
        for (int j = A->p[i]; j < A->p[i + 1]; j++) {
          E1[i] += POWF(ABS(A->x[j]), alpha_pc);
        }
        E1[i] = SQRTF(POWF(E1[i], 1 / alpha_pc));
      }
    }
    if (Q != ABIP_NULL) {
      for (int i = 0; i < n; i++) {
        for (int j = Q->p[i]; j < Q->p[i + 1]; j++) {
          E2[i] += POWF(ABS(Q->x[j]), alpha_pc);
        }
        E2[i] = SQRTF(POWF(E2[i], 1 / alpha_pc));
      }
    }

    for (int i = 0; i < n; i++) {
      E[i] = E1[i] < E2[i] ? E2[i] : E1[i];
    }

    if (k->q) {
      for (int i = 0; i < k->qsize; i++) {
        mean_E = ABIP(vec_mean)(&E[count], k->q[i]);
        for (int j = 0; j < k->q[i]; j++) {
          E[j + count] = mean_E;
        }
        count += k->q[i];
      }
    }

    if (k->rq) {
      for (int i = 0; i < k->rqsize; i++) {
        mean_E = ABIP(vec_mean)(&E[count], k->rq[i]);
        for (int j = 0; j < k->rq[i]; j++) {
          E[j + count] = mean_E;
        }
        count += k->rq[i];
      }
    }

    if (A != ABIP_NULL) {
      for (int i = 0; i < A->p[n]; i++) {
        D[A->i[i]] += POWF(ABS(A->x[i]), 2 - alpha_pc);
      }
      for (int i = 0; i < m; i++) {
        D[i] = SQRTF(POWF(D[i], 1 / (2 - alpha_pc)));
        if (D[i] < min_row_scale)
          D[i] = 1;
        else if (D[i] > max_row_scale)
          D[i] = max_row_scale;
      }

      for (int i = 0; i < n; i++) {
        if (E[i] < min_col_scale)
          E[i] = 1;
        else if (E[i] > max_col_scale)
          E[i] = max_col_scale;
        for (int j = A->p[i]; j < A->p[i + 1]; j++) {
          A->x[j] /= E[i];
        }
      }
    }

    if (Q != ABIP_NULL) {
      for (int i = 0; i < n; i++) {
        for (int j = Q->p[i]; j < Q->p[i + 1]; j++) {
          Q->x[j] /= E[i];
        }
      }
      for (int i = 0; i < Q->p[n]; i++) {
        Q->x[i] /= E[Q->i[i]];
      }
    }

    if (A != ABIP_NULL) {
      for (int i = 0; i < A->p[n]; i++) {
        A->x[i] /= D[A->i[i]];
      }
    }

    for (int i = 0; i < n; i++) {
      E_hat[i] *= E[i];
    }

    for (int i = 0; i < m; i++) {
      D_hat[i] *= D[i];
    }
  }

  abip_float sc =
      SQRTF(SQRTF(ABIP(norm_sq)(self->c, n) + ABIP(norm_sq)(self->b, m)));

  if (self->b != ABIP_NULL) {
    for (int i = 0; i < m; i++) {
      self->b[i] /= D_hat[i];
    }
  }

  for (int i = 0; i < n; i++) {
    self->c[i] /= E_hat[i];
  }

  if (sc < MIN_SCALE)
    sc = 1;
  else if (sc > MAX_SCALE)
    sc = MAX_SCALE;
  self->sc_b = 1 / sc;
  self->sc_c = 1 / sc;

  if (self->b != ABIP_NULL) {
    ABIP(scale_array)(self->b, self->sc_b * self->stgs->scale, m);
  }
  ABIP(scale_array)(self->c, self->sc_c * self->stgs->scale, n);

  abip_free(E);
  abip_free(E1);
  abip_free(E2);
  abip_free(D);
}

/**
@brief Get the unscaled solution of the general qcp problem
*/
void un_scaling_qcp_sol(qcp *self, ABIPSolution *sol) {
  abip_int i;

  abip_float *D = self->D;
  abip_float *E = self->E;

  for (i = 0; i < self->n; ++i) {
    sol->x[i] /= (E[i] * self->sc_b);
  }

  for (i = 0; i < self->m; ++i) {
    sol->y[i] /= (D[i] * self->sc_c);
  }

  for (i = 0; i < self->n; ++i) {
    sol->s[i] *= E[i] / (self->sc_c * self->stgs->scale);
  }
}

/**
@brief Check whether the inner loop of the genral qcp problem has converged
*/
abip_float qcp_inner_conv_check(qcp *self, ABIPWork *w) {
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
      (1 + ABIP(norm)(Qu, m + n + 1) + ABIP(norm)(w->v_origin, m + n + 1));

  abip_free(Qu);
  abip_free(Mu);
  abip_free(tem);

  return error_inner;
}

/**
@brief Calculate the residuals of the general qcp problem
*/
void calc_qcp_residuals(qcp *self, ABIPWork *w, ABIPResiduals *r,
                        abip_int ipm_iter, abip_int admm_iter) {
  DEBUG_FUNC

  abip_int n = w->n;
  abip_int m = w->m;

  abip_float *y = (abip_float *)abip_malloc(m * sizeof(abip_float));
  abip_float *x = (abip_float *)abip_malloc(n * sizeof(abip_float));
  abip_float *s = (abip_float *)abip_malloc(n * sizeof(abip_float));

  abip_float this_pr;
  abip_float this_dr;
  abip_float this_gap;

  if (admm_iter && r->last_admm_iter == admm_iter) {
    RETURN;
  }

  r->last_ipm_iter = ipm_iter;
  r->last_admm_iter = admm_iter;

  r->tau = ABS(w->u[n + m]);
  r->kap =
      ABS(w->v_origin[n + m]) /
      (self->stgs->normalize ? (self->stgs->scale * self->sc_c * self->sc_b)
                             : 1);

  memcpy(y, w->u, m * sizeof(abip_float));
  memcpy(x, &w->u[m], n * sizeof(abip_float));
  memcpy(s, &w->v_origin[m], n * sizeof(abip_float));

  ABIP(scale_array)(y, 1 / r->tau, m);
  ABIP(scale_array)(x, 1 / r->tau, n);
  ABIP(scale_array)(s, 1 / r->tau, n);

  abip_float *Ax = (abip_float *)abip_malloc(m * sizeof(abip_float));
  abip_float *Ax_b = (abip_float *)abip_malloc(m * sizeof(abip_float));

  memset(Ax, 0, m * sizeof(abip_float));
  self->spe_A_times(self, x, Ax);

  memcpy(Ax_b, Ax, m * sizeof(abip_float));
  ABIP(add_scaled_array)(Ax_b, self->b, m, -1);

  // for scs cg tol
  r->Ax_b_norm = ABIP(norm_inf)(Ax_b, m);

  ABIP(c_dot)(Ax, self->D, m);
  ABIP(c_dot)(Ax_b, self->D, m);

  this_pr = ABIP(norm_inf)(Ax_b, m) /
            (self->sc_b + MAX(ABIP(norm_inf)(Ax, m), self->sc_b * w->nm_inf_b));

  abip_float *Qx = (abip_float *)abip_malloc(n * sizeof(abip_float));
  abip_float *ATy = (abip_float *)abip_malloc(n * sizeof(abip_float));
  abip_float *Qx_ATy_c_s = (abip_float *)abip_malloc(n * sizeof(abip_float));

  memset(Qx, 0, n * sizeof(abip_float));
  abip_float xQx_2 = 0;

  if (self->Q != ABIP_NULL) {
    ABIP(accum_by_A)(self->Q, x, Qx);
    xQx_2 = ABIP(dot)(x, Qx, n) / (2 * self->sc_b * self->sc_c);
  }

  memset(ATy, 0, n * sizeof(abip_float));
  self->spe_AT_times(self, y, ATy);

  memcpy(Qx_ATy_c_s, Qx, n * sizeof(abip_float));
  ABIP(add_scaled_array)(Qx_ATy_c_s, ATy, n, -1);
  ABIP(add_scaled_array)(Qx_ATy_c_s, self->c, n, 1);
  ABIP(add_scaled_array)(Qx_ATy_c_s, s, n, -1);

  r->Qx_ATy_c_s_norm = ABIP(norm_inf)(Qx_ATy_c_s, n);

  ABIP(c_dot)(Qx, self->E, n);
  ABIP(c_dot)(ATy, self->E, n);
  ABIP(c_dot)(Qx_ATy_c_s, self->E, n);
  ABIP(c_dot)(s, self->E, n);

  this_dr = ABIP(norm_inf)(Qx_ATy_c_s, n) /
            (self->sc_c + MAX(self->sc_c * w->nm_inf_c, ABIP(norm_inf)(Qx, n)));

  abip_float cTx = ABIP(dot)(self->c, x, n) / (self->sc_b * self->sc_c);
  abip_float bTy = ABIP(dot)(self->b, y, m) / (self->sc_b * self->sc_c);

  this_gap = ABS(2 * xQx_2 + cTx - bTy) /
             (1 + MAX(2 * xQx_2, MAX(ABS(cTx), ABS(bTy))));

  r->pobj = xQx_2 + cTx;
  r->dobj = -xQx_2 + bTy;

  r->res_dif = MAX(MAX(ABS(this_pr - r->res_pri), ABS(this_dr - r->res_dual)),
                   ABS(this_gap - r->rel_gap));
  r->res_pri = this_pr;
  r->res_dual = this_dr;
  r->rel_gap = this_gap;
  r->error_ratio =
      MAX(r->res_pri / self->stgs->eps_p,
          MAX(r->res_dual / self->stgs->eps_d, r->rel_gap / self->stgs->eps_g));

  if (ABIP(dot)(self->c, &w->u[m], n) < 0) {
    ABIP(scale_array)(Qx, r->tau, n);
    ABIP(scale_array)(Ax, r->tau, m);
    r->res_unbdd = MAX(ABIP(norm)(Qx, n), ABIP(norm)(Ax, m)) /
                   (-ABIP(dot)(self->c, &w->u[m], n));
  } else {
    r->res_unbdd = INFINITY;
  }

  if (ABIP(dot)(self->b, w->u, m) > 0) {
    ABIP(scale_array)(ATy, r->tau, n);
    ABIP(scale_array)(s, r->tau, n);
    ABIP(add_scaled_array)(ATy, s, n, 1);

    r->res_infeas = ABIP(norm)(ATy, n) / ABIP(dot)(self->b, w->u, m);
  } else {
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

/**
@brief Formulate the qcp KKT matrix of the general qcp problem
*/
/*K =  -rho_dr(1:m)*I           -A
                        Q + rho_dr(m+1:m+n)*I
*/
cs *form_qcp_kkt(qcp *self) {
  abip_int m = self->m;
  abip_int n = self->n;
  ABIPMatrix *A = self->A;
  ABIPMatrix *Q = self->Q;

  abip_int nnzA = A == ABIP_NULL ? 0 : A->p[A->n];
  ;
  abip_int nnzQ = Q == ABIP_NULL ? 0 : Q->p[Q->n];
  abip_int nnzK = m + nnzA + nnzQ + n;
  abip_int i;
  abip_int j;
  abip_float tem;

  cs *K = cs_spalloc(m + n, m + n, nnzK, 1, 1);

  for (i = 0; i < m; i++) {
    cs_entry(K, i, i, -self->rho_dr[i]);
  }

  if (A != ABIP_NULL) {
    for (i = 0; i < n; i++) {
      for (j = A->p[i]; j < A->p[i + 1]; j++) {
        cs_entry(K, A->i[j], m + i, -A->x[j]);
      }
    }
  }

  for (i = 0; i < n; i++) {
    if (Q == ABIP_NULL || Q->p[i] == Q->p[i + 1]) {
      cs_entry(K, m + i, m + i, self->rho_dr[m + i]);
    } else {
      for (j = Q->p[i]; j < Q->p[i + 1]; j++) {
        if (Q->i[j] > i)
          tem = 0;
        else if (Q->i[j] == i)
          tem = Q->x[j] + self->rho_dr[m + i];
        else
          tem = Q->x[j];
        cs_entry(K, m + Q->i[j], m + i, tem);
      }
    }
  }

  cs *K_csc = cs_compress(K);
  cs_spfree(K);
  cs_dropzeros(K_csc);

  return K_csc;
}

/**
@brief Initialize the preconditioner of conjugate gradient method for the
general qcp problem
*/
void init_qcp_precon(qcp *self) {
  self->L->M = (abip_float *)abip_malloc(self->q * sizeof(abip_float));
  memset(self->L->M, 0, self->q * sizeof(abip_float));

  for (int i = 0; i < self->q; i++) {
    for (int j = self->A->p[i]; j < self->A->p[i + 1]; j++) {
      self->L->M[i] +=
          self->A->x[j] * self->A->x[j] / self->rho_dr[self->A->i[j]];
    }
  }

  if (self->Q != ABIP_NULL) {
    for (int i = 0; i < self->q; i++) {
      for (int j = self->Q->p[i]; j < self->Q->p[i + 1]; j++) {
        if (i == self->Q->i[j]) {
          self->L->M[i] += self->Q->x[j];
          break;
        }
      }
    }
  }

  ABIP(add_scaled_array)(self->L->M, &self->rho_dr[self->p], self->q, 1);
  for (int i = 0; i < self->q; i++) {
    self->L->M[i] = 1.0 / self->L->M[i];
  }
}

/**
@brief Get the tolerance of the conjugate gradient method for the general qcp
problem
*/
abip_float get_qcp_pcg_tol(abip_int k, abip_float error_ratio,
                           abip_float norm_p) {
  if (k == -1) {
    return 1e-9 * norm_p;
  } else {
    return MAX(1e-9, 1e-5 * norm_p / POWF((k + 1), 2));
  }
}

/**
@brief Initialize the linear system solver work space for the general qcp
problem
*/
abip_int init_qcp_linsys_work(qcp *self) {
  if (self->stgs->linsys_solver == 0) {  // mkl-dss need lower triangle
    self->L->K = cs_transpose(form_qcp_kkt(self), 1);
  } else if (self->stgs->linsys_solver == 1) {  // qdldl need upper triangle
    self->L->K = form_qcp_kkt(self);
  } else if (self->stgs->linsys_solver == 2) {  // cholesky need upper triangle
    self->L->K = form_qcp_kkt(self);
  } else if (self->stgs->linsys_solver == 3) {  // pcg doesn't need kkt matrix
    init_qcp_precon(self);
    self->L->K = ABIP_NULL;
  } else if (self->stgs->linsys_solver ==
             4) {  // mkl-pardiso need lower triangle
    self->L->K = cs_transpose(form_qcp_kkt(self), 1);
  } else if (self->stgs->linsys_solver ==
             5) {  //  dense cholesky need upper triangle
    self->L->K = form_qcp_kkt(self);
  } else {
    printf("\nlinsys solver type error\n");
    return -1;
  }

  return ABIP(init_linsys_work)(self);
}

/**
@brief Linear system solver for the general qcp problem
*/
abip_int solve_qcp_linsys(qcp *self, abip_float *b, abip_float *pcg_warm_start,
                          abip_int iter, abip_float error_ratio) {
  ABIP(timer) linsys_timer;
  ABIP(tic)(&linsys_timer);

  if (self->stgs->linsys_solver == 3) {  // pcg

    abip_int n = self->n;
    abip_int m = self->m;
    abip_int i;

    abip_float norm_p = ABIP(norm)(&b[m], n);

    abip_float *tem = (abip_float *)abip_malloc(sizeof(abip_float) * m);
    memcpy(tem, b, m * sizeof(abip_float));
    for (i = 0; i < m; i++) {
      tem[i] /= self->rho_dr[i];
    }
    self->spe_AT_times(self, tem, &b[m]);

    abip_free(tem);

    abip_float pcg_tol = get_qcp_pcg_tol(iter, error_ratio, norm_p);

    abip_int cg_its;
    if (iter == -1) {
      cg_its = ABIP(solve_linsys)(self, &b[m], n, pcg_warm_start, error_ratio);
    } else {
      cg_its =
          ABIP(solve_linsys)(self, &b[m], n, &pcg_warm_start[m], error_ratio);
    }

    if (iter >= 0) {
      self->L->total_cg_iters += cg_its;
    }

    ABIP(scale_array)(b, -1, m);
    self->spe_A_times(self, &b[m], b);
    for (i = 0; i < m; i++) {
      b[i] /= -self->rho_dr[i];
    }

  } else {  // direct methods

    abip_int n = self->n + self->m;

    // difference here
    ABIP(scale_array)(b, -1, self->m);

    ABIP(solve_linsys)(self, b, n, ABIP_NULL, 0);
  }

  self->L->total_solve_time += ABIP(tocq)(&linsys_timer);

  return 0;
}

/**
@brief Free the linear system solver work space for the general qcp problem
*/
void free_qcp_linsys_work(qcp *self) { ABIP(free_linsys)(self); }
