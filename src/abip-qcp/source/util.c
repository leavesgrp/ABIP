#define _CRT_SECURE_NO_WARNINGS
#include "util.h"

#include "glbopts.h"
#include "linsys.h"

#if (defined NOTIMER)

void ABIP(tic)(ABIP(timer) * t) {}

abip_float ABIP(tocq)(ABIP(timer) * t) { return NAN; }

#elif (defined _WIN32 || _WIN64 || defined _WINDLL)

void ABIP(tic)(ABIP(timer) * t) {
  QueryPerformanceFrequency(&t->freq);
  QueryPerformanceCounter(&t->tic);
}

abip_float ABIP(tocq)(ABIP(timer) * t) {
  QueryPerformanceCounter(&t->toc);
  return (1e3 * (t->toc.QuadPart - t->tic.QuadPart) /
          (abip_float)t->freq.QuadPart);
}

#elif (defined __APPLE__)

void ABIP(tic)(ABIP(timer) * t) {
  /* read current clock cycles */
  t->tic = mach_absolute_time();
}

abip_float ABIP(tocq)(ABIP(timer) * t) {
  uint64_t duration;

  t->toc = mach_absolute_time();
  duration = t->toc - t->tic;

  mach_timebase_info(&(t->tinfo));
  duration *= t->tinfo.numer;
  duration /= t->tinfo.denom;

  return (abip_float)duration / 1e6;
}

#else

void ABIP(tic)(ABIP(timer) * t) { clock_gettime(CLOCK_MONOTONIC, &t->tic); }

abip_float ABIP(tocq)(ABIP(timer) * t) {
  struct timespec temp;

  clock_gettime(CLOCK_MONOTONIC, &t->toc);

  if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) {
    temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
    temp.tv_nsec = 1e9 + t->toc.tv_nsec - t->tic.tv_nsec;
  } else {
    temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
    temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
  }

  return (abip_float)temp.tv_sec * 1e3 + (abip_float)temp.tv_nsec / 1e6;
}

#endif

abip_float ABIP(toc)(ABIP(timer) * t) {
  abip_float time = ABIP(tocq)(t);
  abip_printf("time: %8.4f milli-seconds.\n", time);
  return time;
}

abip_float ABIP(str_toc)(char *str, ABIP(timer) * t) {
  abip_float time = ABIP(tocq)(t);
  abip_printf("%s - time: %8.4f milli-seconds.\n", str, time);
  return time;
}

void ABIP(print_work)(const ABIPWork *w) {
  abip_int i;
  abip_int l = w->n + w->m;

  abip_printf("\n u_t is \n");
  for (i = 0; i < l; i++) {
    abip_printf("%f\n", w->u_t[i]);
  }

  abip_printf("\n u is \n");
  for (i = 0; i < l; i++) {
    abip_printf("%f\n", w->u[i]);
  }

  abip_printf("\n v is \n");
  for (i = 0; i < l; i++) {
    abip_printf("%f\n", w->v[i]);
  }
}

void ABIP(print_data)(const ABIPData *d) {
  abip_printf("m = %i\n", (int)d->m);
  abip_printf("n = %i\n", (int)d->n);

  abip_printf("max_ipm_iters = %i\n", (int)d->stgs->max_ipm_iters);
  abip_printf("max_admm_iters = %i\n", (int)d->stgs->max_admm_iters);

  abip_printf("verbose = %i\n", (int)d->stgs->verbose);
  abip_printf("normalize = %i\n", (int)d->stgs->normalize);

  abip_printf("eps_p = %4f\n", d->stgs->eps_p);
  abip_printf("eps_d = %4f\n", d->stgs->eps_d);
  abip_printf("eps_g = %4f\n", d->stgs->eps_g);
  abip_printf("eps_inf = %4f\n", d->stgs->eps_inf);
  abip_printf("eps_unb = %4f\n", d->stgs->eps_unb);
  abip_printf("alpha = %4f\n", d->stgs->alpha);
  abip_printf("rho_y = %4f\n", d->stgs->rho_y);
}

void ABIP(print_array)(const abip_float *arr, abip_int n, const char *name) {
  abip_int i;
  abip_int j;
  abip_int k = 0;

  abip_int num_on_one_line = 10;

  abip_printf("\n");
  for (i = 0; i < n / num_on_one_line; ++i) {
    for (j = 0; j < num_on_one_line; ++j) {
      abip_printf("%s[%li] = %4f, ", name, (long)k, arr[k]);
      k++;
    }
    abip_printf("\n");
  }

  for (j = k; j < n; ++j) {
    abip_printf("%s[%li] = %4f, ", name, (long)j, arr[j]);
  }

  abip_printf("\n");
}

void ABIP(free_info)(ABIPInfo *info) {
  if (info) {
    abip_free(info);
  }
}

void ABIP(free_cone)(ABIPCone *k) {
  if (k) {
    if (k->q) {
      abip_free(k->q);
    }
    if (k->rq) {
      abip_free(k->rq);
    }
    abip_free(k);
  }
}

void ABIP(free_data)(ABIPData *d) {
  if (d) {
    if (d->b) {
      abip_free(d->b);
    }

    if (d->c) {
      abip_free(d->c);
    }

    if (d->stgs) {
      abip_free(d->stgs);
    }

    if (d->A) {
      ABIP(free_A_matrix)(d->A);
    }

    abip_free(d);
  }
}

void ABIP(free_sol)(ABIPSolution *sol) {
  if (sol) {
    if (sol->x) {
      abip_free(sol->x);
    }

    if (sol->y) {
      abip_free(sol->y);
    }

    if (sol->s) {
      abip_free(sol->s);
    }

    abip_free(sol);
  }
}

/**
@brief Default parameter settings
*/
void ABIP(set_default_settings)(ABIPData *d) {
  abip_int n = d->n;
  abip_int m = d->m;
  abip_int nz = d->A == ABIP_NULL ? 0 : d->A->p[n];
  abip_float sparsity = (abip_float)nz / (m * n);

  d->stgs->normalize = 1;
  d->stgs->scale_E = 1;
  d->stgs->scale_bc = 1;
  d->stgs->max_ipm_iters = MAX_IPM_ITERS;
  d->stgs->max_admm_iters = MAX_ADMM_ITERS;
  d->stgs->eps = EPS;
  d->stgs->eps_p = EPS;
  d->stgs->eps_d = EPS;
  d->stgs->eps_g = EPS;
  d->stgs->eps_inf = EPS;
  d->stgs->eps_unb = EPS;
  d->stgs->alpha = ALPHA;
  d->stgs->cg_rate = CG_RATE;

  d->stgs->use_indirect = 0;

  d->stgs->scale = SCALE;
  d->stgs->rho_y = 1e-6;
  d->stgs->rho_x = 1;
  d->stgs->rho_tau = 1;
  d->stgs->verbose = VERBOSE;

  d->stgs->err_dif =
      0;  // tol between max(dres,pres,dgap) of two consecutive inters

  d->stgs->inner_check_period = 500;
  d->stgs->outer_check_period = 1;

  // 0:mkl_dss, 1:qdldl, 2:sparse cholesky, 3:pcg,  4:pardiso, 5:dense cholesky
  if (m * n > 1e12) {
    d->stgs->linsys_solver = 3;
  } else if (sparsity > 0.4) {
    d->stgs->linsys_solver = 5;
  } else {
    d->stgs->linsys_solver = 1;
  }

  d->stgs->prob_type = 3;          // 0:general_qp, 1:lasso, 2:svm, 3:QCP
  d->stgs->time_limit = INFINITY;  // in s

  d->stgs->psi = 1;  // for qp&socp
  // d->stgs->psi = 1.5; //for ml

  d->stgs->origin_scaling = 1;
  d->stgs->ruiz_scaling = 1;
  d->stgs->pc_scaling = 0;
}
