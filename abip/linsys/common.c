#include "common.h"
#include "linsys.h"

#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)

abip_int ABIP(copy_A_matrix)
(
      ABIPMatrix **dstp, 
      const ABIPMatrix *src
) 
{
      abip_int Annz = src->p[src->n];
      ABIPMatrix *A = (ABIPMatrix *)abip_calloc(1, sizeof(ABIPMatrix));
      if (!A) 
      {
            return 0;
      }
      A->n = src->n;
      A->m = src->m;
      A->x = (abip_float *)abip_malloc(sizeof(abip_float) * Annz);
      A->i = (abip_int *)abip_malloc(sizeof(abip_int) * Annz);
      A->p = (abip_int *)abip_malloc(sizeof(abip_int) * (src->n + 1));
     
      if (!A->x || !A->i || !A->p) 
      {
            return 0;
      }
      
      memcpy(A->x, src->x, sizeof(abip_float) * Annz);
      memcpy(A->i, src->i, sizeof(abip_int) * Annz);
      memcpy(A->p, src->p, sizeof(abip_int) * (src->n + 1));
      
      *dstp = A;
      return 1;
}

abip_int ABIP(validate_lin_sys)
(
      const ABIPMatrix *A
) 
{
      abip_int i;  
      abip_int r_max;  
      abip_int Annz;
      
      if (!A->x || !A->i || !A->p) 
      {
            abip_printf("ERROR: incomplete data!\n");
            return -1;
      }
      
      for (i = 0; i < A->n; ++i) 
      {
            if (A->p[i] == A->p[i + 1]) 
            {
                  abip_printf("WARN: the %li-th column empty!\n", (long)i);
            } 
            else if (A->p[i] > A->p[i + 1]) 
            {
                  abip_printf("ERROR: the column pointers decreases!\n");
                  return -1;
            }
      }
      
      Annz = A->p[A->n];
      if (((abip_float)Annz / A->m > A->n) || (Annz <= 0)) 
      {
            abip_printf("ERROR: the number of nonzeros in A = %li, outside of valid range!\n", (long) Annz);
            return -1;
      }
      
      r_max = 0;
      for (i = 0; i < Annz; ++i) 
      {
            if (A->i[i] > r_max) 
            {
                  r_max = A->i[i];
            }
      }
      if (r_max > A->m - 1) 
      {
            abip_printf("ERROR: the number of rows in A is inconsistent with input dimension!\n");
            return -1;
      }

      return 0;
}

void ABIP(free_A_matrix)
(
      ABIPMatrix *A
) 
{
      if (A->x) 
      {
            abip_free(A->x);
      }
      if (A->i) 
      {
            abip_free(A->i);
      }
      if (A->p) 
      {
            abip_free(A->p);
      }
      
      abip_free(A);
}

#if EXTRA_VERBOSE > 0

static void print_A_matrix
(
      const ABIPMatrix *A
) 
{
      abip_int i;  
      abip_int j;

      if (A->p[A->n] < 2500) 
      {
            abip_printf("\n");
            for (i = 0; i < A->n; ++i) 
            {
                  abip_printf("Col %li: ", (long)i);
                  for (j = A->p[i]; j < A->p[i + 1]; j++) 
                  {
                          abip_printf("A[%li,%li] = %4f, ", (long)A->i[j], (long)i, A->x[j]);
                  }
                  abip_printf("norm col = %4f\n", ABIP(norm)(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
            }
            abip_printf("norm A = %4f\n", ABIP(norm)(A->x, A->p[A->n]));
      }
}
#endif

void ABIP(_normalize_A)
(
      ABIPMatrix *A, 
      const ABIPSettings *stgs, 
      ABIPScaling *scal
) 
{
      abip_float *D = (abip_float *)abip_malloc(A->m * sizeof(abip_float));
      abip_float *E = (abip_float *)abip_malloc(A->n * sizeof(abip_float));
      abip_float *Dt = (abip_float *)abip_malloc(A->m * sizeof(abip_float));
      abip_float *Et = (abip_float *)abip_malloc(A->n * sizeof(abip_float));
      abip_float *nms = (abip_float *)abip_calloc(A->m, sizeof(abip_float));

      abip_float min_row_scale = MIN_SCALE * SQRTF((abip_float)A->n); 
      abip_float max_row_scale = MAX_SCALE * SQRTF((abip_float)A->n);
      abip_float min_col_scale = MIN_SCALE * SQRTF((abip_float)A->m); 
      abip_float max_col_scale = MAX_SCALE * SQRTF((abip_float)A->m);
      
      abip_int i; 
      abip_int j;
      abip_int c1; 
      abip_int c2;

      abip_float wrk; 
      abip_float e;

      #if EXTRA_VERBOSE > 0
            ABIP(timer) normalize_timer;
            ABIP(tic)(&normalize_timer);
            abip_printf("normalizing A\n");
            print_A_matrix(A);
      #endif

      memset(D, 0, A->m * sizeof(abip_float));
      memset(E, 0, A->n * sizeof(abip_float));

      for (i = 0; i < A->n; ++i) 
      {
            c1 = A->p[i + 1] - A->p[i];
            e = ABIP(norm)(&(A->x[A->p[i]]), c1);
            if (e < min_col_scale) 
            {
                  e = 1;
            } 
            else if (e > max_col_scale) 
            {
                  e = max_col_scale;
            }
            ABIP(scale_array)(&(A->x[A->p[i]]), 1.0 / e, c1);
            E[i] = e;
      }

      for (i = 0; i < A->n; ++i) 
      {
            c1 = A->p[i];
            c2 = A->p[i + 1];
            for (j = c1; j < c2; ++j) 
            {
                  wrk = A->x[j];
                  D[A->i[j]] += wrk * wrk;
            }
      }

      for (i = 0; i < A->m; ++i) 
      {
            D[i] = SQRTF(D[i]);
            if (D[i] < min_row_scale) 
            {
                  D[i] = 1;
            } 
            else if (D[i] > max_row_scale) 
            {
                  D[i] = max_row_scale;
            }
      }

      for (i = 0; i < A->n; ++i) 
      {
            for (j = A->p[i]; j < A->p[i + 1]; ++j) 
            {
                  A->x[j] /= D[A->i[j]];
            }
      }

      for (i = 0; i < A->n; ++i) 
      {
            for (j = A->p[i]; j < A->p[i + 1]; ++j) 
            {
                  wrk = A->x[j];
                  nms[A->i[j]] += wrk * wrk;
            }
      }
      
      scal->mean_norm_row_A = 0.0;
      for (i = 0; i < A->m; ++i) 
      {
            scal->mean_norm_row_A += SQRTF(nms[i]) / A->m;
      }
      
      abip_free(nms);

      scal->mean_norm_col_A = 0.0;
      for (i = 0; i < A->n; ++i) 
      {
            c1 = A->p[i + 1] - A->p[i];
            scal->mean_norm_col_A+= ABIP(norm)(&(A->x[A->p[i]]), c1) / A->n;
      }

      if (stgs->scale != 1) 
      {
            ABIP(scale_array)(A->x, stgs->scale, A->p[A->n]);
      }

      scal->D = D; 
      scal->E = E;

      #if EXTRA_VERBOSE > 0
            abip_printf("finished normalizing A, time: %1.2e seconds. \n", ABIP(tocq)(&normalize_timer) / 1e3);
            print_A_matrix(A);
      #endif

}

void ABIP(_un_normalize_A)
(
      ABIPMatrix *A, 
      const ABIPSettings *stgs, 
      const ABIPScaling *scal
) 
{
      abip_int i; 
      abip_int j;
      
      abip_float *D = scal->D;
      abip_float *E = scal->E;
      
      for (i = 0; i < A->n; ++i) 
      {
            for (j = A->p[i]; j < A->p[i + 1]; ++j) 
            {
                  A->x[j] *= D[A->i[j]];
            }
      }

      for (i = 0; i < A->n; ++i) 
      {
            ABIP(scale_array)(&(A->x[A->p[i]]), E[i] / stgs->scale, A->p[i + 1] - A->p[i]);
      }
}

void ABIP(_accum_by_Atrans)
(
      abip_int n, 
      abip_float *Ax, 
      abip_int *Ai, 
      abip_int *Ap, 
      const abip_float *x, 
      abip_float *y
) 
{
      abip_int p; 
      abip_int j;
      
      abip_int c1; 
      abip_int c2;
      abip_float yj;

      #if EXTRA_VERBOSE > 0
            ABIP(timer) mult_by_Atrans_timer;
            ABIP(tic)(&mult_by_Atrans_timer);
      #endif
      
      #ifdef _OPENMP
            #pragma omp parallel for private(p, c1, c2, yj)
      #endif
      
      for (j = 0; j < n; j++) 
      {
            yj = y[j];
            c1 = Ap[j];
            c2 = Ap[j + 1];
            for (p = c1; p < c2; p++) 
            {
                  yj += Ax[p] * x[Ai[p]];
            }
            y[j] = yj;
      }

      #if EXTRA_VERBOSE > 0
            abip_printf("mult By A trans time: %1.2e seconds. \n", ABIP(tocq)(&mult_by_Atrans_timer) / 1e3);
      #endif
}

void ABIP(_accum_by_A)
(
      abip_int n, 
      abip_float *Ax, 
      abip_int *Ai, 
      abip_int *Ap, 
      const abip_float *x, 
      abip_float *y
) 
{
      abip_int p; 
      abip_int j;
  
      abip_int c1; 
      abip_int c2;
      abip_float xj;
      
      #if EXTRA_VERBOSE > 0
            ABIP(timer) mult_by_A_timer;
            ABIP(tic)(&mult_by_A_timer);
      #endif
      
      #ifdef _OPENMP
            #pragma omp parallel for private(p, c1, c2, xj)
            for (j = 0; j < n; j++) 
            {
                  xj = x[j];
                  c1 = Ap[j];
                  c2 = Ap[j + 1];
                  for (p = c1; p < c2; p++) 
                  {
                        #pragma omp atomic
                        y[Ai[p]] += Ax[p] * xj;
                  }
            }
      #endif

      for (j = 0; j < n; j++) 
      {
            xj = x[j];
            c1 = Ap[j];
            c2 = Ap[j + 1];
            for (p = c1; p < c2; p++) 
            {
                  y[Ai[p]] += Ax[p] * xj;
            }
      }

      #if EXTRA_VERBOSE > 0
            abip_printf("mult By A time: %1.2e seconds \n", ABIP(tocq)(&mult_by_A_timer) / 1e3);
      #endif
}

abip_float ABIP(cumsum)
(
      abip_int *p, 
      abip_int *c, 
      abip_int n
) 
{
      abip_int i; 
      abip_float nz = 0;
      abip_float nz2 = 0;
      
      if (!p || !c) 
      {
            return (-1);
      }
      
      for (i = 0; i < n; i++) 
      {
            p[i] = nz;
            nz += c[i];
            nz2 += c[i];
            c[i] = p[i];
      }
      
      p[n] = nz;
      return nz2;
}
