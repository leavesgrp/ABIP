#include "cs.h"
#include "abip.h"

/* NB: this is a subset of the routines in the CSPARSE package by Tim Davis et. al., for the full package please visit
 * http://www.cise.ufl.edu/research/sparse/CSparse/ */

/* wrapper for free */
static void *cs_free
(
      void *p
) 
{
      if (p) 
      {
            abip_free(p);              
      }                  
      return (ABIP_NULL);       
}

static cs *cs_done
(
      cs *C, 
      void *w, 
      void *x, 
      abip_int ok
) 
{
      cs_free(w);                                                     
      cs_free(x);
      return (ok ? C : ABIP(cs_spfree)(C));          
}

cs *ABIP(cs_compress)
(
      const cs *T
) 
{
      abip_int m; 
      abip_int n;  
      abip_int nnz;  
      abip_int p; 
      abip_int k; 

      abip_int *Cp;  
      abip_int *Ci; 
      abip_int *w; 
      abip_int *Ti; 
      abip_int *Tj;

      abip_float *Cx; 
      abip_float *Tx;
      
      cs *C;
      
      m = T->m;
      n = T->n;
      Ti = T->i;
      Tj = T->p;
      Tx = T->x;
      nnz = T->nnz;

      C = ABIP(cs_spalloc)(m, n, nnz, Tx != ABIP_NULL, 0);
      w = (abip_int *) abip_calloc(n, sizeof(abip_int));

      if (!C || !w) 
      {
            return (cs_done(C, w, ABIP_NULL, 0));
      }

      Cp = C->p;
      Ci = C->i;
      Cx = C->x;

      for (k = 0; k < nnz; k++)
      {
            w[Tj[k]]++;                                              
      }

      ABIP(cs_cumsum)(Cp, w, n);                         
      
      for (k = 0; k < nnz; k++) 
      {
            Ci[p = w[Tj[k]]++] = Ti[k];                  
            if (Cx) 
            {
                  Cx[p] = Tx[k];
            }
      }

      return (cs_done(C, w, ABIP_NULL, 1));           
}

cs *ABIP(cs_spalloc)
(
      abip_int m, 
      abip_int n, 
      abip_int nnzmax, 
      abip_int values,
      abip_int triplet
) 
{
      cs *A = (cs *) abip_calloc(1, sizeof(cs));        
      if (!A) 
      {
              return (ABIP_NULL);
      }                                                                       
      
      A->m = m;                                                       
      A->n = n;
      
      A->nnzmax = nnzmax = MAX(nnzmax, 1);
      A->nnz = triplet ? 0 : -1;                                
      
      A->p = (abip_int *) abip_malloc((triplet ? nnzmax : n + 1) * sizeof(abip_int));
      A->i = (abip_int *) abip_malloc(nnzmax * sizeof(abip_int));
      A->x = values ? (abip_float *) abip_malloc(nnzmax * sizeof(abip_float)) : ABIP_NULL;
      
      return ((!A->p || !A->i || (values && !A->x)) ? ABIP(cs_spfree)(A) : A);
}

cs *ABIP(cs_spfree)
(
      cs *A
) 
{
      if (!A) 
      {
            return (ABIP_NULL);
      } 
                                                                   
      cs_free(A->p);
      cs_free(A->i);
      cs_free(A->x);

      return ((cs *)cs_free(A));                      
}

abip_float ABIP(cs_cumsum)
(
      abip_int *p, 
      abip_int *c, 
      abip_int n
) 
{
      abip_int i; 
      abip_int nnz = 0;
      abip_float nnz2 = 0;
      
      if (!p || !c) 
      {
            return (-1);
      }
      
      for (i = 0; i < n; i++) 
      {
            p[i] = nnz;
            nnz += c[i];
            nnz2 += c[i];                                                                    
            c[i] = p[i];                                                                         
      }
      
      p[n] = nnz;
      return (nnz2);                                                                         
}

abip_int *ABIP(cs_pinv)
(
      abip_int const *p, 
      abip_int n
) 
{
      abip_int k; 
      abip_int *pinv;
      
      if (!p) 
      {
            return (ABIP_NULL);
      }                                                                                               

      pinv = (abip_int *)abip_malloc(n * sizeof(abip_int));            
      if (!pinv) 
      {
            return (ABIP_NULL);
      }
      
      for (k = 0; k < n; k++)
      { 
            pinv[p[k]] = k;
      }     
           
      return (pinv); 
}

cs *ABIP(cs_symperm)
(
      const cs *A, 
      const abip_int *pinv, 
      abip_int values
) 
{
      abip_int i;  
      abip_int j; 
      abip_int p; 
      abip_int q;  
      abip_int i2;  
      abip_int j2; 
      abip_int n; 
      abip_int *Ap; 
      abip_int *Ai; 
      abip_int *Cp; 
      abip_int *Ci; 
      abip_int *w;

      abip_float *Cx; 
      abip_float *Ax;
      
      cs *C;
      
      n = A->n;
      Ap = A->p;
      Ai = A->i;
      Ax = A->x;
      
      C = ABIP(cs_spalloc)(n, n, Ap[n], values && (Ax != ABIP_NULL), 0);                
      w = (abip_int *) abip_calloc(n, sizeof(abip_int));                                                  
      
      if (!C || !w) 
      {
            return (cs_done(C, w, ABIP_NULL, 0));
      }
      
      Cp = C->p;
      Ci = C->i;
      Cx = C->x;
      
      for (j = 0; j < n; j++)                                       
      {
            j2 = pinv ? pinv[j] : j;                                 
            
            for (p = Ap[j]; p < Ap[j + 1]; p++) 
            {
                   i = Ai[p];
                   if (i > j) 
                   {
                         continue;
                   } 
                   i2 = pinv ? pinv[i] : i;                                                          
                   w[MAX(i2, j2)]++;
            }
      }
      
      ABIP(cs_cumsum)(Cp, w, n);                              
      
      for (j = 0; j < n; j++) 
      {
            j2 = pinv ? pinv[j] : j;                                  
            
            for (p = Ap[j]; p < Ap[j + 1]; p++) 
            {
                   i = Ai[p];
                   
                   if (i > j) 
                   {
                         continue;
                   }                                                             
                   
                   i2 = pinv ? pinv[i] : i;                            
                   Ci[q = w[MAX(i2, j2)]++] = MIN(i2, j2);
                   
                   if (Cx) 
                   {
                          Cx[q] = Ax[p];
                   }
            }
      }
            
      return (cs_done(C, w, ABIP_NULL, 1));             
}
