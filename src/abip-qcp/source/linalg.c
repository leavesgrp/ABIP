#include "linalg.h"
#include <math.h>

/* x = x .* y*/
void ABIP(c_dot)
(
      abip_float *x,
      const abip_float *y,
      const abip_int len
)
{
      for(int i=0; i<len; i++){
            x[i] *= y[i];
      }
}

/* y = mean(x) */
abip_float ABIP(vec_mean)
(
      abip_float *x,
      abip_int len
)
{
      if(len<=0 || x == ABIP_NULL){
            printf("invalid ABIP(vec_mean) parameter");
            return -1;
      }
      abip_float y=0;

      for(int i=0;i<len;i++){
            y += x[i];
      }

      return y/len;
}

/* x = b*a */
void ABIP(set_as_scaled_array)
(
      abip_float *x, 
      const abip_float *a, 
      const abip_float b,
      abip_int len
) 
{
      abip_int i;
      for (i = 0; i < len; ++i)
      {
            x[i] = b * a[i];
      }
}

/* x = sqrt(v) */
void ABIP(set_as_sqrt)
(
      abip_float *x, 
      const abip_float *v, 
      abip_int len
)
{
      abip_int i; 
      for (i = 0; i < len; ++i)
      {
            x[i] = SQRTF(v[i]); 
      }
}

/* x = v.^2 */
void ABIP(set_as_sq)
(
      abip_float *x, 
      const abip_float *v, 
      abip_int len
)
{
      abip_int i; 
      for (i = 0; i < len; ++i)
      {
            x[i] = v[i]*v[i]; 
      }
}  

/* a *= b */
void ABIP(scale_array)
(
      abip_float *a, 
      const abip_float b, 
      abip_int len
) 
{
      if(a == ABIP_NULL){
            return;
      }
      abip_int i;
      for (i = 0; i < len; ++i)
      {
            a[i] *= b;
      }
}

/* x'*y */
abip_float ABIP(dot)
(
      const abip_float *x, 
      const abip_float *y, 
      abip_int len
) 
{

      if(x == ABIP_NULL || y == ABIP_NULL){
            return 0;
      }

      abip_int i;
      abip_float ip = 0.0;
      for (i = 0; i < len; ++i) 
      {
            ip += x[i] * y[i];
      }
      return ip;
}

/* ||v||_2^2 */
abip_float ABIP(norm_sq)
(
      const abip_float *v, 
      abip_int len
) 
{    
      if(v == ABIP_NULL){
            return 0;
      }

      abip_int i;
      abip_float nmsq = 0.0;
      for (i = 0; i < len; ++i) 
      {
            nmsq += v[i] * v[i];
      }
      return nmsq;
}

/* ||v||_2 */
abip_float ABIP(norm)
(
      const abip_float *v, 
      abip_int len
) 
{
      if(v == ABIP_NULL){
            return 0;
      }
      return SQRTF(ABIP(norm_sq)(v, len));
}

/* ||x||_1 */
abip_float ABIP(norm_1)
(
	const abip_float *x, 
	const abip_int len
)
{
      if(x == ABIP_NULL){
            return 0;
      }
      abip_float result = 0;
      for(int i=0; i<len; i++){
            result += ABS(x[i]);
      }
      return result;
}

/*the absolute value of the largest component of x*/
abip_float ABIP(cone_norm_1)
(
	const abip_float *x,
	const abip_int len
)
{
    abip_int i;
    abip_float tmp; 
    abip_float max = 0.0;
    for (i = 0; i < len; ++i) 
    {
        tmp = x[i];
        if (tmp > max) 
        {
                max = tmp;
        }
    }
    return ABS(max);
}

/* max(|v|) */
abip_float ABIP(norm_inf)
(
      const abip_float *a, 
      abip_int len
) 
{

      if(a == ABIP_NULL || len == 0){
            return 0;
      }
      abip_int i;
      abip_float tmp; 
      abip_float max = 0.0;
      for (i = 0; i < len; ++i) 
      {
            tmp = ABS(a[i]);
            if (tmp > max) 
            {
                  max = tmp;
            }
      }
      return max;
}

/* a .+= b */
void ABIP(add_array)
(
      abip_float *a, 
      const abip_float b, 
      abip_int len
)
{
      abip_int i; 
      for (i = 0; i <len; ++i) 
      {
            a[i] += b;
      }
} 

/* saxpy a += sc*b */
void ABIP(add_scaled_array)
( 
      abip_float *a, 
      const abip_float *b, 
      abip_int len,
      const abip_float sc
) 
{
      if(b == ABIP_NULL){
            return;
      }
      abip_int i;
      for (i = 0; i < len; ++i) 
      {
            a[i] += sc * b[i];
      }
}

/* ||a-b||_2^2 */
abip_float ABIP(norm_diff)
(
      const abip_float *a, 
      const abip_float *b, 
      abip_int len
) 
{
      abip_int i;
      abip_float tmp;
      abip_float nm_diff = 0.0; 
      for (i = 0; i < len; ++i) 
      {
            tmp = (a[i] - b[i]);
            nm_diff += tmp * tmp;
      }
      return SQRTF(nm_diff);
}

/* max(|a-b|) */
abip_float ABIP(norm_inf_diff)
(
      const abip_float *a, 
      const abip_float *b,
      abip_int len
) 
{
      abip_int i;
      abip_float tmp; 
      abip_float max = 0.0;
      for (i = 0; i < len; ++i) 
      {
            tmp = ABS(a[i] - b[i]);
            if (tmp > max) 
            {
                  max = tmp;
            }
      }
      return max;
}

abip_int arr_ind(const abip_int i_col, const abip_int i_row, const abip_int nrows, const abip_int ncols, const abip_int format) {
	return (format == RowMajor) ? (i_col + i_row * ncols) : (i_row + i_col * nrows);
}

abip_float * ABIP(csc_to_dense)(const cs * in_csc, const abip_int out_format) {
	abip_int i_row, i_col, nnz_in_col, i_val = 0, i_nnz;
	const abip_int nrows = in_csc->m;
	const abip_int ncols = in_csc->n;
	abip_float * out_matrix = (abip_float *)abip_malloc(nrows * ncols * sizeof(abip_float));
      memset(out_matrix,0,nrows * ncols * sizeof(abip_float));
	abip_int * col_nnz = in_csc->p;
	abip_int * rows = in_csc->i;
	abip_float * values = in_csc->x;

	for (i_col = 0; i_col < ncols; i_col++) {
		nnz_in_col = col_nnz[i_col + 1] - col_nnz[i_col];
		if (nnz_in_col > 0) {
			for (i_nnz = 0; i_nnz < nnz_in_col; i_nnz++) {
				i_row = rows[i_val];
				out_matrix[arr_ind(i_col, i_row, nrows, ncols, out_format)] = values[i_val];
				i_val++;
			}
		}
	}

	return out_matrix;
}