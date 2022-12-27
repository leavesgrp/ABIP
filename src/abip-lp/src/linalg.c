#include "linalg.h"
#include <math.h>

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
    return SQRTF(ABIP(norm_sq)(v, len));
}

abip_float ABIP(min_abs_sqrt)
(
 const abip_float *a,
 abip_int len,
 abip_float ref
 )
{
    abip_int i;
    abip_float tmp;
    for (i = 0; i < len; ++i)
    {
        tmp = ABS(a[i]);
        if (tmp <= ref && tmp > 0)
        {
            ref = tmp;
        }
    }
    return SQRTF(ref);
}

// add by Kurt
abip_float ABIP(norm_one)
(
 const abip_float *v,
 abip_int len
 )
{    
    abip_int i;
    abip_float nmone = 0.0;
    for (i = 0; i < len; ++i)
    {
        nmone += ABS(v[i]);
    }
    return nmone;
}
abip_float ABIP(norm_one_sqrt)
(
 const abip_float *v,
 abip_int len
 )
{
    return SQRTF(ABIP(norm_one)(v, len));
}
// ------



/* max(|v|) */
abip_float ABIP(norm_inf)
(
 const abip_float *a,
 abip_int len
 )
{
    abip_int i;
    abip_float tmp;
    abip_float max = 0.0;
    for (i = 0; i < len; ++i)
    {
        tmp = ABS(a[i]);
        if (tmp >= max)
        {
            max = tmp;
        }
    }
    return max;
}

// add by Kurt
abip_float ABIP(norm_inf_sqrt)
(
 const abip_float *v,
 abip_int len
 )
{
    return SQRTF(ABIP(norm_inf)(v, len));
}
// ------

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
