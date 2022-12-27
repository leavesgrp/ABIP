/* ========================================================================= */
/* === AMD_order =========================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: DrTimothyAldenDavis@gmail.com                                      */
/* ------------------------------------------------------------------------- */

/* User-callable AMD minimum degree ordering routine.  See amd.h for
 * documentation.
 */

#include "amd_internal.h"

/* ========================================================================= */
/* === AMD_order =========================================================== */
/* ========================================================================= */

GLOBAL Int AMD_order
(
            Int n,
            const Int Ap [ ],
            const Int Ai [ ],
            Int P [ ],
            abip_float Control [ ],
            abip_float ABIPInfo [ ]
)
{
            Int *Len; 
            Int *S;
            Int nz;
            Int i;
            Int *Pinv;
            Int info; 
            Int status; 

            Int *Rp; 
            Int *Ri; 
            Int *Cp; 
            Int *Ci; 
            Int ok;
            
            size_t nzaat; 
            size_t slen;
            
            abip_float mem = 0 ;

            #ifndef NDEBUG
            AMD_debug_init ("amd") ;
            #endif

            /* clear the ABIPInfo array, if it exists */
            info = ABIPInfo != (abip_float *) ABIP_NULL ;
            
            if (info)
            {
	           for (i = 0 ; i < AMD_INFO ; i++)
	           {
	                       ABIPInfo [i] = EMPTY ;
	           }
	           ABIPInfo [AMD_N] = n ;
	           ABIPInfo [AMD_STATUS] = AMD_OK ;
            }

            /* make sure inputs exist and n is >= 0 */
            if (Ai == (Int *) ABIP_NULL || Ap == (Int *) ABIP_NULL || P == (Int *) ABIP_NULL || n < 0)
            {
	           if (info) ABIPInfo [AMD_STATUS] = AMD_INVALID ;
	           return (AMD_INVALID) ;	    /* arguments are invalid */
            }

            if (n == 0)
            {
	           return (AMD_OK) ;	    /* n is 0 so there's nothing to do */
            }

            nz = Ap [n] ;
            
            if (info)
            {
	           ABIPInfo [AMD_NZ] = nz ;
            }
            
            if (nz < 0)
            {
	           if (info) ABIPInfo [AMD_STATUS] = AMD_INVALID ;
	           return (AMD_INVALID) ;
            }

            /* check if n or nz will cause size_t overflow */
            if (((size_t) n) >= SIZE_T_MAX / sizeof (Int) || ((size_t) nz) >= SIZE_T_MAX / sizeof (Int))
            {
	           if (info) ABIPInfo [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	           return (AMD_OUT_OF_MEMORY) ;	    /* problem too large */
            }

            /* check the input matrix:	AMD_OK, AMD_INVALID, or AMD_OK_BUT_JUMBLED */
            status = AMD_valid (n, n, Ap, Ai) ;

            if (status == AMD_INVALID)
            {
	           if (info) ABIPInfo [AMD_STATUS] = AMD_INVALID ;
	           return (AMD_INVALID) ;	    /* matrix is invalid */
            }

            /* allocate two size-n integer workspaces */
            Len  = amd_malloc (n * sizeof (Int)) ;
            Pinv = amd_malloc (n * sizeof (Int)) ;
            mem += n ;
            mem += n ;
            
            if (!Len || !Pinv)
            {
	           /* :: out of memory :: */
	           amd_free (Len) ;
	           amd_free (Pinv) ;
	           if (info) ABIPInfo [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	           return (AMD_OUT_OF_MEMORY) ;
            }

            if (status == AMD_OK_BUT_JUMBLED)
            {
	           /* sort the input matrix and remove duplicate entries */
	           AMD_DEBUG1 (("Matrix is jumbled\n")) ;
	           Rp = amd_malloc ((n+1) * sizeof (Int));
	           Ri = amd_malloc (nz * sizeof (Int));
	           mem += (n+1) ;
	           mem += MAX (nz,1) ;
	           
                        if (!Rp || !Ri)
	           {
	                       /* :: out of memory :: */
	                       amd_free (Rp) ;
	                       amd_free (Ri) ;
	                       amd_free (Len) ;
	                       amd_free (Pinv) ;
	    
                                    if (info) ABIPInfo [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	                       return (AMD_OUT_OF_MEMORY) ;
	           }

	           /* use Len and Pinv as workspace to create R = A' */
	           AMD_preprocess (n, Ap, Ai, Rp, Ri, Len, Pinv) ;
	           Cp = Rp ;
	           Ci = Ri ;
            }
            else
            {
	           /* order the input matrix as-is.  No need to compute R = A' first */
	           Rp = ABIP_NULL ;
	           Ri = ABIP_NULL ;
	           Cp = (Int *) Ap ;
	           Ci = (Int *) Ai ;
            }

            /* --------------------------------------------------------------------- */
            /* determine the symmetry and count off-diagonal nonzeros in A+A' */
            /* --------------------------------------------------------------------- */

            nzaat = AMD_aat (n, Cp, Ci, Len, P, ABIPInfo) ;
            AMD_DEBUG1 (("nzaat: %g\n", (abip_float) nzaat)) ;
            ASSERT ((MAX (nz-n, 0) <= nzaat) && (nzaat <= 2 * (size_t) nz)) ;

            /* --------------------------------------------------------------------- */
            /* allocate workspace for matrix, elbow room, and 6 size-n vectors */
            /* --------------------------------------------------------------------- */

            S = ABIP_NULL ;
            slen = nzaat ;			             /* space for matrix */
            ok = ((slen + nzaat/5) >= slen) ;	/* check for size_t overflow */
            slen += nzaat/5 ;			/* add elbow room */
    
            for (i = 0 ; ok && i < 7 ; i++)
            {
	           ok = ((slen + n) > slen) ;	/* check for size_t overflow */
	           slen += n ;			/* size-n elbow room, 6 size-n work */
            }
            
            mem += slen ;
            ok = ok && (slen < SIZE_T_MAX / sizeof (Int)) ;     /* check for overflow */
            ok = ok && (slen < Int_MAX) ;	                                /* S[i] for Int i must be OK */
            
            if (ok)
            {
	           S = amd_malloc (slen * sizeof (Int)) ;
            }
            
            AMD_DEBUG1 (("slen %g\n", (abip_float) slen)) ;
            
            if (!S)
            {
	           /* :: out of memory :: (or problem too large) */
	           amd_free (Rp) ;
	           amd_free (Ri) ;
	           amd_free (Len) ;
	           amd_free (Pinv) ;
	           if (info) ABIPInfo [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	           return (AMD_OUT_OF_MEMORY) ;
            }
            
            if (info)
            {
	           /* memory usage, in bytes. */
	           ABIPInfo [AMD_MEMORY] = mem * sizeof (Int) ;
            }

            /* --------------------------------------------------------------------- */
            /* order the matrix */
            /* --------------------------------------------------------------------- */

            AMD_1 (n, Cp, Ci, P, Pinv, Len, slen, S, Control, ABIPInfo) ;

            /* --------------------------------------------------------------------- */
            /* free the workspace */
            /* --------------------------------------------------------------------- */

            amd_free (Rp) ;
            amd_free (Ri) ;
            amd_free (Len) ;
            amd_free (Pinv) ;
            amd_free (S) ;
            
            if (info) ABIPInfo [AMD_STATUS] = status ;
            return (status) ;	    /* successful ordering */
}
