#ifndef QDLDL_TYPES_H
# define QDLDL_TYPES_H

# ifdef __cplusplus
extern "C" {
# endif /* ifdef __cplusplus */

#include <limits.h> //for the QDLDL_INT_TYPE_MAX
#include"glbopts.h"

// QDLDL integer and float types

typedef abip_int   QDLDL_int;   /* for indices */
typedef abip_float  QDLDL_float; /* for numerical values  */
typedef abip_int   QDLDL_bool;  /* for boolean values  */

//Maximum value of the signed type QDLDL_int.
#define QDLDL_INT_MAX INT_MAX

# ifdef __cplusplus
}
# endif /* ifdef __cplusplus */

#endif /* ifndef QDLDL_TYPES_H */