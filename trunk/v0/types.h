#ifndef ICP_TYPES_INCLUDED
#define ICP_TYPES_INCLUDED 1
#ifdef __cplusplus
extern "C" {
#endif
    
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <assert.h>
#include <time.h>

typedef double      real_t;
#define ICP_REAL    MPI_DOUBLE

typedef ptrdiff_t   indx_t;

/*! return a value in ]0:1[ */
real_t  icp_alea();

/*! return a value in [imin:imax] */
indx_t  icp_rand_index( const indx_t imin, const indx_t imax );
#ifdef __cplusplus
    }
#endif
        
    
#endif

