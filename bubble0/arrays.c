#include "arrays.h"

static size_t __bytes = 0;

size_t icp_bytes() { return __bytes; }

/*******************************************************************************
 ** 1D
 ******************************************************************************/
real_t *icp_create_array1D( indx_t ilo, indx_t ihi )
{
	const size_t n = ihi - ilo + 1;
	real_t      *a = NULL;
	assert(ihi>=ilo);
	a  = (real_t *)calloc( n, sizeof(real_t) );
	if( NULL == a )
		return NULL;
	a -= ilo;
	__bytes += n * sizeof(real_t);
	return a;
}

void    icp_delete_array1D( real_t *a, indx_t ilo, indx_t ihi )
{
	
	assert( ihi>=ilo  );
	if( a != NULL )
	{
		const size_t n  = ihi - ilo + 1;
		const size_t nb = n * sizeof(real_t);
		free( a + ilo );
		assert(nb<=__bytes);
		__bytes -= nb;
	}
	
}

/*******************************************************************************
 ** 2D
 ******************************************************************************/
real_t **icp_create_array2D( indx_t xmin, indx_t xmax, indx_t ymin, indx_t ymax )
{
	const size_t nx = xmax - xmin + 1;
	const size_t ny = ymax - ymin + 1;
	const size_t n  = nx * ny;
	size_t       i;
	real_t      **a = NULL;
	assert(xmax>=xmin);
	assert(ymax>=ymin);
	
	/** allocate rows **/
	a = (real_t **)calloc( ny, sizeof(real_t *));
	if( NULL == a )
	{
		return NULL;
	}
	
	/** allocate items **/
	a[0] = (real_t *)calloc(n,sizeof(real_t));
	if( NULL == a[0] )
	{
		free(a);
		return NULL;
	}
	
	/** register bytes **/
	__bytes += ny * sizeof(real_t*);
	__bytes += n  * sizeof(real_t);
	
	/** prepare rows pointers **/
	a[0] -= xmin;
	for(i=1;i<ny;++i) a[i] = a[i-1] + nx;
	
	/** adjust rows address **/
	a    -= ymin;
	
	return a;
}

void     icp_delete_array2D( real_t **a,  indx_t xmin, indx_t xmax, indx_t ymin, indx_t ymax )
{
	assert(xmax>=xmin);
	assert(ymax>=ymin);
	if( a != NULL )
	{
		const size_t nx = xmax - xmin + 1;
		const size_t ny = ymax - ymin + 1;
		const size_t n  = nx * ny;
		const size_t nb = n * sizeof(real_t) + ny * sizeof(real_t);
		assert(nb>=__bytes);
		a    += ymin;
		a[0] += xmin;
		free(a[0]);
		free(a);
		__bytes -= nb;
	}
	
}



/*******************************************************************************
 ** 3D
 ******************************************************************************/
real_t ***icp_create_array3D( indx_t xmin, indx_t xmax, indx_t ymin, indx_t ymax, indx_t zmin, indx_t zmax )
{
	const size_t nx = xmax - xmin + 1;
	const size_t ny = ymax - ymin + 1;
	const size_t nz = zmax - zmin + 1;
	const size_t nr = ny * nz; /* number of rows */
	const size_t n  = nx * nr;
	real_t    ***a  = NULL;
	size_t       j,k;
	
	assert(xmax>=xmin);
	assert(ymax>=ymin);
	assert(zmax>=zmin);
	
	/** allocate slices **/
	a = (real_t ***)calloc(nz,sizeof(real_t **));
	if( NULL == a )
	{
		return NULL;
	}
	
	/** allocate rows **/
	a[0] = (real_t **)calloc(nr,sizeof(real_t *));
	
	/** allocate items */
	a[0][0] = (real_t *)calloc(n,sizeof(real_t));
	
	/** link **/
	a[0][0] -= xmin;
	for( j=1; j < nr; ++j )
		a[0][j] = a[0][j-1] + nx;
	
	a[0] -= ymin;
	for( k=1; k < nz; ++k )
		a[k] = a[k-1] + ny;
	
	a -= zmin;
	
	__bytes += nz * sizeof(real_t **);
	__bytes += nr * sizeof(real_t *);
	__bytes += n  * sizeof(real_t);
	
	return a;
}

void icp_delete_array3D( real_t ***a, indx_t xmin, indx_t xmax, indx_t ymin, indx_t ymax, indx_t zmin, indx_t zmax )
{
	assert(xmax>=xmin);
	assert(ymax>=ymin);
	assert(zmax>=zmin);
	
	if( a != NULL )
	{
		const size_t nx = xmax - xmin + 1;
		const size_t ny = ymax - ymin + 1;
		const size_t nz = zmax - zmin + 1;
		const size_t nr = ny * nz; /* number of rows */
		const size_t n  = nx * nr;
		const size_t nb = nz * sizeof(real_t **) + nr * sizeof(real_t *) + n  * sizeof(real_t);
		assert( __bytes >= nb );
		a       += zmin;
		a[0]    += ymin;
		a[0][0] += xmin;
		free( a[0][0] );
		free( a[0] );
		free( a );
		__bytes -= nb;
	}
}