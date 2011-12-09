#include "types.h"

real_t icp_alea()
{
	static const real_t fac = 1.0/(1.0+(real_t)RAND_MAX);
	return (0.5+(real_t)rand()) * fac;
}

indx_t  icp_rand_index( const indx_t imin, const indx_t imax )
{
	assert(imin<=imax);
	{
		const indx_t n = imax - imin + 1;
		const real_t h = icp_alea();
		const indx_t r = ( (indx_t)(h*n) ) % n;
		return imin+r;
	}
}
