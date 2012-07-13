#ifndef  TYPES_INCLUDED
#define  TYPES_INCLUDED 1

#include "yocto/geom/v2d.hpp"
#include "yocto/math/types.hpp"

using namespace yocto;
using namespace math;

typedef double          Real;
#define MPI_REAL_TYPE   MPI_DOUBLE

typedef geom::v2d<Real> V2D;

//! initialize random number
void AleaInit() throw();

//! in ]0:1[
Real Alea() throw();

//! return the nearest integer
Real Anint( Real x ) throw();

//! return the Periodic Boundary Condition value of x
/**
 -L/2 <= x <= L/2
 */
Real PBC1( Real x, const Real L, const Real invL ) throw();



#endif
