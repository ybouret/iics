#ifndef  TYPES_INCLUDED
#define  TYPES_INCLUDED 1

#include "yocto/geom/v2d.hpp"
#include "yocto/math/types.hpp"

using namespace yocto;
using namespace math;

typedef double          Real;
typedef geom::v2d<Real> V2D;


void AleaInit() throw();
Real Alea() throw();



#endif
