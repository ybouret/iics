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
Real Anint( Real u ) throw();

//! return the Periodic Boundary Condition value of y
/**
 -L/2 <= y <= L/2
 */
Real PBC1( Real y, const Real L, const Real invL ) throw();

class PBC
{
public:
    const Real L;
    const Real invL;
    PBC( Real length ) throw();
    ~PBC() throw();
    PBC( const PBC &other ) throw();
    
    Real operator()( Real y ) throw(); //!< make -L/2 <= y <= L/2
    void operator()( V2D &v ) throw(); //!< act on y
    
private:
    YOCTO_DISABLE_ASSIGN(PBC);
};



#endif
