#ifndef TYPES_INCLUDED
#define TYPES_INCLUDED 1

#include "yocto/swamp/array2d.hpp"
#include "yocto/math/types.hpp"

using namespace yocto;
using namespace swamp;
using namespace math;

typedef double                Real;
typedef vertex2D<Real>::type  Vertex;
typedef coord2D               Coord;

void AleaInit() throw();
Real Alea() throw();

inline Real PBC1( Real x, const Real L, const Real invL ) throw()
{
    static const Real __half = 0.5;
    return x - L * Floor( (invL*x) + __half );
}


class PBC
{
public:
    const Real L;
    const Real invL;
    const Real lo;
    const Real up;
    PBC( Real length ) throw();
    ~PBC() throw();
    PBC( const PBC &other ) throw();
    
    Real apply( Real y ) const throw();
    
    void operator()( Vertex &v ) const throw(); //!< act on y
    
private:
    YOCTO_DISABLE_ASSIGN(PBC);
};


#endif
