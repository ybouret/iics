#include "types.hpp"


#include "yocto/code/rand.hpp"

void AleaInit() throw()
{
    _rand.wseed();
}

Real Alea() throw()
{
    return alea<Real>();
}

Real Anint( Real x ) throw()
{
    static const Real __half = 0.5;
    return Floor( x + __half );
}

Real PBC1( Real x, const Real L, const Real invL ) throw()
{
     return x - L * Anint( invL * x );
}

PBC:: PBC( Real length ) throw() : L(length), invL(1/L)
{
}

PBC::~PBC() throw() {}

PBC:: PBC( const PBC &other ) throw() : L(other.L), invL(other.invL) 
{
}

Real PBC:: operator()( Real y ) throw()
{
    return PBC1(y,L,invL);
}

void PBC:: operator()( V2D &v ) throw()
{
    v.y = (*this)(v.y);
}
