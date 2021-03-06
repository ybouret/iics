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

PBC:: PBC( Real length ) throw() : L(length), invL(1/L), lo(-L/2), up(L/2)
{
}

PBC::~PBC() throw() {}

PBC:: PBC( const PBC &other ) throw() : L(other.L), invL(other.invL), lo(other.lo), up(other.up)
{
}

Real PBC:: operator()( Real y ) const throw()
{
    const Real ans = PBC1(y,L,invL);
    assert( ans >= lo);
    assert( ans <  up);
    return ans;
}

void PBC:: operator()( V2D &v ) const throw()
{
    v.y = (*this)(v.y);
}
