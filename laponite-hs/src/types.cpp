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