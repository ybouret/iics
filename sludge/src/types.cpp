#include "./types.hpp"
#include "yocto/code/rand.hpp"

void AleaInit() throw()
{
    _rand.wseed();
}

Real Alea() throw()
{
    return alea<Real>();
}

PBC:: PBC( Real length ) throw() : L(length), invL(1/L), lo(-L/2), up(L/2)
{
}

PBC::~PBC() throw() {}

PBC:: PBC( const PBC &other ) throw() : L(other.L), invL(other.invL), lo(other.lo), up(other.up)
{
}


Real PBC:: apply( Real y ) const throw()
{
    return PBC1(y, L, invL);
}

void PBC:: operator()(Vertex &v) const throw()
{
    v.y = apply(v.y);
}