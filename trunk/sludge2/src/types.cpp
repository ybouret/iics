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

PBC:: PBC( Real length ) throw() :
L(length),
invL(1/L),
lo(-L/2),
up(L/2)
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

#include "yocto/ios/ocstream.hpp"

void SaveGrid( const Grid &G, const string &filename )
{
    ios::ocstream fp( filename, false );
    const Array1D &X = G.X();
    const Array1D &Y = G.Y();
    
    for( unit_t j= G.lower.y; j<= G.upper.y;++j)
    {
        if( (j&1) )
        {
            for( unit_t i= G.lower.x; i <= G.upper.x; ++i )
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
        else
        {
            for( unit_t i= G.upper.x; i >= G.lower.x; --i )
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
    }
    
    for( unit_t i=G.lower.x; i <= G.upper.x; ++i )
    {
        if( (i&1) )
        {
            for( unit_t j=G.lower.y; j<=G.upper.y;++j)
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
        else
        {
            for( unit_t j=G.upper.y; j>=G.lower.y;--j)
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
    }
    
}
