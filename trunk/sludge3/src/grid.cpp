#include "grid.hpp"
#include "yocto/ios/ocstream.hpp"

void __Grid:: SaveDat( const Grid &grid, const string &fn )
{
    ios::ocstream fp(fn,false);
    const Array1D &X = grid.X();
    const Array1D &Y = grid.Y();
    
    // gnuplot format
    for( unit_t i=X.lower; i <= X.upper; ++i )
    {
        const Real x = X[i];
        for( unit_t j=Y.lower; j <= Y.upper; ++j )
        {
            fp("%g %g\n",x, Y[j]);
        }
        fp("\n");
    }
    
    for( unit_t j=Y.lower; j <= Y.upper; ++j )
    {
        const Real y = Y[j];
        for( unit_t i=X.lower; i <= X.upper; ++i )
        {
            fp("%g %g\n",X[i],y);
        }
        fp("\n");
    }
    
}
