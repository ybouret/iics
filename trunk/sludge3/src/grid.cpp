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

#include "yocto/code/utils.hpp"

Real __Grid:: ComputeLambda( const Grid &grid)
{
    
    const Array1D &X = grid.X();
    assert(X.width>1);
    Real dx =Fabs( X[X.lower+1] - X[X.lower] );
    for( unit_t i=X.lower+1;i<X.upper;++i)
    {
        const Real tmp = Fabs(X[i+1]-X[i]);
        if( tmp < dx ) dx = tmp;
    }
    
    const Array1D &Y = grid.Y();
    Real dy = Fabs( Y[Y.lower+1] - Y[Y.lower] );
    for( unit_t j=Y.lower+1;j<Y.upper;++j)
    {
        const Real tmp = Fabs(Y[j+1] - Y[j]);
        if( tmp < dy ) dy = tmp;
    }
    
    return min_of(dx,dy)/2;
}
