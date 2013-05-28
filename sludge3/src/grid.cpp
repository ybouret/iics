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

Real __Grid:: ComputeLambda( const Grid &grid) throw()
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

static inline
unit_t __locate( const Array1D &A, const Real a ) throw()
{
    unit_t jlo = A.lower;
    unit_t jup = A.upper;
    assert(jup>jlo);
    assert(a>=A[jlo]);
    assert(a<=A[jup]);
    
    while(jup-jlo>1)
    {
        const unit_t mid = (jlo+jup)>>1;
        if( a < A[mid] )
        {
            jup = mid;
        }
        else
        {
            jlo = mid;
        }
    }
    
    return jlo;
}

int __Grid::Locate(const Grid &grid, const Vertex &p, Coord &lo) throw()
{
    assert( grid.is_valid() );
    const Array1D &X = grid.X();
    const Array1D &Y = grid.Y();
    int ans = SLUDGE_INSIDE;
    
    if( p.x < X[X.lower])
    {
        ans |= SLUDGE_LEFT;
    }
    else
    {
        if(p.x>X[X.upper])
        {
            ans |= SLUDGE_RIGHT;
        }
        else
        {
            //-- effective location
            lo.x = __locate( X, p.x );
        }
    }
    
    if( p.y < Y[Y.lower])
    {
        ans |= SLUDGE_BOTTOM;
    }
    else
    {
        if(p.y>Y[Y.upper])
        {
            ans |= SLUDGE_TOP;
        }
        else
        {
            //-- effective location
            lo.y = __locate( Y, p.y );
        }
    }

    
    
    
    return ans;
}
