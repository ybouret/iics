#include "cell.hpp"

static int __locate_point( Real x, const Array1D &a )
{
    int ilo = a.lower;
    if( x <= a[ilo] )
        return ilo;
    
    int ihi = a.upper;
    assert(ihi-ilo>=1);
    if( x >= a[ihi] )
        return ihi -1;
    
    while(ihi-ilo>1)
    {
        const int mid = (ilo+ihi)>>1;
        if( x < a[mid] )
        {
            ihi = mid;
        }
        else 
        {
            ilo = mid;
        }
    }
    return ilo;
}

void Cell:: locate_point( Point &p ) const
{
    const V2D v = p.vertex;
    p.i = __locate_point(v.x, X);
    p.j = __locate_point(v.y, Y);
    fprintf( stderr, " (%g,%g) <= (%g,%g) <= (%g,%g)\n", X[p.i],Y[p.j], v.x, v.y, X[p.i+1], Y[p.j+1]);
}
