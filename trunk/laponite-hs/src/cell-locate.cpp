#include "cell.hpp"

static unit_t __locate_point(Real           x, 
                             const Array1D &a,
                             unit_t         ilo,
                             unit_t         ihi)
{
    if( x <= a[ilo] )
        return ilo;
    
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
    //==========================================================================
    // first pass: locate in which grid we are
    //==========================================================================
    
    //--------------------------------------------------------------------------
    // by simple bissection
    //--------------------------------------------------------------------------
    const V2D    P = p.vertex;
    assert(P.x>=0);
    assert(P.x<=Length.x);
    assert(P.y>=-Length.y/2);
    assert(P.y<=Length.y/2);
    const unit_t i = p.pos.x = __locate_point(P.x, X, lower.x, upper.x);
    const unit_t j = p.pos.y = __locate_point(P.y, Y, lower.y, upper.y);
    
    
    //fprintf( stderr, " (%g,%g) <= (%g,%g) <= (%g,%g)\n", X[i],Y[j], v.x, v.y, X[i+1], Y[j+1]);
    
    //--------------------------------------------------------------------------
    // compute coefficient of bilinear interpolations
    //--------------------------------------------------------------------------
    p.w.x = (P.x - X[i])/dX[i];
    p.w.y = (P.y - Y[j])/dY[j];
    
    //==========================================================================
    // second pass: use p_next to determine possible intersections
    //==========================================================================
    const V2D Q = P + p.r_next;
    
}


//! locate all points in all spots
void Cell:: locate_points( )
{
    B.ldz();
    for( size_t i=segments.size();i>0;--i) 
        segments[i].empty();
    
    for( Bubble *b = bubbles.first(); b; b=b->next )
    {
        for( Spot *s = b->spots.head; s; s=s->next )
        {
            locate_point( * (s->point) );
        }
    }
    
   

}

