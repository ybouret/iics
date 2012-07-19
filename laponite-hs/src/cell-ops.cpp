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
    const V2D    v = p.vertex;
    assert(v.x>=0);
    assert(v.x<=Length.x);
    assert(v.y>=-Length.y/2);
    assert(v.y<=Length.y/2);
    const unit_t i = p.pos.x = __locate_point(v.x, X, lower.x, upper.x);
    const unit_t j = p.pos.y = __locate_point(v.y, Y, lower.y, upper.y);
   
    
    //fprintf( stderr, " (%g,%g) <= (%g,%g) <= (%g,%g)\n", X[i],Y[j], v.x, v.y, X[i+1], Y[j+1]);
    
    //--------------------------------------------------------------------------
    // compute coefficient of bilinear interpolations
    //--------------------------------------------------------------------------
    p.w.x = (v.x - X[i])/dX[i];
    p.w.y = (v.y - Y[j])/dY[j];
    
    //==========================================================================
    // second pass: find possible 4 intersections with grid
    //==========================================================================
    
    
}

void Cell:: advect_point( Point &p, double dt ) const
{
    const U2D q = p.pos;
    assert(q.x >= U.lower.x );
    assert(q.y >= U.lower.y );
    assert(q.x <  U.upper.x );
    assert(q.y <  U.upper.y );
    //--------------------------------------------------------------------------
    // collect velocities
    //--------------------------------------------------------------------------
    const V2D u00 = U[q.y  ][q.x  ];
    const V2D u01 = U[q.y+1][q.x  ];
    const V2D u11 = U[q.y+1][q.x+1];
    const V2D u10 = U[q.y  ][q.x+1];
    
    //--------------------------------------------------------------------------
    // prepare coefficients
    //--------------------------------------------------------------------------
    const double x   = p.w.x;
    const double y   = p.w.y;
    const double umx = 1-x;
    const double umy = 1-y;
    
    //--------------------------------------------------------------------------
    // create interpolated field
    //--------------------------------------------------------------------------
    const V2D    u   = umx*umy*u00 + x*y*u11 + umx*y * u01 + x*umy * u10;
    
    //--------------------------------------------------------------------------
    // move !
    //--------------------------------------------------------------------------
    p.vertex += dt * u;
    
    
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

void Cell:: advect_points( double dt )
{
    for( Bubble *b = bubbles.first(); b; b=b->next )
    {
        for( Spot *s = b->spots.head; s; s=s->next )
        {
            advect_point( * (s->point), dt );
        }
    }
}


void Cell:: collect_inside( vector<V2D> &pts ) const
{
    pts.free();
    for( unit_t j=lower.y; j <= upper.y; ++j )
    {
        for( unit_t i=lower.x; i <= upper.x; ++i )
        {
            if( B[j][i] > 0 )
            {
                const V2D v(X[i],Y[j]);
                pts.push_back( v );
            }
        }
    }
}
