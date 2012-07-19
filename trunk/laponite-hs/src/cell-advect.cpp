#include "cell.hpp"


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
