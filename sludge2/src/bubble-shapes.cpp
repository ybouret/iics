#include "bubble.hpp"
#include "yocto/code/utils.hpp"

void Bubble:: map_circle( const Vertex &center, Real radius )
{
    assert(lam>0);
    empty();
    const Real   theta_max = 2 * atan( lam/(radius+radius) );
    const size_t nmin      = max_of<size_t>(3,size_t( ceil( numeric<Real>::two_pi/theta_max) ));
    const Real   dtheta = numeric<Real>::two_pi / nmin;
    const Real   theta0 = numeric<Real>::two_pi * Alea();
    for( size_t i=0; i < nmin; ++i )
    {
        const Real  theta = i * dtheta + theta0;
        Tracer     *p     = append();
        p->vertex.x = center.x + radius * Cos( theta );
        p->vertex.y = center.y + radius * Sin( theta );
    }

}

