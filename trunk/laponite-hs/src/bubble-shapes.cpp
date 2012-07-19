#include "bubble.hpp"
#include "yocto/code/utils.hpp"

void Bubble:: map_circle(const V2D &center, Real radius)
{
    assert(lambda>0);
    empty();
    const double theta_max = 2 * atan( lambda/(radius+radius) );
    const size_t nmin      = max_of<size_t>(3,size_t( ceil( numeric<Real>::two_pi/theta_max) ));
    const double dtheta = numeric<Real>::two_pi / nmin;
    const double theta0 = numeric<Real>::two_pi * Alea();
    for( size_t i=0; i < nmin; ++i )
    {
        const double theta = i * dtheta + theta0;
        Point       *p     = create();
        p->vertex.x = center.x + radius * Cos( theta );
        p->vertex.y = center.y + radius * Sin( theta );
        push_back(p);
    }
}

void Bubble:: map_astroid( const V2D &center, Real radius )
{
    assert(lambda>0);
    empty();
    const size_t ntop = max_of<size_t>(1,0.71*radius/lambda);
    const double dtheta = numeric<Real>::pi / (2*ntop);
    const size_t n = 4*ntop;
    const double theta0 = numeric<Real>::two_pi * Alea();
    for( size_t i=0; i < n; ++i )
    {
        const double theta = i * dtheta + theta0;
        Point       *p     = create();
        const Real C = Cos( theta );
        const Real S = Sin( theta );
        p->vertex.x = center.x + radius * C*C*C;
        p->vertex.y = center.y + radius * S*S*S;
        push_back(p);
    }
}
