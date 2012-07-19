#include "bubble.hpp"
#include "yocto/code/utils.hpp"

void Bubble:: map_circle(const V2D &center, Real radius)
{
    assert(lambda>0);
    empty();
    const Real   theta_max = 2 * atan( lambda/(radius+radius) );
    const size_t nmin      = max_of<size_t>(3,size_t( ceil( numeric<Real>::two_pi/theta_max) ));
    const Real   dtheta = numeric<Real>::two_pi / nmin;
    const Real   theta0 = numeric<Real>::two_pi * Alea();
    for( size_t i=0; i < nmin; ++i )
    {
        const Real  theta = i * dtheta + theta0;
        Point       *p     = create();
        p->vertex.x = center.x + radius * Cos( theta );
        p->vertex.y = center.y + radius * Sin( theta );
        push_back(p);
    }
}

#if 0
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
#endif

void Bubble:: map_peanut( const V2D &center, Real radius, Real alpha )
{
    assert(lambda>0);
    empty();
    alpha = clamp<Real>(0,Fabs(alpha),1);
    const Real b = Fabs(radius);
    const Real a = alpha * b;
    const Real a2 = a*a;
    const Real b2 = b*b;
    const Real a4 = a2*a2;
    const Real b4 = b2*b2;
    const Real   theta_max = 2 * atan( lambda/(radius+radius) );
    const size_t nmin      = max_of<size_t>(3,size_t( ceil( numeric<Real>::two_pi/theta_max) ));
    const Real   dtheta    = numeric<Real>::two_pi / nmin;
    const Real   theta0    = numeric<Real>::two_pi * Alea();
    
    for( size_t i=0; i < nmin; ++i )
    {
        const Real   theta = i * dtheta + theta0;
        const Real   t2    = theta+theta; 
        Point       *p     = create();
        const  Real  c1    = max_of<Real>(0,b4 - a4 * Sin( t2 ) );
        const  Real  C     = a2 * Cos( t2 ) + Sqrt(c1);
        const  Real  rho   = Sqrt(C);
        p->vertex.x = center.x + rho * Cos( theta );
        p->vertex.y = center.y + rho * Sin( theta );
        push_back(p);
    }

    
}