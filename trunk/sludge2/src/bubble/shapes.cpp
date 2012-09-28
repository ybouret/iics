#include "../bubble.hpp"
#include "yocto/code/utils.hpp"

void Bubble:: rotate( Real alpha )
{
    Vertex G;
    Tracer *p  = root;
    Vertex  v0 = p->vertex;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        G  += v0;
        v0 += p->edge;
    }
    G /= Real(size);
    p = root;
    const Real c = Cos(alpha);
    const Real s = Sin(alpha);
    for( size_t i=size;i>0;--i,p=p->next)
    {
        Vertex r(G,p->vertex);
        pbc(r);
        const Vertex q(r.x * c - r.y*s,r.x*s+r.y*c);
        p->vertex = G + q;
        pbc(p->vertex);
    }
    
}

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
    update_area();
}

void Bubble:: map_peanut( const Vertex &center, Real radius, Real alpha )
{
    assert(lam>0);
    empty();
    alpha = clamp<Real>(0,Fabs(alpha),1);
    const Real b  = Fabs(radius);
    const Real a  = alpha * b;
    const Real a2 = a*a;
    const Real b2 = b*b;
    const Real a4 = a2*a2;
    const Real b4 = b2*b2;
    const Real   theta_max = 2 * atan( lam/(radius+radius) );
    const size_t nmin      = 2*max_of<size_t>(3,size_t( ceil( numeric<Real>::two_pi/theta_max) ));
    const Real   dtheta    = numeric<Real>::two_pi / nmin;
    
    //fprintf( stderr, "peanut @(%g,%g), b=%g, a=%g\n", center.x, center.y, b, a );
    
    for( size_t i=0; i < nmin; ++i )
    {
        const Real   theta = i * dtheta;;
        const Real   t2    = theta+theta;
        Tracer      *p     = append();
        const Real   S     = Sin(t2);
        const  Real  c1    = max_of<Real>(0,b4 - a4 * S*S );
        const  Real  C     = max_of<Real>(0,a2 * Cos( t2 ) + Sqrt(c1));
        const  Real  rho   = Sqrt(C);
        p->vertex.x = center.x + rho * Cos( theta );
        p->vertex.y = center.y + rho * Sin( theta );
    }
    update_area();
}
