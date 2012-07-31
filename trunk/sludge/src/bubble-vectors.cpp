#include "bubble.hpp"



void Bubble:: compute_area()
{
    const Tracer *p = root;
    Vertex  v0(0,0);
    Real ans= 0;
    for( size_t i=size;i>0;--i,p=p->next )
    {
        const Vertex v1 = v0 + p->edge;
        ans += v0.x * v1.y - v0.y * v1.x;
        v0 = v1;
    }
    area = Fabs(ans)/2;
}

void Bubble:: compute_geometry()
{
    //--------------------------------------------------------------------------
    // first pass: frenet (t,n) for each point
    //--------------------------------------------------------------------------
    Tracer *p = root;
    for( size_t i=size;i>0;--i,p=p->next )
    {
        p->compute_frenet();
    }
    
    //--------------------------------------------------------------------------
    // second pass: curvature
    //--------------------------------------------------------------------------
    for( size_t i=size;i>0;--i,p=p->next )
    {
        p->compute_curvature();
    }
}