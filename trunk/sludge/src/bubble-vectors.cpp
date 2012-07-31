#include "bubble.hpp"

void Bubble:: raw_initialize()
{
    
    //--------------------------------------------------------------------------
    // first pass: pbc, edge, s
    //--------------------------------------------------------------------------
    Tracer *p = root;
    for( size_t i=size;i>0;--i,p=p->next )
    {
        pbc(p->vertex);
        const Tracer *q = p->next; assert(q!=NULL);
        Vertex        pq(p->vertex,q->vertex);
        pbc(pq);
        p->edge = pq;
        p->s2   = p->edge.norm2();
        p->s    = Sqrt( p->s2 );
    }
    content = 0;
    compute_geometry();
    
}


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
    // first pass: area and frenet (t,n) for each point
    //--------------------------------------------------------------------------
    Tracer *p = root;
    Vertex  v0(0,0);
    area = 0;
    for( size_t i=size;i>0;--i,p=p->next )
    {
        p->compute_frenet();
        const Vertex v1 = v0 + p->edge;
        area += v0.x * v1.y - v0.y * v1.x;
        v0 = v1;
    }
    area     = Fabs(area)/2;
    pressure = content / area;
    //--------------------------------------------------------------------------
    // second pass: curvature
    //--------------------------------------------------------------------------
    for( size_t i=size;i>0;--i,p=p->next )
    {
        p->compute_curvature();
    }
}