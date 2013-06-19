#include "bubble.hpp"



void Bubble:: init_contour() throw()
{
    assert(size>=3);
    G.ldz();
    area = 0;
    Tracer *tr = root;
    for( size_t i=size;i>0;--i,tr=tr->next)
    {
        const Vertex &p = tr->pos;
        G += p;
        const Vertex &q = tr->next->pos;
        tr->edge = q - p;
        tr->dist = tr->edge.norm();
        area += p.x * q.y - p.y * q.x;
    }
    G.x /= size;
    G.y /= size;
    area = Fabs(area)/2;
}



void Bubble:: compute_curvatures()
{
    assert(size>=3);
    Tracer *tr = root;
    for(size_t i=size;i>0;--i,tr=tr->next)
    {
        tr->compute_order1();
    }
    
    assert(root==tr);
    for(size_t i=size;i>0;--i,tr=tr->next)
    {
        tr->compute_order2();
    }

}


void Bubble:: regularize()
{
    adjust_contour();
    compute_curvatures();
}
