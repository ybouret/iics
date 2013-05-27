#include "bubble.hpp"

void Bubble:: init_contour() throw()
{
    assert(size>=3);
    Tracer *tr = root;
    for( size_t i=size;i>0;--i,tr=tr->next)
    {
        tr->edge = tr->next->pos - tr->pos;
        tr->dist = tr->edge.norm();
    }
}