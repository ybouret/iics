#include "tracer.hpp"

const int Tracer::Tag = 1;

Tracer:: Tracer() throw() :
prev(0),next(0),
pos(),
edge(),
dist(0),
t(),
n(),
C(0)
{}


Tracer:: Tracer( const Vertex v ) throw() :
prev(0),
next(0),
pos(v),
edge(),
dist()
{
    
}

Tracer:: ~Tracer() throw() {}

void Tracer:: hash_tracer( Hasher &h ) const throw()
{
    h(pos);
    h(edge);
    h(dist);
    h(t);
    h(n);
    h(C);
}

Tracer:: Ring:: Ring() throw() {}

Tracer:: Ring:: ~Ring() throw() { auto_delete(); }

void Tracer:: Ring:: hash_ring( Hasher &h ) const throw()
{
    h(size);
    const Tracer *tr = root;
    for(size_t i=size;i>0;--i,tr=tr->next)
    {
        tr->hash_tracer(h);
    }
}


void Tracer:: compute_curvature()
{
    assert(prev);
    assert(next);
    assert(dist>0);
    assert(prev->dist>0);
    assert(next->dist>0);
    const Vertex M0Mm = -prev->edge;
    const Vertex M0Mp =  edge;
    const Real   tm   = prev->dist;
    const Real   tp   = next->dist;
    const Real   h    = tm+tp;
    const Vertex dM    = (tm/tp * M0Mp - tp/tm * M0Mm) / h;
    const Real   speed = dM.norm();
    t = dM / speed;
    
    n.x = -t.y;
    n.y =  t.x;
    
    const Vertex d2M = (M0Mp/tp + M0Mm/tm) / h;
    
    C = Vertex::det(dM,d2M) / (speed*speed*speed);
}

