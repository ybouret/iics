#include "tracer.hpp"

const int Tracer::Tag = 1;

#define _TRACER_CTOR(args) \
prev(0),\
next(0),\
pos(args),\
edge(),\
dist(0),\
t(),\
n(),\
C(0),\
speed(0),\
coord(),\
flags(0)

Tracer:: Tracer() throw() :
_TRACER_CTOR()
{}


Tracer:: Tracer( const Vertex &v ) throw() :
_TRACER_CTOR(v)
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


void Tracer:: compute_order1()
{
    assert(prev);
    assert(next);
    assert(dist>0);
    assert(prev->dist>0);
    assert(next->dist>0);
    const Real   tm    = prev->dist;
    const Real   tp    = next->dist;
    const Vertex Vm    = -prev->edge/tm;
    const Vertex Vp    =  edge/tp;
    const Real   h     = 1.0/(tm+tp);
    const Vertex dM    = h * (tm*Vp - tp * Vm);
    speed = dM.norm();
    t = dM / speed;
    
    n.x = -t.y;
    n.y =  t.x;
}

void Tracer:: compute_order2()
{
    assert(prev);
    assert(next);
    assert(dist>0);
    assert(prev->dist>0);
    assert(next->dist>0);
    const Real   tm    = prev->dist;
    const Real   tp    = next->dist;
    const Vertex Vm    = (prev->t -t)/tm;
    const Vertex Vp    = (next->t -t)/tp;
    const Real   h     = 1.0/(tm+tp);
    const Vertex tmp   = h * (tm*Vp-tp*Vm);
    
    C = (tmp*n)/speed;
    //std::cerr << "C=" << C << std::endl;
}


Real Tracer::Ring::  __area() const throw()
{
    Real ans = 0;
    Tracer *tr = root;
    for( size_t i=size;i>0;--i,tr=tr->next)
    {
        const Vertex &p = tr->pos;
        const Vertex &q = tr->next->pos;
        ans += p.x * q.y - p.y * q.x;
    }
    return Fabs(ans)/2;
}

