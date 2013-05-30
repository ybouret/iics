#include "bubble.hpp"
#include "yocto/code/fourcc.hpp"

const int Bubble::Tag = 2;

Bubble:: Bubble( Real &lam ) throw() :
prev(0),
next(0),
lambda( lam ),
G(),
area(0),
pressure(1),
flags(0)
{
    
}

Bubble:: ~Bubble() throw()
{
}

void Bubble::hash_bubble(Hasher &h) const throw()
{
    hash_ring(h);
    h(G);
    h(area);
    h(pressure);
}


Tracer * Bubble:: append()
{
    Tracer *tr = new Tracer();
    push_back(tr);
    return tr;
}