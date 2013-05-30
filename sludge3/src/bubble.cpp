#include "bubble.hpp"
#include "yocto/code/fourcc.hpp"

////////////////////////////////////////////////////////////////////////////////
//
// Marker
//
////////////////////////////////////////////////////////////////////////////////
Marker:: Marker( const Tracer *tr, const size_t s) :
tracer(tr),
shift(s)
{
    assert(tracer);
}

Marker:: ~Marker() throw()
{
}

////////////////////////////////////////////////////////////////////////////////
//
// Marker List
//
////////////////////////////////////////////////////////////////////////////////

Marker:: List:: List() throw() {}

Marker:: List:: ~List() throw()
{
    clear();
}

void Marker::List:: clear() throw()
{
    auto_delete();
}

void Marker:: List:: append( const Tracer *tracer, const size_t shift)
{
    Marker *m = new Marker(tracer,shift);
    push_back(m);
}

////////////////////////////////////////////////////////////////////////////////
//
// Bubble
//
////////////////////////////////////////////////////////////////////////////////
const int Bubble::Tag = 2;

Bubble:: Bubble( Real &lam, Real &gam, size_t uid) throw() :
prev(0),
next(0),
lambda( lam ),
gamma( gam ),
G(),
area(0),
pressure(1),
flags(0),
UID(uid)
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
    h(UID);
}


Tracer * Bubble:: append()
{
    Tracer *tr = new Tracer();
    push_back(tr);
    return tr;
}

void Bubble:: collect_markers( const Real ymin, const Real ymax)
{
    markers.clear();
    const Tracer *tracer = root;
    size_t        old   = 0;
    
    for( size_t i=0; i < size; ++i )
    {
        const Vertex &v = tracer->pos;
        if(v.y>=ymin && v.y<=ymax)
        {
            markers.append(tracer, i-old);
            old=i;
        }
    }
    
    
}
