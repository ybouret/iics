#include "bubble.hpp"


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

void Bubble:: append( const Vertex v)
{
    append()->pos = v;
}

void Bubble:: collect_markers( const Real ymin, const Real ymax)
{
    markers.clear();
    Tracer *tracer = root;
    size_t        old   = 0;
    
    for( size_t i=0; i < size; ++i, tracer=tracer->next )
    {
        const Vertex &v = tracer->pos;
        if(v.y>=ymin && v.y<=ymax)
        {
            markers.append(tracer, i-old);
            old=i;
        }
    }
    
    //std::cerr << "Bubble #" << UID << " has " << markers.size << " markers" << std::endl;
    
    
}
