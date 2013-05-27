#include "tracer.hpp"
#include "yocto/code/fourcc.hpp"

const int Tracer::Tag = int(YOCTO_FOURCC('T','R','A','C'));

Tracer:: Tracer() throw() :
prev(0),next(0),
pos(),
edge(),
dist(0)
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