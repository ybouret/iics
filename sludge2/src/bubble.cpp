#include "bubble.hpp"

Bubble:: ~Bubble() throw()
{
}

Bubble:: Bubble(BubbleID      bubble_id,
                const PBC    &bubble_pbc,
                Real          &bubble_lam,
                Tracer::Cache &tcache,
                Spot::Cache   &scache ) throw() :
Tracers(tcache),
next(0),
prev(0),
id(bubble_id),
pbc(bubble_pbc),
lam(bubble_lam),
spots(scache)
{
    
}

void Bubble:: clear() throw()
{
    empty();
    spots.empty();
}

void  Bubble:: hash( hashing::function &h ) const
{
    const Tracer *tracer = root;
    h.run( &size, sizeof(size) );
    for( size_t i=size;i>0;--i,tracer=tracer->next)
        tracer->hash(h);
}

size_t Bubble:: get_hash( hashing::function &h) const
{
    h.set();
    hash(h);
    return h.key<size_t>();
}
