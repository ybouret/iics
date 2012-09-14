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

