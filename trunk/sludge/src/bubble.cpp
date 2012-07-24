#include "bubble.hpp"

Bubble:: Bubble(Real                &lambda_ref, 
                const PBC           &pbc_ref,
                Tracer::Cache       &tracer_cache,
                Spot::Cache         &spot_cache,
                Marker::Cache       &mcache
                ) throw() :
Tracer::List( tracer_cache ),
id(0),
lambda(lambda_ref),
pbc(pbc_ref),
area(0),
spots( spot_cache ),
markers(mcache),
borders(mcache),
active(false),
next(0),prev(0)
{
    
}

Bubble:: ~Bubble() throw()
{
    
}

void Bubble:: clear() throw()
{
    empty();
    spots.empty();
    markers.empty();
    borders.empty();
}

void Bubble:: collect_spots_within(const Real y_lo, const Real y_up)
{
    spots.empty();
    Tracer *p = root;
    size_t last_index = 0;
    for( size_t i=0;i<size;++i,p=p->next)
    {
        p->bubble = this;
        const double y = p->vertex.y;
        if( y_lo <= y && y < y_up )
        {
            spots.attach(p);           
            spots.tail->jump = i-last_index;
            last_index = i;
            p->is_spot = true;
        }
        else 
        {
            p->is_spot = false;
        }
    }
    active = spots.size > 0;
}
