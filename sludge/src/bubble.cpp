#include "bubble.hpp"

Bubble:: Bubble(Real                &lambda_ref, 
                const PBC           &pbc_ref,
                Tracer::Cache       &tracer_cache,
                Spot::Cache         &spot_cache 
                ) throw() :
Tracer::List( tracer_cache ),
lambda(lambda_ref),
pbc(pbc_ref),
area(0),
spots( spot_cache )
{
    
}

Bubble:: ~Bubble() throw()
{
    
}

void Bubble:: clear() throw()
{
    empty();
    spots.empty();
}