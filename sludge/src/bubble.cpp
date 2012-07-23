#include "bubble.hpp"

Bubble:: Bubble(Real                &lambda_ref, 
                const PBC           &pbc_ref,
                Tracer::Cache       &tracer_cache,
                Tracer::HandleCache &tracer_handle_cache  
                ) throw() :
Tracer::List( tracer_cache ),
lambda(lambda_ref),
pbc(pbc_ref),
area(0),
handles( tracer_handle_cache )
{
    
}

Bubble:: ~Bubble() throw()
{
    
}