#ifndef BUBBLE_INCLUDED
#define BUBBLE_INCLUDED 1

#include "tracer.hpp"

class Bubble  : public Tracer::List
{
public:
    explicit Bubble(Real                &lambda_ref, 
                    const PBC           &pbc_ref,
                    Tracer::Cache       &tracer_cache,
                    Tracer::HandleCache &tracer_handle_cache  
                    ) throw();
    virtual ~Bubble() throw();
    
    const Real         &lambda;
    const PBC          &pbc;
    Real                area;
    Tracer::HandleList  handles;
    
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
};

#endif
