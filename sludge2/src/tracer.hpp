#ifndef TRACER_INCLUDED
#define TRACER_INCLUDED 1

#include "types.hpp"
#include "yocto/sequence/cached-list.hpp"
#include "yocto/core/clist.hpp"

class Tracer
{
public:
    Tracer() throw();
    ~Tracer() throw();

    Tracer *next;
    Tracer *prev;
    
    Vertex  vertex; //!< position (should be PBC)
    Vertex  edge;   //!< vector to next->vertex
    
    typedef cache_of<Tracer>                     Cache;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Tracer);
};

typedef cached_list< core::clist_of, Tracer> Tracers; //!< circular list


#endif
