#ifndef TRACER_INCLUDED
#define TRACER_INCLUDED 1

#include "types.hpp"
#include "yocto/sequence/cached-list.hpp"
#include "yocto/core/clist.hpp"
#include "yocto/core/list.hpp"
#include "yocto/hashing/sha1.hpp"

class Tracer
{
public:
    Tracer() throw();
    ~Tracer() throw();

    Tracer *next;
    Tracer *prev;
    
    Vertex  vertex; //!< position (should be PBC) : +2
    Vertex  edge;   //!< vector to next->vertex   : +2
    Real    s2;     //!< |edge|^2                   +1
    
    static const size_t IO_COUNT = 5;
    
    typedef cache_of<Tracer>                     Cache;
    
    void hash( hashing::function &h ) const;
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Tracer);
};

typedef cached_list< core::clist_of, Tracer> Tracers; //!< circular list
typedef cached_list< core::list_of,  Tracer> SubTracers;

#endif
