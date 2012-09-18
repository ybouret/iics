#ifndef JUNCTION_INCLUDED
#define JUNCTION_INCLUDED 1

#include "types.hpp"
#include "yocto/sequence/cached-list.hpp"
#include "yocto/core/list.hpp"

//! bubble/mesh junctions
class Junction
{
public:
    Junction( const Vertex &at ) throw();
    ~Junction() throw();
    
    Junction       *next;
    Junction       *prev;
    const Vertex    pos;
    
    typedef cache_of<Junction> Cache;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junction);
};

typedef cached_list<core::list_of, Junction> Junctions;

#endif

