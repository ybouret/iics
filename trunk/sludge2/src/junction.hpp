#ifndef JUNCTION_INCLUDED
#define JUNCTION_INCLUDED 1

#include "bubble.hpp"
#include "yocto/core/list.hpp"

//! bubble/mesh junctions
class Junction
{
public:
    Junction() throw();
    ~Junction() throw();
    
    Junction       *next;
    Junction       *prev;
    Vertex          pos;
    unit_t          klo;
    unit_t          khi;
    const Bubble   *bubble;
    Real            curvature;
    
    typedef cache_of<Junction> Cache;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junction);
};

typedef cached_list<core::list_of, Junction> Junctions;

#endif

