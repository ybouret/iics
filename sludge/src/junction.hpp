#ifndef JUNCTION_INCLUDED
#define JUNCTION_INCLUDED 1

#include "types.hpp"
#include "yocto/core/cached-list.hpp"
#include "yocto/core/list.hpp"

class Bubble;

class Junction 
{
public:
    Junction() throw();
    ~Junction() throw();
    
    void reset() throw();
    
    Vertex    vertex;
    Bubble   *bubble;
    Junction *next;
    Junction *prev;
    
    typedef cache_of<Junction>                  Cache;
    typedef cached_list<core::list_of,Junction> List;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junction);
};

#endif
