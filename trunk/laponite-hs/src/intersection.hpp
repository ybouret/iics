#ifndef INTERSECTION_INCLUDED
#define INTERSECTION_INCLUDED 1

#include "types.hpp"
#include "yocto/core/cached-list.hpp"
#include "yocto/core/list.hpp"

class Bubble;
class Intersection
{
public:
    Intersection() throw();
    ~Intersection() throw();
    V2D           vertex;
    ptrdiff_t     lo,up;
    Bubble       *bubble;
    Intersection *next;
    Intersection *prev;
    
    void reset() throw();
    
    typedef cache_of<Intersection>                   Pool;
    typedef cached_list<core::list_of, Intersection> List;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Intersection);
};




#endif
