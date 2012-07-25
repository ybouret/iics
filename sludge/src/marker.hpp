#ifndef MARKER_INCLUDED
#define MARKER_INCLUDED 1

#include "types.hpp"
#include "yocto/core/cached-list.hpp"
#include "yocto/core/list.hpp"


//! coordinate if a point inside a bubble
class Marker 
{
public:
    Marker() throw();
    ~Marker() throw();
    
    Coord   coord;
    Marker *next;
    Marker *prev;
    void reset() throw();
    
    typedef cache_of<Marker>                  Cache;
    typedef cached_list<core::list_of,Marker> List;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Marker);
};

#endif
