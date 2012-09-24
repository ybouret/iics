#ifndef MARKER_INCLUDED
#define MARKER_INCLUDED 1

#include "bubble.hpp"

class Marker
{
public:
    Marker() throw();
    ~Marker() throw();
    Marker *next;
    Marker *prev;
    Coord         inside; //!< inside coordinate
    const Bubble *bubble; //!< belonging to
    
    typedef cache_of<Marker> Cache;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Marker);
};

typedef cached_list< core::list_of, Marker > Markers;

#endif

