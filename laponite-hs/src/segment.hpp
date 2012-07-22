#ifndef SEGMENT_INCLUDED
#define SEGMENT_INCLUDED 1

#include "intersection.hpp"

class Segment
{
public:
    Segment() throw();
    ~Segment() throw();
    void reset() throw();
    Intersection *inter;
    Segment      *next;
    Segment      *prev;
    Segment( const Segment &other ) throw();   
    
    typedef cache_of<Segment>                  Pool;
    typedef cached_list<core::list_of,Segment> List;
private:
    YOCTO_DISABLE_ASSIGN(Segment);
};

#endif


