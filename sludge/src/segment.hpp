#ifndef SEGMENT_INCLUDED
#define SEGMENT_INCLUDED 1

#include "junction.hpp"
#include "yocto/core/handle-list.hpp"


class Segment 
{
public:
    enum Kind
    {
        IsInternal,
        IsUpperPBC,
        IsLowerPBC
    };
    
    Segment() throw();
    ~Segment() throw();
    Segment( const Segment &other ) throw();
    
    Junction *handle;
    Segment  *next;
    Segment  *prev;
    
    void reset() throw();
    
    typedef cache_of<Segment>    Cache;
    typedef handle_list<Segment> List;
    
private:
    YOCTO_DISABLE_ASSIGN(Segment);
};

#endif
