
#ifndef GMARKER_INCLUDED
#define GMARKER_INCLUDED 1

#include "./types.hpp"
#include "yocto/core/cached-list.hpp"
#include "yocto/core/list.hpp"

class GridMarker
{
public:
    U2D      pos;
    GridMarker *next;
    GridMarker *prev;
    
    GridMarker() throw();
    ~GridMarker() throw();
    void reset() throw();
    
    typedef cache_of<GridMarker>                  Pool;
    typedef cached_list<core::list_of,GridMarker> List;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(GridMarker);
};


#endif
