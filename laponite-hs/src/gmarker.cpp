#include "./gmarker.hpp"
#include "yocto/core/offset.hpp"
#include <cstring>
GridMarker:: ~GridMarker() throw()
{
}

GridMarker:: GridMarker() throw() :
pos(),
next(0), 
prev(0)
{
    
}

void GridMarker:: reset() throw()
{
    memset( &pos, 0, sizeof(GridMarker) - YOCTO_OFFSET_OF(GridMarker,pos));
}