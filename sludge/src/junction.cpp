#include "junction.hpp"
#include "yocto/core/offset.hpp"
#include <cstring>

Junction:: ~Junction() throw()
{
    
}

Junction:: Junction() throw() :
vertex(),
bubble(0),
lo(0),
up(0),
next(0),
prev(0)
{
    
}


void Junction:: reset() throw()
{
    YOCTO_HARD_RESET(Junction,vertex,this);
}


void Junction:: copy( const Junction *J ) throw()
{
    assert(J!=NULL);
    vertex = J->vertex;
    bubble = J->bubble;
    lo     = J->lo;
    up     = J->up;
}