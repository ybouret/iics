#include "tracer.hpp"
#include "yocto/core/offset.hpp"
#include <cstring>

Tracer:: ~Tracer() throw()
{
    
}

Tracer:: Tracer() throw() :
vertex(),
edge(),
s(0),
t(),
n(),
next(0),
prev(0)
{
    
}

void Tracer:: set_normal() throw()
{
    n.x = -t.y;
    n.y =  t.x;
}

void Tracer:: reset() throw()
{
    YOCTO_HARD_RESET(Tracer,vertex,this);
}