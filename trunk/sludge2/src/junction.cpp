#include "junction.hpp"

Junction:: ~Junction() throw() {}
Junction:: Junction() throw() :
next(0),
prev(0),
vertex(),
klo(0),
khi(0),
bubble(0),
alpha(0),
curvature(0),
t(),
n()
{}

#if 0
Junction::Junction( const Junction &other ) throw() :
next(0),
prev(0),
vertex( other.vertex ),
klo( other.klo ),
khi( other.khi ),
bubble( other.bubble ),
alpha( other.alpha ),
curvature(0),
n(other.n)
{
    
}
#endif
