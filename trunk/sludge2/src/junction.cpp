#include "junction.hpp"

Junction:: ~Junction() throw() {}
Junction:: Junction() throw() :
next(0),
prev(0),
#if JUNCTION_TAG == 1
kind(0),
tag(0),
#endif
vertex(),
klo(0),
khi(0),
bubble(0),
alpha(0),
curvature(0),
pressure(0),
t(),
n(),
gt(0),
visited(false),
gn(0)
{}

#if 0
Real Junction:: Peff( const Vertex &pos ) const throw()
{
    const Vertex dr(vertex,pos);
    return pressure + dr * g;
}
#endif
