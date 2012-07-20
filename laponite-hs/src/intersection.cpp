#include "intersection.hpp"

Intersection:: Intersection() throw() :
vertex(),
next(0),
prev(0)
{
    
}

Intersection:: ~Intersection() throw()
{
}

#if 0
Intersection:: Intersection( const Intersection &other ) throw() :
vertex(other.vertex),
next(0),
prev(0)
{
}
#endif

void Intersection:: reset() throw()
{
    vertex.x = vertex.y = 0;
}


Segment:: Segment() throw() :
inter(0),
next(0),
prev(0)
{
}

Segment:: ~Segment() throw()
{
}

Segment:: Segment( const Segment &other ) throw() :
inter( other.inter ),
next(0),
prev(0)
{
}

void Segment:: reset() throw()
{
    inter = 0;
}
