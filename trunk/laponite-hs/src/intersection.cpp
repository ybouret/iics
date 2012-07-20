#include "intersection.hpp"

Intersection:: Intersection() throw() :
vertex(),
lo(0),
up(0),
next(0),
prev(0)
{
    
}

Intersection:: ~Intersection() throw()
{
}



void Intersection:: reset() throw()
{
    vertex.x = vertex.y = 0;
    lo = up  = 0;
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
