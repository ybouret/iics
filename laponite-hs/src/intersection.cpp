#include "intersection.hpp"
#include "yocto/core/offset.hpp"
#include <cstring>
Intersection:: Intersection() throw() :
vertex(),
lo(0),
up(0),
bubble(0),
next(0),
prev(0)
{
    
}

Intersection:: ~Intersection() throw()
{
}



void Intersection:: reset() throw()
{
	memset( &vertex,0,sizeof(Intersection)-YOCTO_OFFSET_OF(Intersection,vertex));
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
