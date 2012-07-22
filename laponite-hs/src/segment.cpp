#include "segment.hpp"


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
