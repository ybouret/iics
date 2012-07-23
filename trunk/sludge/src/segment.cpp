#include "segment.hpp"

Segment::~Segment() throw() {}

Segment:: Segment() throw() :
handle(0),
next(0),
prev(0)
{
}

Segment:: Segment( const Segment &other ) throw() :
handle( other.handle ),
next(0), prev(0)
{
}


void Segment:: reset() throw()
{
    handle = 0;
}

