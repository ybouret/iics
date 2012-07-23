#include "segment.hpp"

Segment::~Segment() throw() {}

Segment:: Segment() throw() :
handle(0),
next(0),
prev(0)
{
}

void Segment:: reset() throw()
{
    handle = 0;
}

