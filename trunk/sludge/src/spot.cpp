#include "spot.hpp"

Spot:: Spot() throw() :
handle(0),
jump(0),
next(0),
prev(0)
{
}

Spot:: ~Spot() throw() 
{
}

void Spot:: reset() throw()
{
    handle = 0;
    jump   = 0;
}