#include "spot.hpp"
#include <cstring>

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
void Spot:: reset() throw()
{
    point = 0;
    jump  = 0;
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Spot:: List:: List( Spot::Pool &cache ) throw() : CoreList( cache )
{
    
}

Spot:: List:: ~List() throw()
{
}



void Spot::List:: attach( Point *p )
{
    assert(p!=NULL);
    append()->point = p;
}

Spot:: List:: List( const List &other ) :
CoreList( other.cache )
{
    for( const Spot *spot = other.head; spot; spot=spot->next )
    {
        attach( spot->point );
    }
}

