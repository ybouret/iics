#include "spot.hpp"
#include <cstring>


////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Spot:: Pool:: Pool() throw() : CorePool()
{
}

Spot:: Pool:: ~Pool() throw()
{
    while( size ) delete query();
}


////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Spot:: List:: List( Spot::Pool &cache ) throw() : CoreList(), cache_( cache )
{
    
}

Spot:: List:: ~List() throw()
{
    empty();
}

void Spot:: List:: empty() throw()
{
    while( size ) cache_.store( pop_back() );
}

void Spot::List:: append( Point *p )
{
    assert(p!=NULL);
    Spot *spot = cache_.size > 0 ? cache_.query() : new Spot;
    memset(spot,0,sizeof(Spot));
    spot->point = p;
    push_back(spot);
}
