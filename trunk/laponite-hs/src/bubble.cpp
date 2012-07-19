#include "bubble.hpp"
#include "yocto/code/utils.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

Bubble:: Bubble( BubbleID who, Real L, Point::Pool &pcache, Spot::Pool &scache ) throw() : 
Point::List( pcache ),
id(who),
pbc(L),
lambda(1),
area(0),
spots( scache ),
active( true ),
next(0),
prev(0)
{
}

Bubble::~Bubble() throw()
{
}

void Bubble:: mark_and_find_spots_within(const Real y_lo, const Real y_up)
{
    spots.empty();
    Point *p = root;
    size_t last_index = 0;
    for( size_t i=0;i<size;++i,p=p->next)
    {
        p->bubble = this;
        const double y = p->vertex.y;
        if( y_lo <= y && y <= y_up )
        {
            spots.append(p);           
            spots.tail->jump = i-last_index;
            last_index = i;
        }
    }
    active = spots.size > 0 ;
}

