#include "bubbles.hpp"


Bubbles:: Bubbles( Real box_height ) throw() :
Ly(box_height),
b_list(),
b_pool(),
pcache(),
scache()
{
}

Bubbles:: ~Bubbles() throw()
{
    while( b_list.size > 0 ) delete b_list.pop_back();
    while( b_pool.size > 0 ) delete b_pool.query();
}

void Bubbles:: none() throw()
{
    while( b_list.size > 0 )
    {
        Bubble *b = b_list.pop_back();
        b->empty();
        b_pool.store(b);
    }
}

Bubble * Bubbles:: create()
{
    
    Bubble *pB = b_pool.size > 0 ? b_pool.query() : new  Bubble( Ly, pcache, scache );
    
    assert(pB->size == 0 );
    b_list.push_back(pB);
    return pB;
    
}


void Bubbles:: create( size_t n )
{
    while(n-->0) (void) create();
}

size_t Bubbles:: count() const throw()
{
    return b_list.size;
}

Bubble * Bubbles:: first() throw() { return b_list.head; }
const Bubble * Bubbles:: first() const throw() { return b_list.head; }

void Bubbles:: update_all_points()
{
    for( Bubble *b = b_list.head; b; b=b->next )
    {
        b->update_points();
    }
}

void Bubbles:: update_all_values()
{
    for( Bubble *b = b_list.head; b; b=b->next )
    {
        b->update_values();
    }
}

void Bubbles:: update_all_spots( const Real y_lo, const Real y_up )
{
    for( Bubble *b = b_list.head; b; b=b->next )
    {
        b->build_spots(y_lo,y_up);
        if(b->active)
            b->update_values();
    }
}

