#include "bubbles.hpp"


Bubbles:: Bubbles( const V2D &box_dim ) throw() :
Length(box_dim),
lambda(1.0),
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

Bubble * Bubbles:: append()
{
    const BubbleID id = BubbleID(b_list.size);
    Bubble *pB = 0;
    if( b_pool.size > 0 )
    {
        pB = b_pool.query();
        (BubbleID &)(pB->id) = id;
    }
    else
    {
        pB = new  Bubble( id, Length, (Real&)lambda, pcache, scache );
    }
    
   
    
    assert(pB->size == 0 );
    b_list.push_back(pB);
    return pB;
    
}


void Bubbles:: append( size_t n )
{
    while(n-->0) (void) append();
}

size_t Bubbles:: count() const throw()
{
    return b_list.size;
}


Bubble * Bubbles:: first() throw() { return b_list.head; }
const Bubble * Bubbles:: first() const throw() { return b_list.head; }

void Bubbles:: update_topologies()
{
    for( Bubble *b = b_list.head; b; b=b->next )
    {
        b->update_contour();
    }
}

#if 0
void Bubbles:: update_properties()
{
    for( Bubble *b = b_list.head; b; b=b->next )
    {
        b->compute_values();
    }
}
#endif

void Bubbles:: spots_and_values_within( const Real y_lo, const Real y_up )
{
    for( Bubble *b = b_list.head; b; b=b->next )
    {
        b->mark_and_find_spots_within(y_lo,y_up);
        if(b->active)
        {
            b->compute_values();
        }
    }
}

