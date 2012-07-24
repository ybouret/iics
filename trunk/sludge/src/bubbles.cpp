#include "bubbles.hpp"

size_t Bubbles:: count() const throw()
{
    return b_list.size;
}

Bubble * Bubbles::first() throw()
{
    return b_list.head;
}

const Bubble * Bubbles::first() const throw()
{
    return b_list.head;
}


Bubbles:: ~Bubbles() throw() 
{
    while( b_list.size ) delete b_list.pop_back();
    while( b_pool.size ) delete b_pool.query();
}

Bubbles:: Bubbles( const Vertex &box ) throw() :
pbc( box.y ),
tcache(),
scache(),
mcache(),
lambda(1),
b_list(),
b_pool()
{
}


void Bubbles:: empty() throw()
{
    while( b_list.size )
    {
        Bubble *bubble = b_list.pop_back();
        bubble->clear();
        b_pool.store(bubble);
    }
}

Bubble *Bubbles:: create()
{
    Bubble *bubble = 0;
    if( b_pool.size )
    {
        bubble = b_pool.query();
    }
    else 
    {
        bubble = new Bubble(lambda, pbc, tcache, scache, mcache);
    }
    
    b_list.push_back(bubble);
    
    bubble->id = b_list.size;
    
    return bubble;
}


void Bubbles:: fill( Array &B ) const
{
    B.ldz();
    for( const Bubble *bubble = first(); bubble; bubble=bubble->next )
    {
        const Real id = bubble->id; assert(id>0);
        for( const Marker *marker = bubble->markers.head; marker; marker=marker->next )
        {
            const Coord c = marker->coord;
            if( B.has(c) )
                B[c] = id;
        }
        
    }
}

