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
lambda(1)
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
    Bubble *bubble = NULL;
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


