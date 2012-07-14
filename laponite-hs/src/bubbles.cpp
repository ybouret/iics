#include "bubbles.hpp"


Bubbles:: Bubbles() throw() :
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

Bubble & Bubbles:: create()
{
    
    Bubble *pB = b_pool.size > 0 ? b_pool.query() : new  Bubble( pcache, scache );
    
    assert(pB->size == 0 );
    b_list.push_back(pB);
    return *pB;
    
}

size_t Bubbles:: count() const throw()
{
    return b_list.size;
}

Bubble * Bubbles:: first() throw() { return b_list.head; }
const Bubble * Bubbles:: first() const throw() { return b_list.head; }



