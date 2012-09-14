#include "bubbles.hpp"

Bubbles:: ~Bubbles() throw()
{
    while( pool.size )    delete pool.query();
    while( bubbles.size ) delete bubbles.pop_back();
}

Bubbles:: Bubbles( const PBC &bubbles_pbc ) throw() :
pbc( bubbles_pbc ),
lambda(1),
bubbles(),
pool(),
tcache(),
scache()
{
}

Bubble * Bubbles:: first() throw()
{
    return bubbles.head;
}

const Bubble * Bubbles:: first() const throw()
{
    return bubbles.head;

}

size_t Bubbles:: count() const throw()
{
    return bubbles.size;
}

Bubble * Bubbles:: append()
{
    const BubbleID id = bubbles.size+1;
    Bubble        *b  = 0;
    if(pool.size > 0)
    {
        b = pool.query();
        (BubbleID&)(b->id) = id;
    }
    else
    {
        b = new Bubble( id, pbc, lambda, tcache, scache);
    }
    bubbles.push_back(b);
    return b;
}


void Bubbles:: clear() throw()
{
    while( bubbles.size )
    {
        Bubble *b = bubbles.pop_back();
        b->empty();
        pool.store( b );
    }
}

void Bubbles:: create( size_t n )
{
    clear();
    for(size_t i=n;i>0;--i)
    {
        (void) append();
    }
}
