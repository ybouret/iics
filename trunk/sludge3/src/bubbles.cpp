#include "bubbles.hpp"
#include "yocto/code/fourcc.hpp"


const int Bubbles::Tag = 3;

Bubbles:: ~Bubbles() throw() { auto_delete(); }


Bubbles:: Bubbles( Real &lam ) throw() :
lambda(lam)
{
}

Bubble * Bubbles:: append()
{
    Bubble *b = new Bubble((Real&)lambda,size);
    push_back(b);
    return b;
}

void Bubbles:: hash( Hasher &h ) const throw()
{
    h(size);
    for( const Bubble *b = head;b;b=b->next)
    {
        b->hash_bubble(h);
    }
}
