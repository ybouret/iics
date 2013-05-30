#include "bubbles.hpp"
#include "yocto/code/fourcc.hpp"


const int Bubbles::Tag = 3;

Bubbles:: ~Bubbles() throw() { auto_delete(); }


Bubbles:: Bubbles() throw() :
lambda(1),
gamma(0)
{
}

Bubble * Bubbles:: append()
{
    Bubble *b = new Bubble(lambda,gamma,size);
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


void Bubbles:: regularize_all()
{
    for( Bubble *b = head;b;b=b->next)
    {
        b->regularize();
    }
}

void Bubbles:: collect_all_markers(Real ymin,Real ymax)
{
    for( Bubble *b = head;b;b=b->next)
    {
        b->collect_markers(ymin, ymax);
    }
}


