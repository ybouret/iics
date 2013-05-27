#include "bubble.hpp"
#include "yocto/code/fourcc.hpp"

const int Bubble::Tag = int(YOCTO_FOURCC('B', 'u', 'b', 'l'));

Bubble:: Bubble( Real &lam ) throw() :
prev(0),
next(0),
lambda( lam ),
G()
{
    
}

Bubble:: ~Bubble() throw()
{
}

void Bubble::hash_bubble(Hasher &h) const throw()
{
    hash_ring(h);
    h(G);
}
