#include "fish.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

Point::  Point() throw() : i(Index++), r() {}
Point:: ~Point() throw() {}

size_t Point::Index = 0;

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Slice:: ~Slice() throw() {}

Slice:: Slice(double zz,double pr) throw() :
z(zz),
perimeter(pr),
points()
{
}


////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Triangle:: ~Triangle() throw() {}


Triangle:: Triangle(const pPoint &A,
                    const pPoint &B,
                    const pPoint &C) throw():
a(A),b(B),c(C),
G( (1.0/3) * (a->r+b->r+c->r) ),
n(),
S(0)
{
    assert(a->i!=b->i);
    assert(a->i!=c->i);
    assert(b->i!=c->i);
    {
        //const vtx_t G = (1.0/3) * (a->r+b->r+c->r);
        const vtx_t AB(a->r,b->r);
        const vtx_t AC(a->r,c->r);
        const vtx_t NN = vtx_t::cross_(AB, AC);
        const vtx_t middle(0,0,0.5);
        const vtx_t Q(middle,G);
        if( (NN*Q) <= 0 )
        {
            b.swap_with(c);
        }
    }

    {
        const vtx_t AB(a->r,b->r);
        const vtx_t AC(a->r,c->r);
        n = vtx_t::cross_(AB, AC);
        S = n.norm();
        n.normalize();
    }
}


Triangle:: Triangle(const Triangle &other) throw() :
a(other.a), b(other.b), c(other.c),
G(other.G),
n(other.n),
S(other.S)
{
}

void Triangle:: recompute() throw()
{
    G = (1.0/3) * (a->r+b->r+c->r);
    const vtx_t AB(a->r,b->r);
    const vtx_t AC(a->r,c->r);
    n = vtx_t::cross_(AB, AC);
    S = n.norm();
    n.normalize();
}


