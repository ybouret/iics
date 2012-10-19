#include "tracer.hpp"

Tracer:: ~Tracer() throw() {}

Tracer:: Tracer() throw() :
next(0),prev(0),
vertex(),
edge(),
s2(0),
s(0),
t(),
n(),
angle(0),
curvature(0),
pressure(0),
bubble(0),
is_spot(false),
jnext(0),
jprev(0)
{}


void Tracer:: hash( hashing::function &h ) const
{
    h.run( &vertex,IO_COUNT*sizeof(Real) );
}
