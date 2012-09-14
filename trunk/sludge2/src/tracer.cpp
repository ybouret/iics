#include "tracer.hpp"

Tracer:: ~Tracer() throw() {}

Tracer:: Tracer() throw() :
next(0),prev(0),
vertex(), edge() {}


void Tracer:: hash( hashing::function &h ) const
{
    h.run( &vertex,IO_COUNT*sizeof(Real) );
}