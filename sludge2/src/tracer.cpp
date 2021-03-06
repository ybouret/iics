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
curvature(0),
pressure(0),
gt(0),
bubble(0),
is_spot(false),
jnext(0),
jprev(0),
spd(0),
dsc(0)
{}


void Tracer:: hash( hashing::function &h ) const
{
    h.run( &vertex,IO_COUNT*sizeof(Real) );
}


#if 0
void Tracer:: compute_gt()
{
    assert(next);
    assert(prev);
    const Real sp2 = s2;
    const Real sm2 = prev->s2;
    const Real Pp  = next->pressure - pressure;
    const Real Pm  = prev->pressure - pressure;
    
    gt = ((sm2 * Pp - sp2 * Pm) / dsc) / spd;
}
#endif