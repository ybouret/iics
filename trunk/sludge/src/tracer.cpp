#include "tracer.hpp"
#include "yocto/core/offset.hpp"
#include <cstring>

Tracer:: ~Tracer() throw()
{
    
}

Tracer:: Tracer() throw() :
vertex(),
edge(),
s2(0),
s(0),
next(0),
prev(0),
t(),
n(),
curvature(0),
gLower(),
bw(),
bubble(0),
is_spot(false)
{
    
}

void Tracer:: set_normal() throw()
{
    n.x = -t.y;
    n.y =  t.x;
}

void Tracer:: reset() throw()
{
    YOCTO_HARD_RESET(Tracer,vertex,this);
}

void Tracer:: compute_frenet()
{
    assert(prev!=NULL);
    assert(next!=NULL);
    assert(s>0); assert(s2>0);
    assert(prev->s>0); assert(prev->s2>0);
    const Real    hm2  = prev->s2;
    const Real    hp2  = s2;
    const Vertex &dfp  = edge;
    const Vertex &dfm  = prev->edge;
    const Vertex t_dir = hm2 * dfp + hp2 * dfm;
    const Real   t_nrm = t_dir.norm(); // todo: check norm ?
    
    t = (1/t_nrm) * t_dir;
    
    set_normal();
    
}

#include "yocto/code/utils.hpp"

void Tracer:: compute_curvature()
{
    assert(prev!=NULL);
    assert(next!=NULL);
    assert(s>0); assert(s2>0);
    const Real    hm   = prev->s;
    const Real    hm2  = prev->s2;
    const Real    hp   = s;
    const Real    hp2  = s2;
    const Vertex  dfp(t,next->t);
    const Vertex  dfm(prev->t,t);
    const Vertex  dt  = hm2 * dfp + hp2 * dfm;
    const Real    den = hm2 * hp + hp2 * hm;
    curvature = (dt*n) / den;
}


