#include "marker.hpp"

////////////////////////////////////////////////////////////////////////////////
//
// Marker
//
////////////////////////////////////////////////////////////////////////////////
Marker:: Marker( Tracer *tr, const size_t s) :
tracer(tr),
shift(s),
gt(0),
gn(0)
{
    assert(tracer);
}

Marker:: ~Marker() throw()
{
}

#if 0
void Marker:: find_anchor( const Array &B)
{
    assert(tracer);
    
    const Vertex ini = tracer->pos;
    //const Vertex dir = - tracer->n;
    //std::cerr << "Anchor for " << ini << " towards " << dir << std::endl;
    
    assert(0==tracer->flags);         //!< inside
    assert(B.has(tracer->coord));     //!< logical lower coordinates
}
#endif

////////////////////////////////////////////////////////////////////////////////
//
// Marker List
//
////////////////////////////////////////////////////////////////////////////////

Marker:: List:: List() throw() {}

Marker:: List:: ~List() throw()
{
    clear();
}

void Marker::List:: clear() throw()
{
    auto_delete();
}

void Marker:: List:: append( Tracer *tracer, const size_t shift)
{
    Marker *m = new Marker(tracer,shift);
    push_back(m);
}
