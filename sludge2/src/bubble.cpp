#include "bubble.hpp"

Bubble:: ~Bubble() throw()
{
}

Bubble:: Bubble(BubbleID      bubble_id,
                const PBC    &bubble_pbc,
                Real         &bubble_lam,
                Real         &bubble_gam,
                Tracer::Cache &tcache,
                Spot::Cache   &scache ) throw() :
Tracers(tcache),
next(0),
prev(0),
id(bubble_id),
pbc(bubble_pbc),
lam(bubble_lam),
gam(bubble_gam),
spots(scache),
pressure(0),
area(0),
content(0)
{
    
}

void Bubble:: clear() throw()
{
    empty();
    spots.empty();
}

void  Bubble:: hash( hashing::function &h ) const
{
    const Tracer *tracer = root;
    h.run( &size, sizeof(size) );
    for( size_t i=size;i>0;--i,tracer=tracer->next)
        tracer->hash(h);
}

size_t Bubble:: get_hash( hashing::function &h) const
{
    h.set();
    hash(h);
    return h.key<size_t>();
}

void Bubble:: locate_spots( const Real ymin, const Real ymax )
{
    assert( 0 == spots.size );
    Tracer *tracer = root;
    size_t from = 0;
    for( size_t i=0;i<size;++i,tracer=tracer->next)
    {
        const Real y = tracer->vertex.y;
        if( y>= ymin && y < ymax )
        {
            const size_t jump = i - from;
            from = i;
            spots.attach(tracer);
            spots.tail->jump  = jump;
            tracer->is_spot   = true;
        }
        else
            tracer->is_spot = false;
    }
}

void Bubble:: update_area() throw()
{
    const Tracer *p = root;
    Vertex  v0(0,0);
    Real ans= 0;
    for( size_t i=size;i>0;--i,p=p->next )
    {
        const Vertex v1 = v0 + p->edge;
        ans += v0.x * v1.y - v0.y * v1.x;
        v0 = v1;
    }
    area = Fabs(ans)/2;
}

void Bubble:: set_pressure( Real p ) throw()
{
    pressure = p;
    content  = area * p;
}

