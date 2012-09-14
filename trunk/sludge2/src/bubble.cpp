#include "bubble.hpp"

Bubble:: ~Bubble() throw()
{
}

Bubble:: Bubble(BubbleID      bubble_id,
                const PBC    &bubble_pbc,
                Real          &bubble_lam,
                Tracer::Cache &tcache,
                Spot::Cache   &scache ) throw() :
Tracers(tcache),
next(0),
prev(0),
id(bubble_id),
pbc(bubble_pbc),
lam(bubble_lam),
spots(scache)
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
        }
    }
}

void Bubble:: save_dat( const string &filename ) const
{
    ios::ocstream fp( filename, false);
    const Tracer *tracer = root;
    for( size_t i=size;i>0;--i,tracer=tracer->next)
    {
        fp("%g %g\n", tracer->vertex.x, tracer->vertex.y);
    }
    fp("%g %g\n", tracer->vertex.x, tracer->vertex.y);

}

void Bubble:: save_spots( const string &filename ) const
{
    ios::ocstream fp( filename, false);
    for( const Spot *spot = spots.head; spot; spot=spot->next )
    {
        const Vertex &v = spot->handle->vertex;
        fp("%g %g\n", v.x, v.y);
    }
    
}

