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
            tracer->is_spot   = true;
        }
        else
            tracer->is_spot = false;
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

void Bubble:: save_vtk( const string &filename ) const
{
    const unsigned n = size;
    ios::ocstream fp( filename, false );
    fp("# vtk DataFile Version 1.0\n");
    fp("Bubble\n");
    fp("ASCII\n");
    fp("DATASET POLYDATA\n");
    fp("POINTS %u float\n", n );
    const Tracer *p = root;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", i, (i+1) % n );
    }
}

void  Bubble:: save_vtk_t( const string &filename ) const
{
    const unsigned n = size;
    ios::ocstream fp( filename, false );
    fp("# vtk DataFile Version 1.0\n");
    fp("Bubble Tangents\n");
    fp("ASCII\n");
    fp("DATASET POLYDATA\n");
    fp("POINTS %u float\n", 2*n );
    const Tracer *p    = root;
    const Real   scale = lam ;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
        fp("%.15g %.15g 0\n",p->vertex.x + scale * p->t.x,p->vertex.y+ scale * p->t.y);
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }
}

void  Bubble:: save_vtk_n( const string &filename ) const
{
    const unsigned n = size;
    ios::ocstream fp( filename, false );
    fp("# vtk DataFile Version 1.0\n");
    fp("Bubble Normals\n");
    fp("ASCII\n");
    fp("DATASET POLYDATA\n");
    fp("POINTS %u float\n", 2*n );
    const Tracer *p     = root;
    const Real   scale = lam/2;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
        fp("%.15g %.15g 0\n",p->vertex.x + p->curvature * scale * p->n.x,p->vertex.y+ p->curvature * scale * p->n.y);
        //fp("%.15g %.15g 0\n",p->vertex.x + scale * p->n.x,p->vertex.y + scale * p->n.y);
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }
}



