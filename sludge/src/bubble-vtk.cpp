#include "bubble.hpp"
#include "yocto/ios/ocstream.hpp"

void Bubble:: save_dat( const string &filename ) const
{
    ios::ocstream fp( filename, false );
    const Tracer *p = root;
    for( size_t i=size;i>0;--i,p=p->next )
    {
        fp("%.15g %.15g\n",p->vertex.x,p->vertex.y);
    }
    fp("%.15g %.15g\n",p->vertex.x,p->vertex.y);
}


void Bubble::save_spots( const string &filename ) const
{
    ios::ocstream fp( filename, false );
    for( const Spot *spot = spots.head; spot; spot=spot->next )
    {
        const Vertex &v = spot->handle->vertex;
        fp("%.15g %.15g\n",v.x,v.y);
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
    const Tracer *p     = root;
    const Real   scale = lambda * 0.5;
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
    const Real   scale = lambda * 0.5;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
        fp("%.15g %.15g 0\n",p->vertex.x + p->curvature * scale * p->n.x,p->vertex.y+ p->curvature * scale * p->n.y);
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }
}

void Bubble:: save_inside( const string &filename, const Grid &grid ) const
{
    ios::ocstream fp( filename, false );
    const Array1D &X = grid.X();
    const Array1D &Y = grid.Y();
    for( const Marker *m = markers.head; m; m=m->next )
    {
        fp("%g %g\n", X[m->coord.x], Y[m->coord.y]);
    }
    
}

