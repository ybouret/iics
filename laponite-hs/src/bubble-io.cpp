#include "bubble.hpp"


void Bubble:: save_dat( const string &filename ) const
{
    ios::ocstream fp( filename, false );
    const Point *p = root;
    for( size_t i=size;i>0;--i,p=p->next )
    {
        fp("%.15g %.15g\n",p->vertex.x,p->vertex.y);
    }
    fp("%.15g %.15g\n",p->vertex.x,p->vertex.y);
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
    const Point *p = root;
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
    const Point *p     = root;
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
    const Point *p     = root;
    const Real   scale = lambda * 0.5;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
        fp("%.15g %.15g 0\n",p->vertex.x + p->kappa * scale * p->n.x,p->vertex.y+ p->kappa * scale * p->n.y);
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }
}

