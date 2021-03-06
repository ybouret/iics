#include "../bubble.hpp"

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
    const Real   scale  = lam;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
        const Real fac = scale * p->curvature;
        //const Real fac = scale;
        fp("%.15g %.15g 0\n",p->vertex.x + fac * p->n.x,p->vertex.y+ fac * p->n.y);
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }
}

void Bubble:: save_vtk_shell( const string &filename ) const
{
    const unsigned n = size;
    ios::ocstream fp( filename, false );
    fp("# vtk DataFile Version 1.0\n");
    fp("Bubble Normals\n");
    fp("ASCII\n");
    fp("DATASET POLYDATA\n");
    fp("POINTS %u float\n", 2*n );
    const Tracer *p     = root;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
        const Real fac = lam;
        fp("%.15g %.15g 0\n",p->vertex.x - fac * p->n.x,p->vertex.y-  fac * p->n.y);
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }

}

void Bubble:: save_vtk_gt( const string &filename ) const
{
    const unsigned n = size;
    ios::ocstream fp( filename, false );
    fp("# vtk DataFile Version 1.0\n");
    fp("Bubble Normals\n");
    fp("ASCII\n");
    fp("DATASET POLYDATA\n");
    fp("POINTS %u float\n", 2*n );
    const Tracer *p     = root;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
        const Real fac =  p->gt;
        fp("%.15g %.15g 0\n",p->vertex.x + fac * p->t.x,p->vertex.y +  fac * p->t.y);
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }
    
}

void Bubble:: save_vtk_gn( const string &filename ) const
{
    const unsigned n = spots.size;
    ios::ocstream fp( filename, false );
    fp("# vtk DataFile Version 1.0\n");
    fp("Bubble GradP on Spots\n");
    fp("ASCII\n");
    fp("DATASET POLYDATA\n");
    fp("POINTS %u float\n", 2*n );
    for( const Spot *spot = spots.head; spot; spot=spot->next )
    {
        const Tracer *p = spot->handle;
        fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
        const Real fac = spot->gn;
        fp("%.15g %.15g 0\n",p->vertex.x + fac * p->n.x,p->vertex.y + fac * p->n.y);
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }
}


void Bubble:: save_vtk_g( const string &filename ) const
{
    const unsigned n = spots.size;
    ios::ocstream fp( filename, false );
    fp("# vtk DataFile Version 1.0\n");
    fp("Bubble GradP on Spots\n");
    fp("ASCII\n");
    fp("DATASET POLYDATA\n");
    fp("POINTS %u float\n", 2*n );
    for( const Spot *spot = spots.head; spot; spot=spot->next )
    {
        const Tracer *p = spot->handle;
        fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
        const Real fac = 2*lam/gam;
        fp("%.15g %.15g 0\n",p->vertex.x + fac * spot->gradP.x,p->vertex.y + fac * spot->gradP.y);
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }
}

void Bubble:: save_vtk_u( const string &filename ) const
{
    const unsigned n = spots.size;
    ios::ocstream fp( filename, false );
    fp("# vtk DataFile Version 1.0\n");
    fp("Bubble U on Spots\n");
    fp("ASCII\n");
    fp("DATASET POLYDATA\n");
    fp("POINTS %u float\n", 2*n );
    for( const Spot *spot = spots.head; spot; spot=spot->next )
    {
        const Tracer *p = spot->handle;
        fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
        //const Real fac = 2*lam/gam;
        const Real fac = 1;
        fp("%.15g %.15g 0\n",p->vertex.x + fac * spot->U.x,p->vertex.y + fac * spot->U.y);
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }
}



