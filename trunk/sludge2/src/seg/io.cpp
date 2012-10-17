#include "../segmenter.hpp"
void Segmenter:: save( const string &filename ) const
{
    assert(hseg);assert(vseg);
    ios::ocstream fp(filename,false);
    for( size_t i=1; i<=segments.size();++i)
    {
        const Segment &seg = *segments[i];
        for( const Junction *J = seg.head; J; J=J->next)
        {
            fp("%g %g\n", J->vertex.x, J->vertex.y);
        }
    }
}

void Segmenter:: save_vtk_n( const string &filename ) const
{
    const unsigned n = num_junctions();
    ios::ocstream fp( filename, false );
    fp("# vtk DataFile Version 1.0\n");
    fp("Bubble Normals\n");
    fp("ASCII\n");
    fp("DATASET POLYDATA\n");
    fp("POINTS %u float\n", 2*n );
    for( size_t i=segcount;i>0;--i)
    {
        for( const Junction *p=segments[i]->head;p;p=p->next)
        {
            fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
            const Real scale = p->bubble->lam/2;
            const Real fac = scale * p->curvature;
            fp("%.15g %.15g 0\n",p->vertex.x + fac * p->n.x,p->vertex.y+ fac * p->n.y);
        }
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }
}

void Segmenter:: save_vtk_gt( const string &filename ) const
{
    const unsigned n = num_junctions();
    ios::ocstream fp( filename, false );
    fp("# vtk DataFile Version 1.0\n");
    fp("Bubble Normals\n");
    fp("ASCII\n");
    fp("DATASET POLYDATA\n");
    fp("POINTS %u float\n", 2*n );
    for( size_t i=segcount;i>0;--i)
    {
        for( const Junction *p=segments[i]->head;p;p=p->next)
        {
            fp("%.15g %.15g 0\n",p->vertex.x,p->vertex.y);
            const Real scale = 2*p->bubble->lam;
            //const Real scale = 1;
            const Real fac  = scale * p->gradP_t;
            fp("%.15g %.15g 0\n",p->vertex.x + fac * p->t.x,p->vertex.y+ fac * p->t.y);
        }
    }
    fp("\n");
    fp("LINES %u %u\n", n, 3*n );
    for( unsigned i=0; i < n; ++i )
    {
        fp("2 %u %u\n", 2*i, 2*i+1 );
    }
}




void Segmenter:: show_jvert(  ) const
{
    fprintf( stderr, "<JVert>\n");
    for( unit_t i=X.lower;i<=X.upper;++i)
    {
        const Junctions &jseg = Vert(i);
        if( jseg.size )
        {
            fprintf(stderr,"@i=%ld:", i);
            for( const Junction *J = jseg.head;J;J=J->next)
            {
                fprintf( stderr, " %g", J->vertex.y);
            }
            fprintf(stderr,"\n");
        }
    }
    fprintf( stderr, "</JVert>\n");
    
}
