#include "segmenter.hpp"

Segmenter:: ~Segmenter() throw()
{
    
}

Segmenter:: Segmenter( const Grid &g ) :
X( g.X() ),
Y( g.Y() ),
hseg(0),
vseg(0),
jcache(),
mcache(),
segcount( X.items + Y.items ),
segments( segcount, as_capacity ),
duplicates(jcache),
markers(mcache)
{
}

void Segmenter:: create()
{
    segments.free();
    assert(segments.capacity()>=segcount);
    hseg=0;
    vseg=0;
    
    for( unit_t j=Y.lower;j<=Y.upper;++j)
    {
        const Segment::Ptr sp( new Segment(Y[j],jcache));
        segments.push_back(sp);
    }
    
    for( unit_t i=X.lower;i<=X.upper;++i)
    {
        const Segment::Ptr sp( new Segment(X[i],jcache));
        segments.push_back(sp);
    }
    
    hseg  = &segments[1];
    vseg  = hseg + Y.items;
    hseg -= Y.lower;
    vseg -= X.lower;
    
}

Segment & Segmenter:: Horz( unit_t j) throw()
{
    assert(hseg);
    assert(j>=Y.lower);
    assert(j<=Y.upper);
    return *hseg[j];
}

const Segment & Segmenter::  Horz( unit_t j) const throw()
{
    assert(hseg);
    assert(j>=Y.lower);
    assert(j<=Y.upper);
    return *hseg[j];
}


Segment & Segmenter:: Vert( unit_t i) throw()
{
    assert(vseg);
    assert(i>=X.lower);
    assert(i<=X.upper);
    return *vseg[i];
}

const Segment & Segmenter:: Vert( unit_t i) const throw()
{
    assert(vseg);
    assert(i>=X.lower);
    assert(i<=X.upper);
    return *vseg[i];
}


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

size_t Segmenter:: num_junctions() const throw()
{
    size_t nj = 0;
    for( size_t i=segcount;i>0;--i)
    {
        nj += segments[i]->size;
    }
    return nj;
}

const Segments & Segmenter:: operator()(void) const throw()
{
    return segments;
}


void Segmenter:: save_vtk_n( const string &filename, Real scale) const
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


