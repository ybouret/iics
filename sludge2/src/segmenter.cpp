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
    
    hseg = &segments[1];
    vseg = hseg + Y.items;
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

void Segmenter:: save( const string &filename ) const
{
    assert(hseg);assert(vseg);
    ios::ocstream fp(filename,false);
    for( size_t i=1; i<=segments.size();++i)
    {
        const Segment &seg = *segments[i];
        for( const Junction *J = seg.head; J; J=J->next)
        {
            fp("%g %g\n", J->pos.x, J->pos.y);
        }
    }
}


