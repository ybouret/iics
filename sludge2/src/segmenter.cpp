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

Junctions & Segmenter:: duplicated() throw()
{
    return duplicates;
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




