#include "segmenter.hpp"

Segmenter:: ~Segmenter() throw()
{
    
}

Segmenter:: Segmenter( const Grid &g ) :
X( g.X() ),
Y( g.Y() ),
hseg(0),
vseg(0),
segcount( X.items + Y.items ),
segments( segcount, as_capacity ),
jcache()
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

Segment & Segmenter:: Horz( unit_t i) throw()
{
    assert(hseg);
    assert(i>=X.lower);
    assert(i<=X.upper);
    return *hseg[i];
}

Segment & Segmenter:: Vert( unit_t j) throw()
{
    assert(vseg);
    assert(j>=Y.lower);
    assert(j<=Y.upper);
    return *vseg[j];
}

