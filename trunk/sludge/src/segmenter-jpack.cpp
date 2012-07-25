#include "segmenter.hpp"


JPack:: ~JPack() throw() {}

JPack:: JPack() throw() :
vc(0),
lo(0),
id(0)
{
    
}

JPack:: JPack( const JPack &other ) throw() :
vc( other.vc ),
lo( other.lo ),
id( other.id )
{
    
}

JPack & JPack:: operator=( const JPack &other ) throw()
{
    vc = other.vc;
    lo = other.lo;
    id = other.id;
    return *this;
}   


void JPack:: Prepare( vector<JPack> &jpack, size_t n )
{
    const JPack zj;
    jpack.make(n,zj);
}

void JPack:: HEncode( vector<JPack> &jpack, const Segment::List &src )
{
    jpack.free();
    jpack.reserve(src.size);
    
    for(const Segment *s = src.head;s;s=s->next)
    {
        const Junction *J = s->handle;
        assert(J->bubble!=NULL);
        assert(J->bubble->id>0);
        {
            const JPack zj;
            jpack.push_back(zj);
        }
        JPack &j = jpack.back();
        j.vc = J->vertex.x;
        j.lo = J->lo;
        j.id = J->bubble->id;
    }
}

void JPack::HDecode(Segment::List &tgt, const vector<JPack> &jpack, Segmenter &segmenter, const Real y, Bubbles &bubbles)
{
    assert(0==tgt.size);
    for( size_t i=1; i <= jpack.size(); ++i )
    {
        const JPack &j = jpack[i];
        assert(j.id>0);
        Junction *J = segmenter.junctions.append();
        J->vertex.x = j.vc;
        J->vertex.y = y;
        J->lo       = j.lo;
        J->up       = j.lo+1;
        J->bubble   = bubbles[j.id];
        tgt.attach(J);
    }
}
