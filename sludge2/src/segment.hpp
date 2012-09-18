#ifndef SEGMENT_INCLUDED
#define SEGMENT_INCLUDED

#include "junction.hpp"
#include "yocto/shared-ptr.hpp"
#include "yocto/sequence/vector.hpp"


class Segment : public Junctions
{
public:
    explicit Segment( const Real x_or_y, Junction::Cache &jcache ) throw();
    virtual ~Segment() throw();
    
    const Real value;
    
    typedef shared_ptr<Segment> Ptr;
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Segment);
};

typedef vector<Segment::Ptr> Segments;

#endif
