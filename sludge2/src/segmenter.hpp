#ifndef SEGMENTER_INCLUDED
#define SEGMENTER_INCLUDED 1

#include "segment.hpp"


class Segmenter
{
public:
    
    //! reserve memory from grid size
    explicit Segmenter( const Grid &g );
    virtual ~Segmenter() throw();
    
    //! create segments from a valid grid
    void create();
    
    Segment & Horz( unit_t i) throw();
    Segment & Vert( unit_t j) throw();
    
    const Array1D  &X;
    const Array1D  &Y;

    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Segmenter);
      Segment::Ptr   *hseg;
    Segment::Ptr   *vseg;
    const size_t    segcount;
    Segments        segments;
    Junction::Cache jcache;
   
};

#endif
