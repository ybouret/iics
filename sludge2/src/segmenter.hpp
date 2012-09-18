#ifndef SEGMENTER_INCLUDED
#define SEGMENTER_INCLUDED 1

#include "segment.hpp"
#include "bubbles.hpp"

class Segmenter
{
public:
    
    //! reserve memory from grid size
    explicit Segmenter( const Grid &g );
    virtual ~Segmenter() throw();
    
    //! create segments from a valid grid
    void create();
    
    //! get horizontal segment @X[i]
    Segment & Horz( unit_t i) throw();
    
    //! get vertical segment @Y[j]
    Segment & Vert( unit_t j) throw();
    
    const Array1D  &X;
    const Array1D  &Y;
    
    void locate_vertex( const Vertex &v, coord2D &klo, coord2D &khi ) const;
    void process( const Bubbles &bubbles );
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Segmenter);
    Segment::Ptr   *hseg;
    Segment::Ptr   *vseg;
    const size_t    segcount;
    Segments        segments;
    Junction::Cache jcache;
    void process1( const Bubble *bubble );
    
};

#endif
