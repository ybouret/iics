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
    
    //! get horizontal segment @Y[j]
    Segment & Horz( unit_t j) throw();
    
    //! get vertical  segment @X[i]
    Segment & Vert( unit_t i) throw();
    
    const Array1D  &X;
    const Array1D  &Y;
    
    void locate_vertex( const Vertex &v, coord2D &klo, coord2D &khi ) const;
    void process( const Bubbles &bubbles );
    
    void save( const string &filename ) const;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Segmenter);
    Segment::Ptr   *hseg;
    Segment::Ptr   *vseg;
    Junction::Cache jcache;
    const size_t    segcount;
    Segments        segments;
    void process_bubble( const Bubble *bubble );
    void process_spot( const Spot *spot);
    void compute_junctions( const Spot *spot, const Vertex &self, const Vertex &other );
    
};

#endif
