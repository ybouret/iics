#ifndef SEGMENTER_INCLUDED
#define SEGMENTER_INCLUDED 1

#include "segment.hpp"
#include "bubbles.hpp"

#include "yocto/sequence/vector.hpp"

typedef layout2D           Layout;
typedef rmesh<Real,Layout> Grid;
typedef array1D<Real>      Array1D;

class Segmenter 
{
public:
    explicit Segmenter( const Grid &grid );
    virtual ~Segmenter() throw();
    
    const Array1D   &X;
    const Array1D   &Y;
    const Array1D   &dX;
    const Array1D   &dY;
    Junction::Cache j_cache;
    Junction::List  junctions;
    Segment::Cache  s_cache;
    
    void allocate_segments();

    Segment::List *horizontal; //!< Y.lower->Y.upper
    Segment::List *vertical;   //!< X.lower->X.upper
    
    
    //! empty segments and junctions
    void clear() throw();
    
    //! process a Tracer spotted on the grid
    void process_tracer( Tracer *p );
    
    //!
    void process_bubble( Bubble *bubble );
    
    //! should clear_segments() before
    void process_bubbles( Bubbles &bubbles );
    
    //! sort segments by junction->lo
    void sort_segments();
    
    void save_junctions( const string &filename ) const;
    
private:
    vector<Segment::List> segments;
    YOCTO_DISABLE_COPY_AND_ASSIGN(Segmenter);
};


#endif
