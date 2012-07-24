#ifndef SEGMENTER_INCLUDED
#define SEGMENTER_INCLUDED 1

#include "segment.hpp"
#include "bubbles.hpp"

#include "yocto/sequence/vector.hpp"



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
    
    //! (re)allocate segments to match the grid
    void allocate_segments();
    
    //! direct access to segments
    Segment::List *horizontal; //!< Y.lower->Y.upper
    Segment::List *vertical;   //!< X.lower->X.upper
    
    
    //! empty segments and junctions
    void clear() throw();
    
    
    //! clear/process
    /**
     then upgrade PBC if required
     */
    void process_bubbles( Bubbles &bubbles );
    
    //! complete Y.lower and Y.upper in case of local PBC
    void horizontal_pbc();
    
    
    //! sort segments and assign borders by horizontal scanning
    /**
     valid only after a process_bubbles() !
     */
    void assign_markers();
    
    //! for debugging
    void save_junctions( const string &filename ) const;
    
    
    
private:
    vector<Segment::List> segments;
    YOCTO_DISABLE_COPY_AND_ASSIGN(Segmenter);
    //! process a Tracer spotted on the grid
    void process_tracer( Tracer *p );
    
    //! process one bubble
    void process_bubble( Bubble *bubble );
    
    //! sort segments by junction->lo
    void sort_segments();
    
    //! bissection algorithm
    static unit_t locate_point( Real, const Array1D &) throw();

    //! factorized code to detect the junctions
    void find_junctions( const Vertex &P, const Vertex &Q, const Vertex &vmin, const Vertex &vmax, Tracer *p );
    
};


#endif
