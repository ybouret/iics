#ifndef SEGMENTER_INCLUDED
#define SEGMENTER_INCLUDED 1

#include "segment.hpp"
#include "marker.hpp"
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
    
    //! get horizontal segment @Y[j]
    const Segment & Horz( unit_t j) const throw();
    
    
    //! get vertical  segment @X[i]
    Segment & Vert( unit_t i) throw();
    
    const Array1D  &X;
    const Array1D  &Y;
    
    void locate_vertex( const Vertex &v, coord2D &klo, coord2D &khi ) const;
    void process( const Bubbles &bubbles );
    
    //! fill the bubble array and compute the markers
    /**
     then need a MPI sync
     */
    void build( Array &B );
    
    //! fill the pressure array with the bubbles pressure
    void pressurize( Array &P ) const;
    
    //! save junctions coordinates, gnuplot style
    void save( const string &filename ) const;
    
    void save_vtk_n( const string &filename ) const;
    
    
    static
    void locate_value( const Real z, const Array1D &Z, unit_t &klo, unit_t &kup ) throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Segmenter);
    Segment::Ptr   *hseg;
    Segment::Ptr   *vseg;
    Junction::Cache jcache;
    Marker::Cache   mcache;
    const size_t    segcount;
    Segments        segments;
    void process_bubble( const Bubble *bubble, const Real half );
    void process_spot( const Spot *spot, const Real half);
    void compute_junctions(const Spot   *spot,
                           const Vertex &vertex,
                           const Vertex &target,
                           const Tracer *to
                           );
    Markers         markers;
};

#endif
