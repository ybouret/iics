#ifndef TRACER_INCLUDED
#define TRACER_INCLUDED 1

#include "types.hpp"
#include "yocto/core/cached-list.hpp"
#include "yocto/core/clist.hpp"

class Bubble;

class Tracer 
{
public:
    Tracer() throw();
    ~Tracer() throw();
   
    
    Vertex  vertex;    //!< to be           broadcasted +2
    Vertex  edge;      //!< to next Tracer, broadcasted +2
    Real    s2;        //!< |edge|^2,       broadcasted +1 
    Real    s;         //!< |edge|,         broadcasted +1
    
    Tracer *next;
    Tracer *prev;
    
    Vertex  t;         //!< tangent vector, computed by process
    Vertex  n;         //!< normal vector,  computed by process
    Real    curvature; //!< curvature,      computed by process
    Coord   gLower;    //!< position on grid, lower indices
    Coord   gUpper;    //!< position on grid, upper indices
    Vertex  bw;        //!< bilinear interpolation weights
    Bubble *bubble;    //!< whose that ?
    bool    is_spot;   //!< spotted on the grid ?
    
    static const size_t IO_COUNT = 6;
    
    void set_normal() throw(); //!< from t
    void reset() throw();
    
    //! compute t and n once vertex, edges, s and s2 are available
    void compute_frenet();
    
    //! compute curvature once compute_frenet was computed
    void compute_curvature();
    
    typedef cache_of<Tracer>                   Cache;
    typedef cached_list<core::clist_of,Tracer> List;   //! circular chained list
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Tracer);
};


#endif
