#ifndef BUBBLE_INCLUDED
#define BUBBLE_INCLUDED 1

#include "spot.hpp"
#include "marker.hpp"


#if defined(HAS_MPI)
#include "yocto/mpi/mpi.hpp"
#endif

typedef unsigned BubbleID;

class Bubble  : public Tracer::List
{
public:
    explicit Bubble(Real                &lambda_ref, 
                    const PBC           &pbc_ref,
                    Tracer::Cache       &tracer_cache,
                    Spot::Cache         &spot_cache,
                    Marker::Cache       &marker_cache
                    ) throw();
    virtual ~Bubble() throw();
    
    BubbleID            id;
    const Real         &lambda;
    const PBC          &pbc;
    Real                pressure; //!< broadcaster
    Real                area;     //!< broadcasted
    Real                content;  //!< pressure * area = content, broadcasted
    Spot::List          spots;
    Marker::List        markers; //!< computed by spots
    bool                active;
    Bubble             *next;
    Bubble             *prev;
    
    static const size_t IO_COUNT = 3; //!< pressure,area,content
    
    void clear() throw(); //!< empty() and spots.empty()
    
    //! compute initial tracer geometry
    /**
     - pbc for vertices and edges, and compute s2 and s
     - compute_geometry()
     */
    void raw_initialize();
    
    
    //! compute are from valid PBC vertices and egges
    void compute_area();
    
    //! compute area, frenet and curvatures for valid vertices, egdes, s2, s
    /**
     upgrade the pressure as well
     */
    void compute_geometry();
    
    
       
    //! compute content accordingly, area must be valid
    void set_pressure( Real pres);
        
    //! map to a circle and raw initialize
    void map_circle(const Vertex &center, Real radius);
    
    //! map to a peanut and raw initialize
    void map_peanut( const Vertex &center, Real radius, Real alpha );
    
    void save_dat( const string &filename ) const;
    void save_spots( const string &filename ) const;
    void save_vtk( const string &filename ) const;
    void save_vtk_t( const string &filename ) const;
    void save_vtk_n( const string &filename ) const;
    void save_inside( const string &filename, const Grid &grid ) const;
    
    
    //! tag each tracer and put it in spots if possible
    /**
     y_lo <= y < y_hi
     */
    void collect_spots_within( Real y_lo, Real y_hi );
    
#if defined(HAS_MPI)
    //! master topology -> slaves
    void dispatch_topology( const mpi &MPI );

    //! slave spots -> master
    void assemble_topology( const mpi &MPI );
#endif
    
    
    void translate( const Vertex &v );
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
};

#endif
