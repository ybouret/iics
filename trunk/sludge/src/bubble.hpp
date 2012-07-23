#ifndef BUBBLE_INCLUDED
#define BUBBLE_INCLUDED 1

#include "spot.hpp"

class Bubble  : public Tracer::List
{
public:
    explicit Bubble(Real                &lambda_ref, 
                    const PBC           &pbc_ref,
                    Tracer::Cache       &tracer_cache,
                    Spot::Cache         &spot_cache  
                    ) throw();
    virtual ~Bubble() throw();
    
    const Real         &lambda;
    const PBC          &pbc;
    Real                area;
    Spot::List          spots;    
    
    void clear() throw(); //!< empty() and spots.empty()
    
    //! compute initial tracer geometry
    /**
     - pbc for vertices and edges, and compute s2 and s
     - compute_geometry()
     */
    void raw_initialize();
    
    //! upgrade topology
    /**
     auto refinement
     */
    void upgrade_topology();

    
    //! compute area, frenet and curvatures for valid vertices, egdes, s2, s
    void compute_geometry();
    
    //! map to a circle and raw initialize
    void map_circle(const Vertex &center, Real radius);
    
    //! map to a peanut and raw initialize
    void map_peanut( const Vertex &center, Real radius, Real alpha );
    
    void save_dat( const string &filename ) const;
    void save_vtk( const string &filename ) const;
    void save_vtk_t( const string &filename ) const;
    void save_vtk_n( const string &filename ) const;

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
};

#endif
