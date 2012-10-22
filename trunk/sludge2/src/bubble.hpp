#ifndef BUBBLE_INCLUDED
#define BUBBLE_INCLUDED 1

#include "spot.hpp"

typedef unsigned BubbleID;

class Bubble : public Tracers
{
public:
    explicit Bubble(BubbleID       bubble_id,
                    const PBC     &bubble_pbc,
                    Real          &bubble_lam,
                    Real          &bubble_gam,
                    Tracer::Cache &tcache,
                    Spot::Cache   &scache ) throw();
    
    virtual ~Bubble() throw();
    
    Bubble         *next;
    Bubble         *prev;
    const BubbleID  id;       //!< identifier
    const PBC      &pbc;      //!< periodic boundary conditions
    Real           &lam;      //!< spatial resolution
    Real           &gam;      //!< surface tension
    Spots           spots;    //!< associated spots on domain
    
    Real            pressure; //!< internal pressure
    Real            area;     //!< from tracers
    Real            content;  //!< pressure*area=content
    
    static const size_t IO_COUNT = 3; //!< pressure, area, content to dispatch
    
    void            set_pressure( Real p ) throw();
    void            update_area_fast() throw(); //!< use precomputed edges
    void            update_area_full() throw(); //!< don't use edges

    void clear() throw(); //!< no tracers, no spots
    
    //! partial hash
    void  hash( hashing::function &h ) const;
    
    //! full hash
    size_t get_hash( hashing::function &h) const;
    
    //! build_spots
    /**
     assume no initial spots
     */
    void locate_spots( const Real ymin, const Real ymax );
    
    
    //! apply PBC on tracers and use splines
    void compute_contour();
    
    
#if defined(HAS_MPI)
    //! dispatch tracers
    /**
     - empty spots
     - empty slave tracers
     - broadcast master tracers
     */
    void dispatch( const mpi &MPI );
    
    //! assemble tracers vertex from spots
    void assemble( const mpi &MPI );
#endif
    
    void map_circle( const Vertex &center, Real radius );
    void map_peanut( const Vertex &center, Real radius, Real alpha );
    
    void save_dat(   const string &filename ) const;
    void save_spots( const string &filename ) const;
    void save_vtk(   const string &filename )   const;
    void save_vtk_t( const string &filename ) const; //!< tangents
    void save_vtk_n( const string &filename ) const; //!< normals + curvature
    void save_vtk_shell( const string &filename ) const; //!< normals + lambda shell
    void rotate(Real alpha);
    void save_vtk_gt( const string &filename) const; //!< tangential gradient, after compute topology
    void save_vtk_gn( const string &filename) const; //!< normal gradent for the spots
    void save_vtk_g( const string &filename ) const; //!< grad P for the spots
    void save_vtk_u( const string &filename ) const; //!< U for the spots

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
};

#endif
