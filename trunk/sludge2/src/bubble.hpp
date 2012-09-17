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
                    Tracer::Cache &tcache,
                    Spot::Cache   &scache ) throw();
    
    virtual ~Bubble() throw();
    
    Bubble         *next;
    Bubble         *prev;
    const BubbleID  id;    //!< identifier
    const PBC      &pbc;   //!< periodic boundary conditions
    Real           &lam;   //!< spatial resolution
    Spots           spots; //!< associated spots on domain
    
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
    
    
    //! apply PBC on tracers
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
    
    void save_dat(   const string &filename ) const;
    void save_spots( const string &filename ) const;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
};

#endif
