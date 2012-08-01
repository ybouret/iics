#ifndef BUBBLES_INCLUDED
#define BUBBLES_INCLUDED 1


#include "bubble.hpp"
#include "yocto/auto-ptr.hpp"

class Rescaler;
class Bubbles 
{
public:
    explicit Bubbles( const Vertex &box ) throw();
    virtual ~Bubbles() throw();
    
    size_t        count() const throw();
    Bubble *      first() throw();
    const Bubble *first() const throw();
    
    const PBC     pbc;
    Tracer::Cache tcache;
    Spot::Cache   scache;
    Marker::Cache mcache;
    Real          lambda;
    
    Bubble *create();
    void    empty() throw();
    
       
    //! find spots and compute geometries for active bubbles
    void check_geometries_within( Real y_lo, Real y_hi );
    
#if defined(HAS_MPI)
    //! master check and dispatch to slaves
    void check_and_dispatch_all(const mpi &MPI, auto_ptr<Rescaler> &rescaler);
    
    //! each slave -> modified bubbles: update metrics and pressure
    void assemble_all( const mpi &MPI, auto_ptr<Rescaler> &rescaler );
#endif
    
    //! fill with markers from segmentation
    void fill( Array &B ) const;
    
    Bubble * operator[]( BubbleID id ) throw();
    
    
private:
    core::list_of<Bubble> b_list; //!< Bubble list
    core::pool_of<Bubble> b_pool; //!< Bubble pool
    vector<Bubble*>       bp_vec; //!< Bubble pointer vector
    
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubbles);
};


#endif
