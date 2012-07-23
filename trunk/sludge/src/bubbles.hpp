#ifndef BUBBLES_INCLUDED
#define BUBBLES_INCLUDED 1


#include "bubble.hpp"



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
    
    //! check all topologies, non crossing etc...
    /**
     TODO: ....
     */
    void check_topologies();
    
    //! find spots and compute geometries for active bubbles
    void check_geometries_within( Real y_lo, Real y_hi );
    
#if defined(HAS_MPI)
    //! master check and dispatch to slaves
    void check_and_dispatch_all(const mpi &MPI );
    
    //! each slave -> modified bubbles
    void assemble_all( const mpi &MPI );
#endif
    
private:
    core::list_of<Bubble> b_list;
    core::pool_of<Bubble> b_pool;
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubbles);
};


#endif
