#ifndef BUBBLES_INCLUDED
#define BUBBLES_INCLUDED 1

#include "bubble.hpp"

class Bubbles
{
public:
    explicit Bubbles( const PBC &bubbles_pbc ) throw();
    virtual ~Bubbles() throw();
    
    Bubble       *first() throw();
    const Bubble *first() const throw();
    size_t        count() const throw();
    
    Bubble *append();         //!< one new bubble
    void    clear() throw();  //!< empty all
    void    create(size_t n); //!< clear and append n times
    
    size_t get_hash( hashing::function &H ) const;
    
    
    const PBC  &pbc;
    Real        lambda; //!< default is 1
    Real        gamma;  //!< default is 0
    
    void        update_topology();
    
#if defined(HAS_MPI)
    void dispatch( const mpi & MPI ); //!< dispatch all the bubbles
    void assemble( const mpi & MPI ); //!< assemble all the bubbles
#endif
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubbles);
    core::list_of<Bubble> bubbles;
    core::pool_of<Bubble> pool;
    Tracer::Cache         tcache;
    Spot::Cache           scache;
    

};

#endif

