#ifndef BUBBLES_INCLUDED
#define BUBBLES_INCLUDED 1


#include "bubble.hpp"


class Bubbles 
{
public:
    explicit Bubbles(Real box_height) throw();
    virtual ~Bubbles() throw();
    
    size_t         count() const throw();
    Bubble *       first() throw();
    const Bubble * first() const throw();
    
    void   none() throw();
    
    Bubble *create(); //!< return a new bubble, appended to list
    void    create( size_t n ); //!< create n extra new bubbles
    
    void    update_topologies(); //!< for all bubbles
    
    //! find spots and compute values only for active bubbles
    void    spots_and_values_within(const Real y_lo, const Real y_up);
    
#if defined(HAS_MPI)
    void dispatch_all( const mpi &MPI ); //!< broacast master->everybody
    void assemble_all( const mpi &MPI );  //!< collect changed points
#endif
    
    const Real            Ly; //!< for bubble PBC
   
    
private:
    core::list_of<Bubble> b_list;
    core::pool_of<Bubble> b_pool;
    
public:
    Point::Pool           pcache;
    Spot::Pool            scache;
    
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubbles);
};


#endif
