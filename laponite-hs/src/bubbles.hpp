#ifndef BUBBLES_INCLUDED
#define BUBBLES_INCLUDED 1


#include "bubble.hpp"


class Bubbles 
{
public:
    explicit Bubbles(const V2D &box_dim) throw();
    virtual ~Bubbles() throw();
    
    size_t         count() const throw();
    Bubble *       first() throw();
    const Bubble * first() const throw();
    
    void   none() throw();
    
    Bubble *append(); //!< return a new bubble, appended to list
    void    append( size_t n ); //!< append n extra new bubbles
    
    void    update_topologies(); //!< for all bubbles
    
    //! find spots and compute values only for active bubbles
    void    spots_and_values_within(const Real y_lo, const Real y_up);
    
#if defined(HAS_MPI)
    void dispatch_all( const mpi &MPI ); //!< broacast master->everybody
    void assemble_all( const mpi &MPI ); //!< collect changed points
#endif
    
    const V2D  Length;
    const Real lambda; //!< critical length for bubbles, default is 1.0
    
    
private:
    core::list_of<Bubble> b_list;
    core::pool_of<Bubble> b_pool;
    
public:
    Point::Pool           pcache;
    Spot::Pool            scache;
    GridMarker::Pool      gcache;
    
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubbles);
};


#endif
