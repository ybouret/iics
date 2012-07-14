#ifndef BUBBLES_INCLUDED
#define BUBBLES_INCLUDED 1


#include "bubble.hpp"


class Bubbles 
{
public:
    explicit Bubbles() throw();
    virtual ~Bubbles() throw();
    
    size_t         count() const throw();
    Bubble *       first() throw();
    const Bubble * first() const throw();
    
    void   none() throw();
    
    Bubble *create(); //!< return a new bubble, already in list
    void    create( size_t n ); //!< create n extra new bubbles
    
    
#if defined(HAS_MPI)
    void dispatch_all( mpi &MPI );
#endif
    
private:
    core::list_of<Bubble> b_list;
    core::pool_of<Bubble> b_pool;
    Point::Pool           pcache;
    Spot::Pool            scache;
    
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubbles);
};


#endif
