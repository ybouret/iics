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
    
private:
    core::list_of<Bubble> b_list;
    core::pool_of<Bubble> b_pool;
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubbles);
};


#endif
