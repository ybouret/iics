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
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
};

#endif
