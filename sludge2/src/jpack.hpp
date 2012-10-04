#ifndef JPACK_INCLUDED
#define JPACK_INCLUDED 1

#include "junction.hpp"

//! Vertical Junction Packing
class JPack
{
public:
    JPack() throw(); //!< create an invalid Jpack (b=0)
    JPack( unit_t idx, const Junction *J ) throw();
    ~JPack() throw();
    JPack( const JPack &) throw();
    const unit_t   i; //!< abscissa index
    const BubbleID b; //!< which bubble
    const Real     y; //!< junction Y
    
private:
    YOCTO_DISABLE_ASSIGN(JPack);
};

#endif
