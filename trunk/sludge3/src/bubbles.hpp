#ifndef SLUDGE_BUBBLES_INCLUDED
#define SLUDGE_BUBBLES_INCLUDED 1

#include "bubble.hpp"

class Bubbles : public core::list_of<Bubble>
{
public:
    static const int Tag;
    
    explicit Bubbles() throw();
    virtual ~Bubbles() throw();
    
    Real    lambda; //!< shared for bubbles
    Real    gamma;  //!< shared fro bubbles
    
    Bubble *append();
    
    void hash( Hasher &h ) const throw();
    
    void regularize_all();
    void collect_all_markers(Real ymin,Real ymax);
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubbles);
};

#endif

