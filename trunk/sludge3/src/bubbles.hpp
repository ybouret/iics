#ifndef SLUDGE_BUBBLES_INCLUDED
#define SLUDGE_BUBBLES_INCLUDED 1

#include "bubble.hpp"
#include "yocto/core/list.hpp"

class Bubbles : public core::list_of<Bubble>
{
public:
    static const int Tag;
    
    explicit Bubbles( Real &lam ) throw();
    virtual ~Bubbles() throw();
    const Real &lambda;
    Bubble *append();
    
    void hash( Hasher &h ) const throw();
    
    void regularize();
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubbles);
};

#endif

