#ifndef BUBBLES_INCLUDED
#define BUBBLES_INCLUDED 1


#include "bubble.hpp"

typedef core::list_of<Bubble> BubbleList;

class Bubbles : public BubbleList
{
public:
    explicit Bubbles() throw();
    virtual ~Bubbles() throw();
    
    
private:
    Point::Pool pcache;
    Spot::Pool  scache;
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubbles);
};


#endif
