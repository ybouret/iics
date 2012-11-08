#ifndef ARC_INCLUDED
#define ARC_INCLUDED 1

#include "./junction.hpp"

class Arc
{
public:
    Arc(const Junction *jp, const Junction *jn) throw();
    
    ~Arc() throw();
    
       
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Arc);
};


#endif

