#ifndef SPOT_INCLUDED
#define SPOT_INCLUDED 1

#include "tracer.hpp"
#include "yocto/core/list.hpp"

class Spot 
{
public:
    Spot() throw();
    ~Spot() throw();
    
    Tracer *handle; //!< pointing
    size_t  jump;   //!< encoding Tracer in Tracer::List
    Spot   *next;
    Spot   *prev;
    void reset() throw();
    
    typedef cache_of<Spot>                  Cache;
    typedef cached_list<core::list_of,Spot> List;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Spot);
};

#endif
