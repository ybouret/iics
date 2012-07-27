#ifndef SPOT_INCLUDED
#define SPOT_INCLUDED 1

#include "tracer.hpp"
#include "yocto/core/list.hpp"
#include "yocto/core/handle-list.hpp"

//! handle to a tracer
class Spot 
{
public:
    Spot() throw();
    ~Spot() throw();
    
    Tracer *handle; //!< pointing
    size_t  jump;   //!< encoding Tracer in Tracer::List
    Vertex  U;      //!< velocity of the tracer
    Coord   gLower;    //!< position on grid, lower indices
    Coord   gUpper;    //!< position on grid, upper indices
    Vertex  bw;        //!< bilinear interpolation weights
    Spot   *next;
    Spot   *prev;
    void reset() throw();
    
    typedef cache_of<Spot>    Cache;
    typedef handle_list<Spot> List;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Spot);
};

#endif
