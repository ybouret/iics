#ifndef SPOT_INCLUDED
#define SPOT_INCLUDED 1

#include "tracer.hpp"
#include "yocto/sequence/handle-list.hpp"

//! handle to a tracer
class Spot
{
public:
    Spot() throw();
    ~Spot() throw();
    
    Spot   *next;
    Spot   *prev;
    Tracer *handle;
    size_t  jump;             //!< encoding for MPI
    mutable coord2D klo;      //!< lower indices for handle->vertex on grid
    mutable coord2D kup;      //!< upper indices for handle->vertex on grid
    Vertex  gradP;
    Vertex  U;
    typedef cache_of<Spot> Cache;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Spot);
    
};

typedef handle_list<Spot> Spots;

#endif
