#ifndef TRACER_INCLUDED
#define TRACER_INCLUDED 1

#include "types.hpp"
#include "yocto/core/handle-list.hpp"
#include "yocto/core/clist.hpp"

class Tracer 
{
public:
    Tracer() throw();
    ~Tracer() throw();
    
    Vertex  vertex;    //!< to be broadcasted
    Vertex  edge;      //!< to next Tracer, broadcasted
    Real    s;         //!< length>0 to the next Tracer, broadcasted
    Vertex  t;         //!< tangent vector
    Vertex  n;         //!< normal vector
    Real    curvature; //!< curvature
    Tracer *next;
    Tracer *prev;
    
    void set_normal() throw(); //!< from t
    void reset() throw();
    
    typedef handle_of<Tracer>::node_type       Handle;       //!< managed by a bubble
    typedef handle_of<Tracer>::cache_type      HandleCache;  //!< its cache
    typedef handle_of<Tracer>::list_type       HandleList;   //!< its list
    
    typedef cache_of<Tracer>                   Cache;
    typedef cached_list<core::clist_of,Tracer> List;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Tracer);
};


#endif
