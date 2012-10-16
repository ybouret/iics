#ifndef TRACER_INCLUDED
#define TRACER_INCLUDED 1

#include "types.hpp"
#include "yocto/sequence/cached-list.hpp"
#include "yocto/core/clist.hpp"
#include "yocto/hashing/sha1.hpp"

class Bubble;
class Junction;

class Tracer
{
public:
    Tracer() throw();
    ~Tracer() throw();

    Tracer *next;
    Tracer *prev;
    
    Vertex          vertex;    //!< position (should be PBC) : +2
    Vertex          edge;      //!< vector to next->vertex   : +2
    Real            s2;        //!< |edge|^2                 : +1
    Real            s;         //!< |edge|^2                 : +1
    Vertex          t;         //!< local tangent vector     : +2
    Vertex          n;         //!< local normal vector      : +2
    Real            angle;     //!< tangent positive angle   : +1
    Real            curvature; //!< local curvature          : +1
    Bubble         *bubble;    //!< whose that ?
    bool            is_spot;   //!< default: false
    const Junction *jnext;     //!< a junction is between this and next tracer
    const Junction *jprev;     //!< a junction is between this and prev tracer
    
    static const size_t IO_COUNT = 12;
    
    typedef cache_of<Tracer>                     Cache;
    
    void hash( hashing::function &h ) const;
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Tracer);
};

typedef cached_list< core::clist_of, Tracer> Tracers; //!< circular list

#endif
