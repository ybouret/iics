#ifndef JUNCTION_INCLUDED
#define JUNCTION_INCLUDED 1

#include "bubble.hpp"
#include "yocto/core/list.hpp"

//! bubble/mesh junctions
class Junction
{
public:
    Junction() throw();
    ~Junction() throw();
    static const int Vert = 1;
    static const int Horz = 2;
    
    
    Junction       *next;
    Junction       *prev;
    int             kind;      //!< Vert | Horz
    Vertex          vertex;    //!< localization
    unit_t          klo;       //!< lower indices
    unit_t          khi;       //!< upper indices
    const Bubble   *bubble;    //!< owner
    Real            alpha;     //!< weight of the distant tracer for curvature,etc...
    Real            curvature; //!< average curvature
    Real            pressure;  //!< from bubble and curvature
    Vertex          t;         //!< tangent
    Vertex          n;         //!< normal
    Real            gradP_t;   //!< tangential gradP
    typedef cache_of<Junction> Cache;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junction);
};

typedef cached_list<core::list_of, Junction> Junctions;

#endif

