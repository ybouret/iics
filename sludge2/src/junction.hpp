#ifndef JUNCTION_INCLUDED
#define JUNCTION_INCLUDED 1

#include "bubble.hpp"
#include "yocto/core/list.hpp"

#define JUNCTION_TAG 0

//! bubble/mesh junctions
class Junction
{
public:
    typedef cache_of<Junction> Cache;

    Junction() throw();
    ~Junction() throw();
#if JUNCTION_TAG == 1
    static const int Vert = 1;
    static const int Horz = 2;
#endif
    
    Junction       *next;
    Junction       *prev;
#if JUNCTION_TAG == 1
    int             kind;      //!< Vert | Horz
    unit_t          tag;       //!< Vert=>i, Horz=>j
#endif
    Vertex          vertex;    //!< localization
    unit_t          klo;       //!< lower indices
    unit_t          khi;       //!< upper indices
    const Bubble   *bubble;    //!< owner
    Real            alpha;     //!< weight of the distant tracer for curvature,etc...
    Real            curvature; //!< average curvature
    Real            pressure;  //!< from bubble and curvature
    Vertex          t;         //!< tangent
    Vertex          n;         //!< normal
    mutable bool    visited;   //!< for orthonormal gradient
    mutable Vertex  g;         //!< gradP in cartesian coordinate, LOCALLY computed
    
    Real Peff( const Vertex &pos ) const throw(); //! pressure + (pos-vertex) * g
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junction);
};

typedef cached_list<core::list_of, Junction> Junctions;

#endif

