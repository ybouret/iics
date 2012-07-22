#ifndef POINT_INCLUDED
#define POINT_INCLUDED 1

#include "types.hpp"
//#include "yocto/core/pool.hpp"
#include "yocto/core/clist.hpp"
#include "yocto/core/cached-list.hpp"

//! a point on a border of a bubble
class Bubble;
class Point 
{
public:
    V2D      vertex;  //!< x,y (PBC)                         : I/O +2
    Real     s_next;  //!< distance to next point            : I/O +1
    V2D      r_next;  //!< vector to next point (using PBC)  : I/O +2
    V2D      t;       //!< tangent vector
    V2D      n;       //!< normal vector
    Real     kappa;   //!< curvature
    U2D      pos;     //!< lower coordinate of including mesh   X[pos.x] <= vertex.x < X[pos.y] and Y[pos.y] <= vertex.y < Y[pos.y] 
    V2D      w;       //!< bilinear interpolation weights
    Bubble  *bubble;  //!< owner, set by bubble...
    
    static const size_t IO_COUNT = 5;
    
    explicit Point() throw(); //!< x=y=0, domain=-1, d2next=0
    virtual ~Point() throw(); //!< do nothing
    
    Point( const Point & ) throw(); //! copy
    Point & operator=( const Point & ) throw();
    
    void reset() throw();
    
    Point *next; //!< for linked list
    Point *prev; //!< for linked list
    
    typedef cache_of<Point>                   Pool;
    typedef cached_list<core::clist_of,Point> List;
};



#endif
