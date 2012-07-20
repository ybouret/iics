#include "point.hpp"
#include "yocto/core/offset.hpp"
#include <cstring>

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Point:: Point() throw()  :
vertex(),
s_next(0),
r_next(),
t(),
n(),
next(0), 
prev(0)
{
}

Point:: Point( const Point &other ) throw() : 
vertex( other.vertex) ,
s_next( other.s_next ),
t( other.t ),
n( other.n ),
next(0), prev(0)
{
    
}

Point::~Point() throw()
{
}

Point & Point:: operator=( const Point & other ) throw() 
{
    vertex    = other.vertex;
    s_next    = other.s_next;
    return *this;
}

void Point:: reset() throw()
{
    memset( &vertex,0,sizeof(Point)-YOCTO_OFFSET_OF(Point, vertex));
}


