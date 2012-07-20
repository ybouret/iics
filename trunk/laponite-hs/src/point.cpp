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

#if 0
////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Point:: Pool:: Pool() throw() : CorePool()
{
    
}

Point:: Pool:: Pool( size_t cache_size ) : CorePool()
{
    while( size < cache_size ) store( new Point() );
}

Point:: Pool:: ~Pool() throw()
{
    while( size ) delete query();
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Point:: List:: List( Pool &cache ) throw() : CoreList(), cache_(cache)
{
}

Point:: List:: ~List() throw()
{
    empty();
}

Point * Point::List:: create()
{
    if( cache_.size > 0 )
    {
        Point *p = cache_.query();
        memset( &p->vertex,0,sizeof(Point)-YOCTO_OFFSET_OF(Point, vertex));
        return p;
    }
    else 
    {
        return new Point();
    }
}

void Point:: List:: empty() throw()
{
    while( size) cache_.store( pop_back() );
}

void Point:: List:: remove( Point *p ) throw()
{
    assert(p!=NULL);
    assert(this->owns(p));
    cache_.store( unlink(p) );
}
#endif

