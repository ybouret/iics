#include "point.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Point:: Point() throw()  :
V2D(),
s_next(0),
next(0), prev(0)
{
}

Point:: Point( const Point &other ) throw() : 
V2D(other),
s_next( other.s_next ),
next(0), prev(0)
{
    
}

Point::~Point() throw()
{
}

Point & Point:: operator=( const Point & other ) throw() 
{
    V2D &self = *this;
    self      = other;
    s_next    = other.s_next;
    return *this;
}

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
        p->x      =  0;
        p->y      =  0;
        p->s_next =  0;
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
