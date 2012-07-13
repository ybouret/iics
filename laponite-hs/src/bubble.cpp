#include "bubble.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Point:: Point() throw()  : V2D(), domain(-1), next(0), prev(0)
{
}

Point:: Point( const Point &other ) throw() : V2D(other), domain(other.domain), next(0), prev(0)
{
    
}

Point::~Point() throw()
{
}

Point & Point:: operator=( const Point & other ) throw() 
{
    V2D &self = *this;
    self   = other;
    domain = other.domain;
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
    while( size) cache_.store( pop_back() );
}

Point * Point::List:: New()
{
    if( cache_.size > 0 )
    {
        Point *p = cache_.query();
        p->x = p->y = 0;
        p->domain = -1;
        return p;
    }
    else 
    {
        return new Point();
    }
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Polygon:: Polygon( Point::Pool &cache ) throw() : 
Point::List( cache ),
area(0)
{
}

Polygon:: ~Polygon() throw()
{
}


void Polygon:: Update() throw()
{
    area = 0;
    
}


