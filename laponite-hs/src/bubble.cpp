#include "bubble.hpp"
#include "yocto/code/utils.hpp"

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
    self      = other;
    domain    = other.domain;
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
        p->x = p->y = 0;
        p->domain = -1;
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

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Bubble:: Bubble( Point::Pool &cache ) throw() : 
Point::List( cache ),
lambda(1),
area(0)
{
}

Bubble::~Bubble() throw()
{
}


void Bubble:: Update() throw()
{
    area = 0;
}


void Bubble:: map_circle(const V2D &center, Real radius)
{
    assert(lambda>0);
    empty();
    const double theta_max = 2 * atan( lambda/(radius+radius) );
    const size_t nmin      = max_of<size_t>(3,size_t( ceil( numeric<Real>::two_pi/theta_max) ));
    std::cerr << "lambda=" << lambda << ", radius=" << radius << " => nmin=" << nmin << std::endl;
    const double dtheta = numeric<Real>::two_pi / nmin;
    for( size_t i=0; i < nmin; ++i )
    {
        const double theta = i * dtheta;
        Point       *p     = create();
        p->x = center.x + radius * Cos( theta );
        p->y = center.y + radius * Sin( theta );
        push_back(p);
    }
}

