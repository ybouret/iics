#include "bubble.hpp"
#include "yocto/code/utils.hpp"

#include <cstring>

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

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Spot:: Pool:: Pool() throw() : CorePool()
{
}

Spot:: Pool:: ~Pool() throw()
{
    while( size ) delete query();
}


////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Spot:: List:: List( Spot::Pool &cache ) throw() : CoreList(), cache_( cache )
{
    
}

Spot:: List:: ~List() throw()
{
    empty();
}

void Spot:: List:: empty() throw()
{
    while( size ) cache_.store( pop_back() );
}

void Spot::List:: append( Point *p )
{
    assert(p!=NULL);
    Spot *spot = cache_.size > 0 ? cache_.query() : new Spot;
    memset(spot,0,sizeof(Spot));
    spot->point = p;
    push_back(spot);
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Bubble:: Bubble( Point::Pool &cache, Spot::Pool &spot_cache ) throw() : 
Point::List( cache ),
lambda(1),
area(0),
spots( spot_cache )
{
}

Bubble::~Bubble() throw()
{
}


void Bubble:: update()
{
    assert(size>=3);
    assert(root!=NULL);
    area = 0;
    
    Point *p = root;
    Point *q = p->next;
    
    //--------------------------------------------------------------------------
    // first pass : refinement
    //--------------------------------------------------------------------------
    for( size_t i=0; i < size; ++i )
    {
        for(;;)
        {
            //------------------------------------------------------------------
            // compute length to next vertex
            //------------------------------------------------------------------
            const V2D pq(*p,*q);
            p->s_next = pq.norm();
            //std::cerr << "s_next=" << p->s_next << std::endl;
            
            //------------------------------------------------------------------
            // do we refine ?
            //------------------------------------------------------------------
            if( p->s_next > lambda )
            {
                //std::cerr << "..split" << std::endl;
                Point *I = create();
                I->x = p->x + 0.5 * pq.x;
                I->y = p->y + 0.5 * pq.y;
                insert_after(p, I);
                assert(p->next==I);
                q = I;
                continue;
            }
            
            break;
        }
        
        //----------------------------------------------------------------------
        // update area
        //----------------------------------------------------------------------
        area += p->x * q->y - p->y * q->x;
        
        
        //----------------------------------------------------------------------
        // next edge
        //----------------------------------------------------------------------
        p = q;
        q = q->next;
    }
    area = 0.5 * Fabs( area );
    //--------------------------------------------------------------------------
    // differential properties
    //--------------------------------------------------------------------------
        
}


double Bubble:: evaluate_area() const throw()
{
    double ans = 0;
    const Point *p = root;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        const Point *q = p->next;
        ans +=  p->x * q->y - p->y * q->x;
    }
    return 0.5 * Fabs(ans);
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

void Bubble:: build_spots(const Real y_lo, const Real y_up)
{
    spots.empty();
    Point *p = root;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        const double y = p->y;
        if( y_lo <= y && y < y_up )
            spots.append(p);
    }
}

