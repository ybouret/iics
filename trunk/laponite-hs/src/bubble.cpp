#include "bubble.hpp"
#include "yocto/code/utils.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Bubble:: Bubble( Point::Pool &pcache, Spot::Pool &scache ) throw() : 
Point::List( pcache ),
lambda(1),
area(0),
spots( scache )
{
}

Bubble::~Bubble() throw()
{
}


void Bubble:: update_points()
{
    assert(size>=3);
    assert(root!=NULL);
    // area=0
    
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
        //area += p->x * q->y - p->y * q->x;
        
        
        //----------------------------------------------------------------------
        // next edge
        //----------------------------------------------------------------------
        p = q;
        q = q->next;
    }
    //area = 0.5 * Fabs( area );
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

void Bubble:: update_area() throw() { area = evaluate_area(); }

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
    active = spots.size > 0 ;
}

