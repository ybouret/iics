#include "bubble.hpp"
#include "yocto/code/utils.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

Bubble:: Bubble( Real L, Point::Pool &pcache, Spot::Pool &scache ) throw() : 
Point::List( pcache ),
pbc(L),
lambda(1),
area(0),
spots( scache ),
active( true ),
next(0),
prev(0)
{
}

Bubble::~Bubble() throw()
{
}



void Bubble:: update_contour()
{
    assert(size>=3);
    assert(root!=NULL);
    
    //--------------------------------------------------------------------------
    // pass 0: pbc on vertices
    //--------------------------------------------------------------------------
    Point *p = root;
    
    for( size_t i=size;i>0;--i,p=p->next)
        pbc(p->vertex);
    
    p        = root;
    Point *q = p->next;
    
    //--------------------------------------------------------------------------
    // pass : refinement
    //--------------------------------------------------------------------------
    for( size_t i=0; i < size; ++i )
    {
        
        for(;;)
        {
            //------------------------------------------------------------------
            // compute length to next vertex, with pbc
            //------------------------------------------------------------------
            V2D pq(p->vertex,q->vertex);
            pbc(pq);
            p->s_next = pq.norm(); assert(p->s_next>0);
            
            //------------------------------------------------------------------
            // do we refine ?
            //------------------------------------------------------------------
            if( p->s_next > lambda )
            {
                Point *I = create();
                I->vertex.x = p->vertex.x + 0.5 * pq.x;
                I->vertex.y = p->vertex.y + 0.5 * pq.y;
                pbc(I->vertex);
                insert_after(p, I);
                assert(p->next==I);
                assert(I->prev==p);
                q = I;
                continue;
            }
            
            //------------------------------------------------------------------
            // ok, keep that in mind...
            //------------------------------------------------------------------
            p->r_next = pq;
            assert(p->s_next>0);
            break;
        }
        
        //----------------------------------------------------------------------
        // next edge
        //----------------------------------------------------------------------
        p = q;
        q = q->next;
    }
    
    p = root;
    for( size_t i=0; i < size; ++i, p=p->next)
    {
        assert(p->s_next>0);
    }
    
}

#if 0
double Bubble:: evaluate_area() const throw()
{
    double       ans = 0;
    const Point *p   = root;
    V2D          v0(0,0); // translate p to origin
    for( size_t i=size;i>0;--i,p=p->next)
    {
        // construct next point by effective difference vector
        const V2D v1 = v0 + p->r_next;
        ans += v0.x * v1.y - v0.y * v1.x;
        v0 = v1;
    }
    return 0.5 * Fabs(ans);
}
#endif

void Bubble:: compute_values() throw() 
{ 
    
    area         = 0;
    Point   *p   = root;
    V2D          v0(0,0); // translate p to origin
    for( size_t i=size;i>0;--i,p=p->next)
    {
        //----------------------------------------------------------------------
        // construct next point by effective difference vector
        //----------------------------------------------------------------------
        const V2D    r1 = p->r_next;
        const Real   s1 = p->s_next; assert(s1>0);
        const V2D    v1 = v0 + r1;
        
        //----------------------------------------------------------------------
        // update area
        //----------------------------------------------------------------------
        area += v0.x * v1.y - v0.y * v1.x;
        
        //----------------------------------------------------------------------
        // forward current position
        //----------------------------------------------------------------------
        v0 = v1;
        
        //----------------------------------------------------------------------
        // construct tangent vector
        //----------------------------------------------------------------------
        const Point *q  = p->prev;
        const V2D    r0 = q->r_next;
        const Real   s0 = q->s_next; assert(s0>0);
        const V2D    tp = s0 * r1  + s1 * r0; //!< tangent vector
        const Real   tp_norm = tp.norm();     //!< its norm
        p->t   = (1/tp_norm) * tp;
        
        //----------------------------------------------------------------------
        // deduce normal vector
        //----------------------------------------------------------------------
        p->n.x = - p->t.y;
        p->n.y =   p->t.x;
        
    }
    area = 0.5 * Fabs( area );
    
    
    //----------------------------------------------------------------------
    // compute curvature
    //----------------------------------------------------------------------
    p = root;
    for( size_t i=size;i>0;--i,p=p->next)
    {
        const V2D  t1 = p->next->t;
        const Real s1 = p->s_next; assert(s1>0);
        
        const V2D  t0 = p->prev->t;
        const Real s0 = p->prev->s_next; assert(s0>0);
        const V2D  tmp = s0 * t1 - s1 * t0;
        const Real norm_kappa = tmp.norm() / (s0*s1);
        p->kappa = sign_of( V2D::dot_(tmp, p->n) ) * norm_kappa;
    }
}

void Bubble:: map_circle(const V2D &center, Real radius)
{
    assert(lambda>0);
    empty();
    const double theta_max = 2 * atan( lambda/(radius+radius) );
    const size_t nmin      = max_of<size_t>(3,size_t( ceil( numeric<Real>::two_pi/theta_max) ));
    //std::cerr << "lambda=" << lambda << ", radius=" << radius << " => nmin=" << nmin << std::endl;
    const double dtheta = numeric<Real>::two_pi / nmin;
    const double theta0 = numeric<Real>::two_pi * Alea();
    for( size_t i=0; i < nmin; ++i )
    {
        const double theta = i * dtheta + theta0;
        Point       *p     = create();
        p->vertex.x = center.x + radius * Cos( theta );
        p->vertex.y = center.y + radius * Sin( theta );
        push_back(p);
    }
}

void Bubble:: map_astroid( const V2D &center, Real radius )
{
    assert(lambda>0);
    empty();
    const size_t ntop = max_of<size_t>(1,0.71*radius/lambda);
    const double dtheta = numeric<Real>::pi / (2*ntop);
    const size_t n = 4*ntop;
    const double theta0 = numeric<Real>::two_pi * Alea();
    for( size_t i=0; i < n; ++i )
    {
        const double theta = i * dtheta + theta0;
        Point       *p     = create();
        const Real C = Cos( theta );
        const Real S = Sin( theta );
        p->vertex.x = center.x + radius * C*C*C;
        p->vertex.y = center.y + radius * S*S*S;
        push_back(p);
    }
}

void Bubble:: find_spots_within(const Real y_lo, const Real y_up)
{
    spots.empty();
    Point *p = root;
    size_t last_index = 0;
    for( size_t i=0;i<size;++i,p=p->next)
    {
        const double y = p->vertex.y;
        if( y_lo <= y && y <= y_up )
        {
            spots.append(p);           
            spots.tail->jump = i-last_index;
            last_index = i;
        }
    }
    active = spots.size > 0 ;
}

