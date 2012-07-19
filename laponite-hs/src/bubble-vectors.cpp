#include "bubble.hpp"
#include "yocto/code/utils.hpp"

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
