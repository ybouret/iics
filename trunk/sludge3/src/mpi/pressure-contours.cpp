#include "workspace.hpp"

#if 0
#define BUBBLE_NONE    (0)
#define BUBBLE_BEFORE  (1)
#define BUBBLE_AFTER   (2)
#define BUBBLE_BOTH   ( BUBBLE_BEFORE | BUBBLE_AFTER )

void Workspace:: pressurize_horz()
{
    //==========================================================================
    //
    // using horizontal axis => x components of Leave/Enter
    //
    //==========================================================================
    for(unit_t j=outline.lower.y;j<=outline.upper.y;++j)
    {
        const Junction::List &JL = junctions.Horz(j);
        if(JL.size>0)
        {
            //==================================================================
            // There are some bubbles on the line
            //==================================================================
            
            //------------------------------------------------------------------
            //
            // in bulk
            //
            //------------------------------------------------------------------
            for(unit_t i=bulk_imin;i<=bulk_imax;++i)
            {
                if(B[j][i]<0)
                {
                    
                    //----------------------------------------------------------
                    // detect topology
                    //----------------------------------------------------------
                    const Real      Xi     = X[i];
                    int             status = BUBBLE_NONE;
                    const Junction *before = 0;
                    const Junction *after  = 0;
                    Real            phi    = 0;
                    Real            psi    = 0;
                    if(B[j][i-1]>=0)
                    {
                        status |= BUBBLE_BEFORE;
                        before  = JL.before(Xi);
                        if(!before)
                            throw exception("@y=%g: no junction before x=%g", Y[j], Xi);
                        psi = Xi - before->value;
                        assert(psi>=0);
                    }
                    
                    if(B[j][i+1]>=0)
                    {
                        status |= BUBBLE_AFTER;
                        after   = JL.after(Xi);
                        if(!after)
                            throw exception("@y=%g: no junction after x=%g", Y[j], Xi);
                        phi = after->value - Xi;
                        assert(phi>=0);
                    }
                    
                    //----------------------------------------------------------
                    // handle topology
                    //----------------------------------------------------------
                    switch(status)
                    {
                        case BUBBLE_BEFORE:
                        {
                            const Real Pp   = P[j][i+1];
                            const Real Pb   = before->pressure;
                            Leave[j][i-1].x = Pp - (two_delta.x) * (Pp - Pb)/(delta.x+psi);
                        }
                            break;
                            
                        case BUBBLE_AFTER:
                        {
                            const Real Pm   = P[j][i-1];
                            const Real Pb   = after->pressure;
                            Enter[j][i+1].x = Pm + (two_delta.x) * (Pb - Pm) / (delta.x + phi);
                        }
                            break;
                            
                        case BUBBLE_BOTH:
                        {
                            const Real  fac = two_delta.x/(phi+psi);
                            Enter[j][i+1].x = after->pressure  * fac;
                            Leave[j][i-1].x = before->pressure * fac;
                            
                        }
                            break;
                            
                        default:
                            break;
                    }
                }
            }
            
            
            //------------------------------------------------------------------
            //
            // on Horizonal side(s)
            //
            //------------------------------------------------------------------
        }
        
    }
    
}

void Workspace:: pressurize_vert()
{
    //==========================================================================
    //
    // using Vertical axis => y components of Leave/Enter
    //
    //==========================================================================
    for(unit_t i=outline.lower.x; i <= outline.upper.x; ++i )
    {
        const Junction::List &JL = junctions.Vert(i);
        if(JL.size)
        {
            //==================================================================
            // There are some bubbles on the line
            //==================================================================
            
            //------------------------------------------------------------------
            //
            // in bulk
            //
            //------------------------------------------------------------------
            for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
            {
                if(B[j][i]<0)
                {
                    
                    //----------------------------------------------------------
                    // detect topology
                    //----------------------------------------------------------
                    const Real      Yj     = Y[j];
                    int             status = BUBBLE_NONE;
                    const Junction *before = 0;
                    const Junction *after  = 0;
                    Real            phi    = 0;
                    Real            psi    = 0;
                    if(B[j-1][i]>=0)
                    {
                        status |= BUBBLE_BEFORE;
                        before  = JL.before(Yj);
                        if(!before)
                            throw exception("@x=%g: no junction before y=%g", X[i], Yj);
                        psi = Yj - before->value;
                    }
                    
                    if(B[j+1][i]>=0)
                    {
                        status |= BUBBLE_AFTER;
                        after   = JL.after(Yj);
                        if(!after)
                            throw exception("@x=%g: no junction after y=%g", X[i], Yj);
                        phi = after->value - Yj;
                    }
                    
                    //----------------------------------------------------------
                    // handle topology
                    //----------------------------------------------------------
                    switch(status)
                    {
                        case BUBBLE_BEFORE:
                        {
                            const Real Pp   = P[j+1][i];
                            const Real Pb   = before->pressure;
                            Leave[j-1][i].y = Pp - (two_delta.y) * (Pp - Pb)/(delta.y+psi);
                        }
                            break;
                            
                        case BUBBLE_AFTER:
                        {
                            const Real Pm   = P[j-1][i];
                            const Real Pb   = after->pressure;
                            Enter[j+1][i].y = Pm + (two_delta.y) * (Pb - Pm)/(delta.y+phi);
                        }
                            break;
                            
                        case BUBBLE_BOTH:
                        {
                            const Real fac = (two_delta.y)/(phi+psi);
                            Enter[j+1][i].y = after->pressure  * fac;
                            Leave[j-1][i].y = before->pressure * fac;
                            
                        }
                            break;
                            
                        default:
                            break;
                    }
                    
                    
                }
            }
            
            //------------------------------------------------------------------
            //
            // on Vertical side(s)
            //
            //------------------------------------------------------------------
        }
    }
    
}
#endif



static inline
bool __has_prev_inside_points( const Junction *curr )
{
    assert(curr);
    const Junction *prev = curr->prev;
    return prev && curr->inside && ( !prev->inside || prev->upper<=curr->lower);
}

static inline
bool __has_next_inside_points( const Junction *curr)
{
    assert(curr);
    const Junction *next = curr->next;
    return next && curr->inside && ( !next->inside || next->lower>=curr->upper);
}

void Workspace:: pressurize_horz()
{
    //==========================================================================
    //
    // using horizontal axis => x components of Leave/Enter
    //
    //==========================================================================
    for(unit_t j=outline.lower.y;j<=outline.upper.y;++j)
    {
        const Junction::List &JL = junctions.Horz(j);
        if(JL.size>1)
        {
            bool            in_bubble = true;
            const Junction *J = JL.head; assert(J);
            const Junction *K = J->next; assert(K);
            
            //------------------------------------------------------------------
            //
            // first outside junctions
            //
            //------------------------------------------------------------------
            if(__has_next_inside_points(J))
            {
                const unit_t i  = J->lower;
                const unit_t ip = J->upper;
                if(i>=bulk_imin)
                {
                    const Real   Xi  = X[i];
                    const Real   phi = J->value - Xi;
                    const Real   Pb  = J->pressure;
                    const Real   Pm  = P[j][i-1];
                    Enter[j][ip].x   = Pm + two_delta.x * (Pb - Pm)/(delta.x+phi);
                }
                
            }
            
            while(K)
            {
                
                if( !in_bubble )
                {
                    //----------------------------------------------------------
                    //
                    //----------------------------------------------------------
                    if(J->inside)
                    {
                        
                    }
                }
                
                in_bubble = !in_bubble;
                J=K;
                K=K->next;
            }
            
            assert(!in_bubble);
            
            //------------------------------------------------------------------
            //
            // last outside junctions
            //
            //------------------------------------------------------------------
            K=JL.tail;
            if( __has_prev_inside_points(K) )
            {
                
                const unit_t i  = K->upper;
                const unit_t im = K->lower;
                if(i<=bulk_imax)
                {
                    const Real Xi  = X[i];
                    const Real psi = Xi - K->value;
                    const Real Pb  = K->pressure;
                    const Real Pp  = P[j][i+1];
                    Leave[j][im].x = Pp - two_delta.x * (Pp - Pb ) / (delta.x+psi);
                }
                
            }
            
            
        }
        
    }
    
}

void Workspace:: pressurize_vert()
{
    
}


void Workspace:: pressurize_contours()
{
    Enter.ldz();
    Leave.ldz();
    
    pressurize_horz();
    pressurize_vert();
    
}



