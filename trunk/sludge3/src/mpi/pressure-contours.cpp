#include "workspace.hpp"

#define BUBBLE_NONE    (0)
#define BUBBLE_BEFORE  (1)
#define BUBBLE_AFTER   (2)
#define BUBBLE_BOTH   ( BUBBLE_BEFORE | BUBBLE_AFTER )

void Workspace:: pressurize_contours()
{
    Enter.ldz();
    Leave.ldz();
    
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
                    }
                    
                    if(B[j][i+1]>=0)
                    {
                        status |= BUBBLE_AFTER;
                        after   = JL.after(Xi);
                        if(!after)
                            throw exception("@y=%g: no junction after x=%g", Y[j], Xi);
                        phi = after->value - Xi;
                    }
                    
                    //----------------------------------------------------------
                    // handle topology
                    //----------------------------------------------------------
                    switch(status)
                    {
                        case BUBBLE_BEFORE:
                        {
                            const Real Pp = P[j][i+1];
                            const Real Pb = before->pressure;
                            Leave[j][i-1].x = Pp - (2*delta.x) * (Pp - Pb)/(delta.x+psi);
                        }
                            break;
                            
                        case BUBBLE_AFTER:
                        {
                            const Real Pm   = P[j][i-1];
                            const Real Pb   = after->pressure;
                            Enter[j][i+1].x = Pm + (2*delta.x) * (Pb - Pm) / (delta.x + phi);
                        }
                            break;
                            
                        case BUBBLE_BOTH:
                        {
                            const Real fac = (2*delta.x)/(phi+psi);
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
            // on Horizonal sideside(s)
            //
            //------------------------------------------------------------------
        }
        
    }
    
    
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
            
        }
    }
    
    
}


#if 0
void Workspace:: pressurize_contours()
{
    Enter.ldz();
    Leave.ldz();
    
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
            // Let us study where we are
            //==================================================================
            bool            in_bubble = true;
            const Junction *J         = JL.head;
            const Junction *K         = J->next;
            while(K)
            {
                if(in_bubble)
                {
                    if(J->owner!=K->owner)
                        throw exception("Bubble in Bubble!");
                    
                    //----------------------------------------------------------
                    // left side of the bubble
                    //----------------------------------------------------------
                    if(J->inside)
                    {
                        if( (K->inside && K->lower>=J->upper) || (!K->inside) )
                        {
                            const unit_t i   = J->lower; assert(B[j][i]<0);
                            const Real   Xi  = X[i];
                            const Real   phi = J->value - Xi; assert(phi>=0);
                            const Real   Pb  = J->pressure;
                            
                            if(i>=bulk_imin)
                            {
                                const Real Pm  = P[j][i-1];
                                Enter[j][J->upper].x = Pm + (2*delta.x)*(Pb-Pm) / (delta.x+phi);
                            }
                            else
                            {
                                P[j][i] = J->pressure;
                            }
                            // else meaningless for order 1: 0 grad on side => P side = P->bubble
                        }
                    }
                    
                    //----------------------------------------------------------
                    // right side of the bubble
                    //----------------------------------------------------------
                    if(K->inside)
                    {
                        if( (J->inside && J->upper<=K->lower) || (!J->inside) )
                        {
                            const unit_t i   = K->upper;
                            if(B[j][i]>=0)
                            {
                                std::cerr << "bubble right side error at x=" << X[i] <<", y=" << Y[j] << std::endl;
                                abort();
                            }
                            const Real   Xi  = X[i];
                            const Real   psi = Xi - K->value; assert(psi>=0);
                            const Real   Pb  = K->pressure;
                            if(i<=bulk_imax)
                            {
                                const Real Pp = P[j][i+1];
                                Leave[j][K->lower].x = Pp - (2*delta.x)*(Pp-Pb)/(delta.x+psi);
                            }
                            // else ?
                        }
                    }
                    
                    
                }
                in_bubble = !in_bubble;
                J=K;
                K=K->next;
            }
        }
    }
    
    return;
    
    //==========================================================================
    //
    // using vertical axis => y components of Leave/Enter
    //
    //==========================================================================
    for(unit_t i=outline.lower.x;i<=outline.upper.x;++i)
    {
        const Junction::List &JL = junctions.Vert(i);
        
        if(JL.size>0)
        {
            bool in_bubble = true;
            const Junction *J = JL.head;
            const Junction *K = J->next;
            while(K)
            {
                if(in_bubble)
                {
                    if(J->owner!=K->owner)
                        throw exception("Bubble in Bubble!");
                    
                    if(J->inside||K->inside)
                    {
                        const unit_t ini = J->inside ? J->upper : B.lower.y;
                        const unit_t end = K->inside ? K->lower : B.upper.y;
                        if(end>=ini)
                        {
                            Enter[ini][i].y = J->pressure;
                            Leave[end][i].y = K->pressure;
                        }
                    }
                    
                }
                in_bubble = !in_bubble;
                J=K;
                K=K->next;
            }
        }
    }
    
}
#endif


